/** \file l1d.cxx
 * \ingroup l1d3
 * \brief Lagrangian 1-Dimensional Code -- Main Program.
 *
 * This is the main module for the Lagrangian quasi-one-dimensional 
 * flow solver.  The program considers the dynamics of several gas
 * slugs, each consisting of a number of adjacent Lagrangian cells.
 *
 * The gas-dynamics is based on a control-mass formulation which tracks
 * interface positions of the gas cells.  The pressure and velocity 
 * at each interface is computed via a local Riemann problem.
 * Time-stepping is performed by a second-order explicit (predictor-
 * corrector) scheme.  
 * 
 * Piston dynamics is coupled to the gas dynamics by using the piston 
 * velocity as a boundary condition to the gas slugs and then using 
 * the computed gas pressure as the forcing term in the piston 
 * dynamics.  
 * 
 * Diaphragms are implemented as flags which allow changes to the 
 * boundary conditions on the gas slugs.
 * 
 * Running the code...
 *
 * The code reads initialization data from a parameter file and
 * an initial solution from an input-solution-file.  
 * This solution may have been prepared by l_prep or by a previous run.
 * The solution is advanced in time until a specified number of time
 * steps are taken or until a specified simulation time is reached.
 * Solutions are written to an output-solution-file periodically and 
 * upon termination.
 *
 * \version 1.0  - 02-Dec-91, basic code skeleton
 * \version 17.0 - 04-Jun-00, Adaptivity added.
 *
 * \author PA Jacobs
 * \author David Buttsworth  (USQ) 
 * \author Tony Gardner      (UQ and DLR)
 * \author Vincent Wheatley  (UQ, Caltech, Adelaide and back to UQ).
 * \author Rowan Gollan      (UQ)
 */

//-----------------------------------------------------------------

#include <vector>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
extern "C" {
#   include "../../../lib/nm/source/roberts.h"
}
#include "l_kernel.hh"
#include "l1d.hh"
#include "l_diaph.hh"
#include "l_piston.hh"
#include "l_tube.hh"
#include "l_slug.hh"
#include "l_adapt.hh"
#include "l_bc.hh"
#include "l_io.hh"
#include "l_rivp.hh"
#include "l_tstep.hh"
#include "l_bc.hh"

//-----------------------------------------------------------------

int main(int argc, char **argv)
{
    int js, jp, jd;                    /* slug, piston, diaphragm index */
    std::vector<GasSlug> A;            /* several gas slugs        */
    std::vector<PistonData> Pist;      /* room for several pistons */
    std::vector<DiaphragmData> Diaph;  /* diaphragms            */
    DiaphragmData *dp;
    int step, print_count;             /* global iteration count     */
    int adjust_end_cell_count;
    double filter_start_time;          /* start end-cell adjustment after this time */
    int halt_now;                      /* flag for premature halt    */
    double tplot;                      /* time to write next soln    */
    double thistory;                   /* time to write next sample  */
    int cfl_count;                     /* check CFL occasionally     */
    double cfl_max = 0.0;              /* current CFL maximum        */
    double cfl_tiny;                   /* smallest cfl so far        */
    double time_tiny;                  /* time at which it occurred  */
    time_t start, now;                 /* wall-clock timer           */
    int adaptive_count;                /* adapt cells occasionally   */
    int newly_adapted;
    double max_piston_V[10];
    int max_piston_V_past[10];
    char base_file_name[32];
    FILE *infile,                    /* beginning flow state         */
        *outfile,                    /* computed solution            */
        *hisfile1,                   /* single cell history          */
        *hisfile2;                   /* x-location  history          */
    int left_slug, right_slug, end_id;
    int left_slug_end, right_slug_end;
    double pressure, left_p, right_p, end_dx;
    int bc, piston_id, other_slug;
    int echo_input;
    int prepare_only;
    int i, command_line_error;
    int attempt_number, step_failed, bad_cells;
    char msg_string[256];

    printf("\n------------------------------");
    printf("\nLagrangian 1D flow simulation.");
    printf("\n------------------------------");
    printf("\n\n");

    // Initialise and set defaults.
    cfl_tiny = 1.0; /* Smallest CFL so far    */
    time_tiny = 1.0e6;
    strcpy(base_file_name, "default");
    echo_input = 0; /* by default, don't echo input */
    print_count = 50; /* Print status occasionally */
    cfl_count = 5;  /* Check CFL occasionally    */
    adaptive_count = 5; /* Adapt cells occasionally  */
    adjust_end_cell_count = 50; /* Relax the properties in the end cells */
    filter_start_time = -1.0; /* But don't do it unless instructed to do so. */
    prepare_only = 0; // By default, run a simulation.

    // Decode command line arguments.
    command_line_error = 0;
    if (argc < 2) {
	command_line_error = 1;
	goto usage;
    }
    i = 1;
    while (i < argc) {
        // Process the next command-line argument.
        if (strcmp(argv[i], "-f") == 0) {
            /* Set the base file name */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            strcpy(base_file_name, argv[i]);
            i++;
            printf("Setting base_file_name = %s\n", base_file_name);
        } else if (strcmp(argv[i], "-echo") == 0) {
            echo_input = 1;
            i++;
            printf("Will echo input parameters...\n");
        } else if (strcmp(argv[i], "-prep") == 0) {
	    prepare_only = 1;
            i++;
            printf("Will prepare initial solution files only.\n");
        } else if (strcmp(argv[i], "-print_count") == 0) {
            /* Set the interval between writing log text. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%d", &print_count);
            i++;
            printf("Setting print_count = %d\n", print_count);
        } else if (strcmp(argv[i], "-cfl_count") == 0) {
            /* Set the interval between checking cfl limit. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%d", &cfl_count);
            i++;
            printf("Setting cfl_count = %d\n", cfl_count);
        } else if (strcmp(argv[i], "-adaptive_count") == 0) {
            /* Set the interval between adapting cells. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%d", &adaptive_count);
            i++;
            printf("Setting adaptive_count = %d\n", adaptive_count);
        } else if (strcmp(argv[i], "-filter_start_time") == 0) {
            /* Set the time at which the end-cell adjustment can happen. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%lf", &filter_start_time);
            i++;
            printf("Setting filter_start_time = %g\n", filter_start_time);
        } else if (strcmp(argv[i], "-help") == 0) {
            command_line_error = 1;
            printf("Printing usage message...\n\n");
            goto usage;
        } else {
            printf("Unknown option: %s\n", argv[i]);
            command_line_error = 1;
            goto usage;
        }
    } // end while

    // If the command line arguments are incorrect,
    // write some advice to the console then finish.
    usage:
    if (command_line_error == 1) {
        printf("Command-line options: (defaults are shown in parentheses)\n");
        printf("-f <base_file_name>       (default)\n");
        printf("-print_count <n>          (%d)\n", print_count);
        printf("-cfl_count <n>            (%d)\n", cfl_count);
        printf("-adaptive_count <n>       (%d)\n", adaptive_count);
        printf("-filter_start_time <t>    (%f) negative == None\n", filter_start_time);
        printf("-echo                     (no echo)\n");
        printf("-prep                     (run)\n");
        printf("-help          (print this message)\n");
        exit(1);
    } // end if command_line_error

    // Build the file names.
    string pname = string(base_file_name) + ".Lp";
    string iname = string(base_file_name) + ".L0";
    string aname = string(base_file_name) + ".La";
    string oname = string(base_file_name) + ".Ls";
    string hname1 = string(base_file_name) + ".hc";
    string hname2 = string(base_file_name) + ".hx";
    string efname = string(base_file_name) + ".event";
    string dname = string(base_file_name) + ".dump";

    // Pick up the problem description data from the parameter (INI) file.
    SimulationData SD = SimulationData(pname, echo_input);
    SD.sim_time = 0.0;  /* Global simulation time */
    TubeModel tube = TubeModel(pname, echo_input);
    for (jp = 0; jp < SD.npiston; ++jp) {
        Pist.push_back(PistonData(jp, SD.dt_init, pname, echo_input));
        Pist[jp].sim_time = 0.0;
    }
    for (jd = 0; jd < SD.ndiaphragm; ++jd) {
        Diaph.push_back(DiaphragmData(jd, pname, echo_input));
        Diaph[jd].sim_time = 0.0;
    }
    SD.hncell = 0;
    for (js = 0; js < SD.nslug; ++js) {
        A.push_back(GasSlug(js, SD, pname, echo_input));
        SD.hncell += A[js].hncell;
        A[js].sim_time = 0.0;
    }
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();

    if ( prepare_only ) {
	double x[11000]; 
	// Not that we have made a guess for the required size; 
	// we should use a std::vector but x is used in the call 
	// to distribute_points().
	double beta1, beta2;
	printf("Prepare initial solution files only.\n");
	printf("Set up gas slugs.\n");
	for (js = 0; js < SD.nslug; ++js) {
	    /*
	     * Distribute the gas cells along the gas slug.
	     * Modified to suit sm_3d stretching functions.
	     */
	    if ( A[js].cluster_to_end_1 == 1 ) {
		beta1 = A[js].cluster_strength;
	    } else {
		beta1 = 0.0;
	    }
	    if ( A[js].cluster_to_end_2 == 1 ) {
		beta2 = A[js].cluster_strength;
	    } else {
		beta2 = 0.0;
	    }
	    // FIX-ME: one day, it will be good to use a vector for x
	    // such that its size can be adjested as needed.
	    distribute_points(A[js].xbegin, A[js].xend, A[js].nnx, x, beta1, beta2); 
	    int i = 0;
	    for ( int ix = A[js].ixmin-1; ix <= A[js].ixmax; ++ix) {
		A[js].Cell[ix].x = x[i];
		++i;
	    }
	    A[js].compute_areas(&tube);
	    A[js].fill_data();
	    A[js].encode_conserved();
	    /* 
	     * The following line should fill in all of the 
	     * extra variables.
	     */
	    if ( A[js].decode_conserved() != 0 ) {
		printf( "Failure decoding conserved quantities for slug[%d].\n", js );
		exit( -1 );
	    }
	}   /* end for js... */
	tube.write_area(aname);
	tube.write_dump_file(dname);
	printf("Write starting solution file.\n");
	if ((outfile = fopen(iname.c_str(), "w")) == NULL) {
	    printf("\nCould not open %s; BAILING OUT\n", oname.c_str());
	    return FAILURE;
	}
	for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].write_state(outfile);
	for (jd = 0; jd < SD.ndiaphragm; ++jd) Diaph[jd].write_state(outfile);
	for (js = 0; js < SD.nslug; ++js) A[js].write_state(outfile);
	if (outfile != NULL) fclose(outfile);
	return SUCCESS;
    }

    // Pick up the initial data that was previously generated.
    if ((infile = fopen(iname.c_str(), "r")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", iname.c_str());
        exit(1);
    }
    for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].read_state(infile);
    for (jd = 0; jd < SD.ndiaphragm; ++jd) Diaph[jd].read_state(infile);
    for (js = 0; js < SD.nslug; ++js) {
        A[js].read_state(infile);
        A[js].compute_areas(&tube);
        A[js].encode_conserved();
        // Fill in all of the extra variables.
        if ( A[js].decode_conserved() != 0 ) {
	    printf( "Failure decoding conserved quantities for slug[%d].\n", js );
	    printf( "This occured just after reading starting solution.\n" );
	    exit(-1);
        }
	A[js].set_chemistry_timestep(-1.0);
	A[js].set_thermal_timestep(-1.0);
    } // end for js...
    if ( infile != NULL ) fclose(infile);

    SD.sim_time = A[0].sim_time; // Pick up the old time.
    tplot = SD.sim_time + SD.get_dt_plot();
    thistory = SD.sim_time + SD.get_dt_history();
    SD.dt_global = SD.dt_init;
    for (js = 0; js < SD.nslug; ++js) {
        A[js].dt = SD.dt_global;
        A[js].cfl_target = SD.CFL;
    }

    // Open output files to catch flow and history data.
    if ((outfile = fopen(oname.c_str(), "w")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", oname.c_str());
        return FAILURE;
    }
    if ((hisfile1 = fopen(hname1.c_str(), "w")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", hname1.c_str());
        return FAILURE;
    }
    if ((hisfile2 = fopen(hname2.c_str(), "w")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", hname2.c_str());
        return FAILURE;
    }

    newly_adapted = 0;
    step = 0;   /* Global Iteration Count    */
    halt_now = 0;   /* no other reason to stop   */
    for ( jp=0; jp < SD.npiston; ++jp ) {
        max_piston_V[jp] = Pist[jp].V;
        max_piston_V_past[jp] = 0;
    }

    //-----------------------------------
    printf("\nBeginning MAIN LOOP...\n");
    //-----------------------------------
    start = time(NULL); // start of wallclock timing
    fflush(stdout);
    sprintf( msg_string, "\nStart time stepping... " );
    strcat( msg_string, "\n" );
    log_event( efname.c_str(), msg_string );

    while (SD.sim_time <= SD.max_time && step <= SD.max_step && halt_now == 0) {

	// --------------------------------
	// 0. Adjust the cell distribution.
	// --------------------------------
        if (step > 0 && (step / adaptive_count) * adaptive_count == step) {
            for (js = 0; js < SD.nslug; ++js) {
                if (A[js].adaptive > ADAPT_NONE) {
                    L_adapt_cells(&(A[js]));
                    newly_adapted = 1;
                }
	    } // end for js...
        }

	// --------------------------------
        // 1. Set the size of the time step.
        // --------------------------------
        if (step == 0 ||
            (step / cfl_count) * cfl_count == step || newly_adapted == 1) {
            // Check the CFL number and determine an allowable time-step.
            for (js = 0; js < SD.nslug; ++js) A[js].check_cfl();
            SD.dt_allow = A[0].dt_allow;
            cfl_max = A[0].cfl_max;
            if (SD.nslug > 1) {
                for (js = 1; js < SD.nslug; ++js) {
                    if (A[js].dt_allow < SD.dt_allow) SD.dt_allow = A[js].dt_allow;
                    if (A[js].cfl_max > cfl_max) cfl_max = A[js].cfl_max;
                } // end for
            } // end if nslug > 1

            if (step < 3) {
                // Use initial time step.
                SD.dt_global = SD.dt_init;
            } else {
                // Adjust time step.
                if (SD.dt_allow > SD.dt_global) {
                    // cautious increase
                    SD.dt_global += 0.5 * (SD.dt_allow - SD.dt_global);
                } else {
                    // rapid decrease
                    SD.dt_global = SD.dt_allow;
                }
            } // end if step < 3

            // Propagate the time step to each piston and slug.
            for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].dt = SD.dt_global;
            for (js = 0; js < SD.nslug; ++js) A[js].dt = SD.dt_global;
        } // end if step...

        for (js = 0; js < SD.nslug; ++js) {
            if ((step > 5) && (A[js].cfl_max < cfl_tiny)) {
                // Record the worst CFL.
                cfl_tiny = A[js].cfl_max;
                time_tiny = SD.sim_time;
            }
            if (A[js].cfl_max > 100.0) {
                // If the CFL has gone crazy, bail out safely.
                printf("\n\n---------------------------\n");
                printf("Slug %d, CFL = %e: Bailing out!\n", js,
                       A[js].cfl_max);
                halt_now = 1;
                goto BailOut;
            }
        } // end for js... 

        // --------------------------------
	// 2. Prepare to take the time step.
	// --------------------------------
        if ((step / print_count) * print_count == step) {
            // Print the current time-stepping, piston and peak pressures.  
            now = time(NULL);
            printf("\n");
            printf("-------- Wall-Clock seconds = %d --------\n",
                   (int)(now - start) );
            print_simulation_status(stdout, NULL, step, SD, A, Diaph, Pist,
				    cfl_max, cfl_tiny, time_tiny);
            fflush(stdout); // make it appear now
        }

        if ((step / adjust_end_cell_count) * adjust_end_cell_count == step &&
	    filter_start_time >= 0.0 && SD.sim_time > filter_start_time) {
	    // Occasionally, relax the gas properties in the end cells
	    // towards the values in their immediate neighbours in the
	    // same slug.
	    // This will, hopefully, eliminate the glitches seen in the
	    // beginning test gas in Ben's expansion tube simulations.
            for (js = 0; js < SD.nslug; ++js) {
		A[js].adjust_end_cells();
	    }
	}

	// ---------------------------
	// 3. Deal with the diaphragms.
	// ---------------------------
	// This amounts to rupturing (so-far unruptured) diaphragms
	// when the specified pressure difference is exceeded.
        for (jd = 0; jd < SD.ndiaphragm; ++jd) {
	    dp = &( Diaph[jd] );

            if ( dp->is_burst == 0 ) {
                // Check only unruptured diaphragms for burst conditions. 

                // Pressure of the left-side of the diaphragm.
                js = dp->left_slug_id;
                if (js >= 0) {
                    end_id = dp->left_slug_end_id;
                    end_dx = dp->left_slug_dx;
                    left_p = A[js].end_pressure(end_id, end_dx);
                } else {
                    left_p = 0.0;
                }

                // Pressure of the right-side of the diaphragm.
                js = dp->right_slug_id;
                if (js >= 0) {
                    end_id = dp->right_slug_end_id;
                    end_dx = dp->right_slug_dx;
                    right_p = A[js].end_pressure(end_id, end_dx);
                } else {
                    right_p = 0.0;
                }

                // Check for excess pressure as the trigger for rupture.
                if ((fabs(left_p - right_p) >= dp->P_burst)
                    && dp->trigger_time < 0.0) {
                    dp->trigger_time = SD.sim_time;
                    sprintf( msg_string,
                             "\nEvent: diaphragm[%d] trigger at t= %e\n",
                             jd, SD.sim_time );
                    log_event( efname.c_str(), msg_string );
		    print_simulation_status(NULL, efname.c_str(), step, SD, A, Diaph, Pist,
					    cfl_max, cfl_tiny, time_tiny);
                }

                // Wait the hold time before rupturing diaphragm.
                if ( dp->trigger_time >= 0.0 && 
		     (SD.sim_time - dp->trigger_time) > dp->hold_period ) {
                    dp->is_burst = 1;
                    sprintf( msg_string, "\nEvent: diaphragm[%d] rupture at t= %e\n",
                             jd, SD.sim_time );
                    log_event( efname.c_str(), msg_string );
		    print_simulation_status(NULL, efname.c_str(), step, SD, A, Diaph, Pist,
					    cfl_max, cfl_tiny, time_tiny);
                } // end if dp->trigger_time >= 0.0 &&...
            } else {
		// For ruptured diaphragms, check to see if we should blend
		// the gas-slug data after a period of time.

		if ( dp->blend_delay > 0.0 && !(dp->already_blended) &&
		     SD.sim_time > (dp->trigger_time + dp->hold_period + 
				    dp->blend_delay) ) {
		    L_blend_slug_ends( &(A[dp->left_slug_id]), dp->left_slug_end_id,
				       &(A[dp->right_slug_id]), dp->right_slug_end_id,
				       dp->blend_dx );
		    A[dp->left_slug_id].compute_areas(&tube);
		    A[dp->left_slug_id].encode_conserved();
		    A[dp->right_slug_id].compute_areas(&tube);
		    A[dp->right_slug_id].encode_conserved();
		    dp->already_blended = 1;
                    sprintf( msg_string,
                             "\nEvent: blend slugs [%d] and [%d] at t= %e\n",
                             dp->left_slug_id, dp->right_slug_id, SD.sim_time );
                    log_event( efname.c_str(), msg_string );
		    print_simulation_status( NULL, efname.c_str(), step, SD, A, Diaph, Pist,
				  cfl_max, cfl_tiny, time_tiny );
		}
            } // end if dp->is_burst == 0 ... else ...
        } // end for jd... 

	// ----------------------
	// 4. Update the dynamics 
	// ----------------------
        for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].record_state(); 
        for (js = 0; js < SD.nslug; ++js) {
            A[js].record_state();
        }

        attempt_number = 0;
        do {
            ++attempt_number;
            step_failed = 0;

	    // 4a. Predictor Stage.
	    // Boundary conditions for the gas slugs.
            for (js = 0; js < SD.nslug; ++js) {
                // Deal with left end condition first. 
                bc = A[js].left_bc_type;
                if (bc == FREE_END) {
                    L_bc_left_free(&(A[js]));
                } else if (bc == SOLID_BOUNDARY) {
                    // The appropriate velocity (possibly zero) was set earlier.
                    L_bc_left_velocity(&(A[js]), A[js].left_ustar);
                } else if (bc == PISTON) {
                    piston_id = A[js].left_piston_id;
                    L_bc_left_velocity(&(A[js]), Pist[piston_id].V);
                } else if (bc == SLUG) {
                    // This is double handling but is OK.
                    other_slug = A[js].left_slug_id;
                    L_exchange_bc_data(&(A[other_slug]), &(A[js]));
                } else if (bc == SLUG_DIAPHRAGM) {
                    jd = A[js].left_diaphragm_id;
                    if (Diaph[jd].is_burst == 0) {
                        L_bc_left_reflect(&(A[js]));
                    } else {
                        // This is double handling but is OK.
                        other_slug = A[js].left_slug_id;
                        L_exchange_bc_data(&(A[other_slug]), &(A[js]));
                    }
                } else {
                    printf("L1d: slug[%d] invalid left BC %d\n", js, bc);
                } // end if bc... 

		// Deal with right end condition second. 
                bc = A[js].right_bc_type;
                if (bc == FREE_END) {
                    L_bc_right_free(&(A[js]));
                } else if (bc == SOLID_BOUNDARY) {
                    // The appropriate velocity (possibly zero) was set earlier.
                    L_bc_right_velocity(&(A[js]), A[js].right_ustar);
                } else if (bc == PISTON) {
                    piston_id = A[js].right_piston_id;
                    L_bc_right_velocity(&(A[js]), Pist[piston_id].V);
                } else if (bc == SLUG) {
                    // This is double handling but is OK.
                    other_slug = A[js].right_slug_id;
                    L_exchange_bc_data(&(A[js]), &(A[other_slug]));
                } else if (bc == SLUG_DIAPHRAGM) {
                    jd = A[js].right_diaphragm_id;
                    if (Diaph[jd].is_burst == 0) {
                        L_bc_right_reflect(&(A[js]));
                    } else {
                        // This is double handling but is OK.
                        other_slug = A[js].right_slug_id;
                        L_exchange_bc_data(&(A[js]), &(A[other_slug]));
                    }
                } else {
                    printf("L1d: slug[%d] invalid right BC %d\n", js, bc);
                } // end if bc, for right end...
            } // end for js...

            // 4b. Gas-dynamic predictor step
            for (js = 0; js < SD.nslug; ++js) {
                A[js].apply_rivp();
                A[js].source_vector();
                A[js].axial_heat_flux(SD.k);
                A[js].time_derivatives(0);
                A[js].predictor_step();
                A[js].compute_areas(&tube);
                if ( A[js].decode_conserved() != 0 ) {
		    printf( "decode_conserved() failed at predictor step\n" );
		    printf( "   for slug[%d], attempt=%d, time-step=%d\n", 
			    js, attempt_number, step );
		}
            } // end for

            // 4c. Boundary conditions for the pistons.
            for (jp = 0; jp < SD.npiston; ++jp) {
                left_slug = Pist[jp].left_slug_id;
                left_slug_end = Pist[jp].left_slug_end_id;
                if (left_slug >= 0) {
                    // Apply gas-slug pressure to left (back) face.
                    if (left_slug_end == RIGHT) {
                        pressure = A[left_slug].right_pstar;
                    } else {
                        pressure = A[left_slug].left_pstar;
                    }
                } else {
                    // There is no gas slug against this face.
                    pressure = 0.0;
                }
                Pist[jp].Pb = pressure;

                right_slug = Pist[jp].right_slug_id;
                right_slug_end = Pist[jp].right_slug_end_id;
                if (right_slug >= 0) {
                    // Apply gas-slug pressure to right (front) face.
                    if (right_slug_end == RIGHT) {
                        pressure = A[right_slug].right_pstar;
                    } else {
                        pressure = A[right_slug].left_pstar;
                    }
                } else {
                    // There is no gas slug against this face.
                    pressure = 0.0;
                }
                Pist[jp].Pf = pressure;
            } // end for jp...

            // 4d. Piston-dynamic predictor step
            for (jp = 0; jp < SD.npiston; ++jp) {
                Pist[jp].time_derivatives(0, SD.sim_time);
                Pist[jp].predictor_step();
            }

	    // -------------------
	    // 4e. Corrector Stage 
	    // -------------------
            if (SD.Torder == 2) {
		// Apply boundary conditions to the gas slugs.
		// Leave the diaphragms as they were for the predictor
		// stage. 
                for (js = 0; js < SD.nslug; ++js) {
                    // Deal with left end condition first. 
                    bc = A[js].left_bc_type;
                    if (bc == FREE_END) {
                        L_bc_left_free(&(A[js]));
                    } else if (bc == SOLID_BOUNDARY) {
                        // The appropriate velocity (possibly zero) was set earlier.
                        L_bc_left_velocity(&(A[js]), A[js].left_ustar);
                    } else if (bc == PISTON) {
                        piston_id = A[js].left_piston_id;
                        L_bc_left_velocity(&(A[js]), Pist[piston_id].V);
                    } else if (bc == SLUG) {
                        // This is double handling but is OK.
                        other_slug = A[js].left_slug_id;
                        L_exchange_bc_data(&(A[other_slug]), &(A[js]));
                    } else if (bc == SLUG_DIAPHRAGM) {
                        jd = A[js].left_diaphragm_id;
                        if (Diaph[jd].is_burst == 0) {
                            L_bc_left_reflect(&(A[js]));
                        } else {
                            // This is double handling but is OK.
                            other_slug = A[js].left_slug_id;
                            L_exchange_bc_data(&(A[other_slug]), &(A[js]));
                        }
                    } else {
                        printf("L1d: slug[%d] invalid left BC %d\n", js, bc);
                    } // end if bc...

                    // Deal with right end condition second. 
                    bc = A[js].right_bc_type;
                    if (bc == FREE_END) {
                        L_bc_right_free(&(A[js]));
                    } else if (bc == SOLID_BOUNDARY) {
                        // The appropriate velocity (possibly zero) was set earlier.
                        L_bc_right_velocity(&(A[js]), A[js].right_ustar);
                    } else if (bc == PISTON) {
                        piston_id = A[js].right_piston_id;
                        L_bc_right_velocity(&(A[js]), Pist[piston_id].V);
                    } else if (bc == SLUG) {
                        // This is double handling but is OK.
                        other_slug = A[js].right_slug_id;
                        L_exchange_bc_data(&(A[js]), &(A[other_slug]));
                    } else if (bc == SLUG_DIAPHRAGM) {
                        jd = A[js].right_diaphragm_id;
                        if (Diaph[jd].is_burst == 0) {
                            L_bc_right_reflect(&(A[js]));
                        } else {
                            // This is double handling but is OK.
                            other_slug = A[js].right_slug_id;
                            L_exchange_bc_data(&(A[js]), &(A[other_slug]));
                        }
                    } else {
                        printf("L1d: slug[%d] invalid right BC %d\n", js, bc);
                    } // end if bc, for right end...
                } // end for js...

                // 4f. Gas-dynamic corrector step
                for (js = 0; js < SD.nslug; ++js) {
                    A[js].apply_rivp();
                    A[js].source_vector();
                    A[js].axial_heat_flux(SD.k);
                    A[js].time_derivatives(1);
                    A[js].corrector_step();
                    A[js].compute_areas(&tube);
		    if ( A[js].decode_conserved() != 0 ) {
			printf( "decode_conserved() failed at corrector step\n" );
			printf( "   for slug[%d], attempt=%d, time-step=%d\n", 
				js, attempt_number, step );
		    }
                }

                // 4g. Boundary conditions for the pistons.
                for (jp = 0; jp < SD.npiston; ++jp) {
                    left_slug = Pist[jp].left_slug_id;
                    left_slug_end = Pist[jp].left_slug_end_id;
                    if (left_slug >= 0) {
                        // Apply gas-slug pressure to left (back) face.
                        if (left_slug_end == RIGHT) {
                            pressure = A[left_slug].right_pstar;
                        } else {
                            pressure = A[left_slug].left_pstar;
                        }
                    } else {
                        // There is no gas slug against this face.
                        pressure = 0.0;
                    }
                    Pist[jp].Pb = pressure;

                    right_slug = Pist[jp].right_slug_id;
                    right_slug_end = Pist[jp].right_slug_end_id;
                    if (right_slug >= 0) {
                        // Apply gas-slug pressure to right (front) face.
                        if (right_slug_end == RIGHT) {
                            pressure = A[right_slug].right_pstar;
                        } else {
                            pressure = A[right_slug].left_pstar;
                        }
                    } else {
                        // There is no gas slug against this face.
                        pressure = 0.0;
                    }
                    Pist[jp].Pf = pressure;
                } // end for jp...

                // 4h. Piston-dynamic corrector step
                for (jp = 0; jp < SD.npiston; ++jp) {
                    Pist[jp].time_derivatives(1, SD.sim_time);
                    Pist[jp].corrector_step();
                }
            } // end of corrector step.

	    // Now, with fixed volume and total energy in each cell,
	    // update the chemical species using Rowan's finite-rate
	    // chemistry module.
	    if ( SD.fr_chem == 1 ) {
		for (js = 0; js < SD.nslug; ++js) {
		    A[js].chemical_increment(SD.dt_global);
		}
	    }

            // Check for bad cells. 
            bad_cells = 0;
            for (js = 0; js < SD.nslug; ++js) {
                bad_cells += A[js].check_cells(js);
            }
            step_failed = (bad_cells > 0);

	    // If the attempt has failed, reduce the 
	    // time step for the next attempt at the
	    // current step.
            if ( step_failed == 1 ) {
                // The reduction factor is somewhat arbitrary.
                SD.dt_global *= 0.2;
                printf("WARNING: time-step attempt %d failed at t=%g\n",
                       attempt_number, SD.sim_time);
                printf("Reducing to dt=%g\n", SD.dt_global);

                // Propagate the time step to each piston and slug.
                for (jp = 0; jp < SD.npiston; ++jp)
                    Pist[jp].dt = SD.dt_global;
                for (js = 0; js < SD.nslug; ++js)
                    A[js].dt = SD.dt_global;

                // Restore the state which existed before the attempt.
                for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].restore_state();

                for (js = 0; js < SD.nslug; ++js) {
                    A[js].restore_state();
                    A[js].compute_areas(&tube);
                    if ( A[js].decode_conserved() != 0 ) {
			printf("decode_conserved() failed while trying\n");
			printf("   to restore original state for slug[%d]\n", js);
			printf("   at time-step=%d\n", step );
		    }
                } // end for js...
            } // end if step_failed == 1...
        } while ( attempt_number < 3 && step_failed == 1 );

        if ( step_failed == 1 ) {
            printf("\n");
            printf("Time step failed at t = %g.\n", SD.sim_time);
            break; // leave the main time-stepping loop.
        }

        // At this point in the main time-stepping loop, 
	// the time step should have been successful so
	// we update the time-step number, etc.
        ++step;
        SD.sim_time += SD.dt_global;
        for ( jp = 0; jp < SD.npiston; ++jp ) {
            Pist[jp].sim_time = SD.sim_time;
        }
        for ( jd = 0; jd < SD.ndiaphragm; ++jd ) {
            Diaph[jd].sim_time = SD.sim_time;
        }
        for ( js = 0; js < SD.nslug; ++js ) {
            A[js].sim_time = SD.sim_time;
        }

        // --------------------------------------------
	// 5. Intermediate solution data to be written? 
	// --------------------------------------------
	// 5a. Full flow along tube, diaphragm and piston states
        if ( SD.sim_time >= tplot ) {
            tplot += SD.get_dt_plot();
            for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].write_state(outfile);
            for (jd = 0; jd < SD.ndiaphragm; ++jd) Diaph[jd].write_state(outfile);
            for (js = 0; js < SD.nslug; ++js) A[js].write_state(outfile);
        }
	// 5b. Selected history points.
        if ( SD.sim_time >= thistory ) {
            thistory += SD.get_dt_history();
            fprintf(hisfile1, "%e %d %d %d # sim_time, hncell, nsp, nmodes\n", 
		    SD.sim_time, SD.hncell, nsp, nmodes);
            for (js = 0; js < SD.nslug; ++js)
                L_write_cell_history(&(A[js]), hisfile1);
            fprintf(hisfile2, "%e %d %d %d  # sim_time, hnloc, nsp, nmodes\n", 
		    SD.sim_time, SD.hnloc, nsp, nmodes);
            for (js = 0; js < SD.hnloc; ++js)
                L_write_x_history(SD.hxloc[js], A, SD.nslug, hisfile2);
        }

	// -------------------------
        // 6. Piston special events.
	// -------------------------
        for ( jp = 0; jp < SD.npiston; ++jp ) {
            if ( Pist[jp].V_old * Pist[jp].V < 0.0 ) {
                sprintf( msg_string,
                         "\nEvent: piston[%d] reversal at t= %e x= %e\n",
                         jp, SD.sim_time, Pist[jp].x );
                log_event( efname.c_str(), msg_string );
		print_simulation_status(NULL, efname.c_str(), step, SD, A, Diaph, Pist,
					cfl_max, cfl_tiny, time_tiny);
                // After reversal, look for new maximum.
                max_piston_V[jp] = 0.0;
                max_piston_V_past[jp] = 0;
            } // end if piston reversal

            if ( Pist[jp].brakes_on_old == 0 && Pist[jp].brakes_on == 1 ) {
                sprintf( msg_string,
                         "\nEvent: piston[%d] brakes on at t= %e x= %e\n",
                         jp, SD.sim_time, Pist[jp].x );
                log_event( efname.c_str(), msg_string );
		print_simulation_status( NULL, efname.c_str(), step, SD, A, Diaph, Pist,
			      cfl_max, cfl_tiny, time_tiny );
            } // end if piston brake applied

            if ( Pist[jp].is_restrain_old == 1 && Pist[jp].is_restrain == 0 ) {
                sprintf( msg_string,
                         "\nEvent: piston[%d] released at t= %e x= %e\n",
                         jp, SD.sim_time, Pist[jp].x );
                log_event( efname.c_str(), msg_string );
		print_simulation_status( NULL, efname.c_str(), step, SD, A, Diaph, Pist,
			      cfl_max, cfl_tiny, time_tiny );
            } // end if piston released

            if ( fabs(Pist[jp].V) > max_piston_V[jp] ) {
                max_piston_V[jp] = fabs(Pist[jp].V);
            }
            if ( fabs(Pist[jp].V) < (max_piston_V[jp] - 1.0e-6) && max_piston_V_past[jp] == 0 ) {
                max_piston_V_past[jp] = 1;
                sprintf( msg_string, 
                         "\nEvent: piston[%d] peak speed at t= %e V= %e\n",
                         jp, SD.sim_time, Pist[jp].V );
                log_event( efname.c_str(), msg_string );
		print_simulation_status( NULL, efname.c_str(), step, SD, A, Diaph, Pist,
			      cfl_max, cfl_tiny, time_tiny );
            } // end if piston max speed
        } // end for jp...

	// 7. Loop termination criteria:
	//    (1) reaching a maximum simulation time
	//    (2) reaching a maximum number of steps
	//    (3) finding that the "halt" file exists in the working directory
	//        This provides a semi-interactive way to terminate the 
	//        simulation and save the data.
	//    Criteria 1 & 2 are tested at the top of the loop.
        if (access("l1d.halt", F_OK) == 0) {
            halt_now = 1;
            printf("Simulation stopped: Halt file exists.\n");
        }

        newly_adapted = 0;
    } // end while (i.e. time step / iteration)

    // ----------
    // Conclusion
    // ----------
    BailOut:
    log_event( efname.c_str(), (const char *)"\nEnd time stepping.\n" );
    print_simulation_status(NULL, efname.c_str(), step, SD, A, Diaph, Pist,
			    cfl_max, cfl_tiny, time_tiny );
    printf("\nTotal number of steps = %d\n", step);

    for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].write_state(outfile);
    for (jd = 0; jd < SD.ndiaphragm; ++jd) Diaph[jd].write_state(outfile);
    for (js = 0; js < SD.nslug; ++js) A[js].write_state(outfile);
    for (js = 0; js < SD.nslug; ++js) {
        delete &(A[js]);
    }
    delete gmodel;

    if (outfile != NULL) fclose(outfile);
    if (hisfile1 != NULL) fclose(hisfile1);
    if (hisfile2 != NULL) fclose(hisfile2);

    return SUCCESS;
} // end function main()
