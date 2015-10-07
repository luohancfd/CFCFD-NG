/** \file main.cxx
 * \ingroup eilmer3
 * \brief Multiple-Block Complete Navier-Stokes -- main program.
 *
 * \author PA Jacobs for the core stepper, geometry description and
 *         the boring house-keeping bits such as macros in the viscous derivatives.
 *         Various postgrad students and post docs for the interesting bits.
 *         Andrew McGhee and Richard Goozee -- MPI parallel computing.
 *         Chris Craddock and Ian Johnston -- first crack at thermochemistry, 3D
 *         Paul Petrie, Ian Johnston and Andrew Pastrello -- flux calculators, grid management
 *         Rowan Gollan  -- thermochemistry, radiation
 *         Andrew Denman -- 3D viscous effects and LES turbulence
 *         Joseph Tang -- hierarchical grids
 *         Brendan O'Flaherty -- piston, thermochemistry, real gas models
 *         Jan-Pieter Nap -- help with the k-omega turbulence
 *         Wilson Chan -- more work on the turbulence models
 *         Dan Potter -- radiation, thermochemistry
 *
 *         ...and quite a few people have spent lots of effort "kicking the tyres",
 *         suffering through buggy versions of new features and generally helping
 *         by making new and interesting applications.
 *         Mike Wendt -- lots of enthusiasm in the early days when the code was
 *                       in its most buggy state
 *         Rainer Kirchhartz -- patience with the turbulent scramjet models
 *         Katsu Tanimizu -- the 3D scramjet work
 *
 * \version 20-Mar-2006 -- New-generation using some C++.
 * \version    Jan-2008 -- split out material to cns_kernel.cxx, much like Elmer2.
 * \version 25-Nov-2008 -- ported last of the big pieces (3D viscous terms)
 */

/*-----------------------------------------------------------------*/
// Intel MPI requires mpi.h included BEFORE stdio.h
#ifdef _MPI
#   include <mpi.h>
#endif
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
extern "C" {
#include <popt.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "../../../lib/util/source/time_to_go.h"
}
#include "../../../lib/util/source/lua_service.hh"
#include "../../wallcon/source/e3conn.hh"
#include "cell.hh"
#include "block.hh"
#include "kernel.hh"
#include "init.hh"
#include "bc.hh"
#include "bc_extrapolate_out.hh"
#include "lua_service_for_e3.hh"
#include "bc_menter_correction.hh"
#include "exch2d.hh"
#include "exch_mapped_cell_shmem.hh"
#include "visc.hh"
#include "visc3D.hh"
#ifdef _MPI
#   include "exch_mpi.hh"
#   include "exch_mapped_cell_mpi.hh"
#endif
#include "piston.hh"
#include "implicit.hh"
#include "conj-ht-interface.hh"
#ifdef GPU_CHEM
#    include "gpu-chem-update.hh"
#endif
#include "main.hh"

//-----------------------------------------------------------------
// Global data
//
bool history_just_written, output_just_written, av_output_just_written;
bool write_at_step_has_been_done = false;
int program_return_flag = 0;
size_t output_counter = 0; // counts the number of flow-solutions written
bool zip_files = true; // flag to indicate if flow and grid files are to be gzipped
bool with_heat_flux_files = false; // flag to indicate that we want heat-flux files
bool with_surface_files = false; // flag to indicate that we want surface files
bool master;
int max_wall_clock = 0; // seconds
time_t start, now; // wall-clock timer

lua_State *L; // for the uder-defined procedures

#ifdef GPU_CHEM
gpu_chem *gchem = 0;
#endif

//-----------------------------------------------------------------

void usage(poptContext optCon, int exitcode, char *error, char *addl) {
    // Print program usage message, copied from popt example.
    if ( master ) {
	poptPrintUsage(optCon, stderr, 0);
	if (error) fprintf(stderr, "%s: %s0", error, addl);
    }
    exit(exitcode);
}

//-----------------------------------------------------------------
// Begin here...
//
int main(int argc, char **argv)
{
    global_data &G = *get_global_data_ptr();
    G.verbosity_level = 1; // emit messages useful for monitoring a long-running process
    bool do_run_simulation = false;
    size_t start_tindx = 0;
    int run_status = SUCCESS;
    char c, job_name[132], text_buf[132], *dotpy;
    char log_file_name[132], mpimap_file_name[132];

#   ifdef _MPI
    MPI_Status status;
    int node_name_len;
    char node_name[MPI_MAX_PROCESSOR_NAME];
#   endif

    poptContext optCon;   /* context for parsing command-line options */
    struct poptOption optionsTable[] = {
	{ "job", 'f', POPT_ARG_STRING, NULL, 'f',
	  "job_name is typically the same as root_file_name", 
	  "<job_name>" },
	{ "run", 'r', POPT_ARG_NONE, NULL, 'r',
	  "run the main simulation time-stepper", 
	  NULL },
	{ "tindx", 't', POPT_ARG_STRING, NULL, 't',
	  "start with this set of flow data", 
	  "<int>" },
	{ "zip-files", 'z', POPT_ARG_NONE, NULL, 'z',
	  "use gzipped flow and grid files", 
	  NULL },
	{ "no-zip-files", 'a', POPT_ARG_NONE, NULL, 'a',
	  "use ASCII (not gzipped) flow and grid files", 
	  NULL },
	{ "heat-flux-files", 'q', POPT_ARG_NONE, NULL, 'q',
	  "write heat-flux files", 
	  NULL },
	{ "surface-files", 's', POPT_ARG_NONE, NULL, 's',
	  "write surface files", 
	  NULL },
	{ "verbosity", 'v', POPT_ARG_STRING, NULL, 'v',
	  "set verbosity level for messages", 
	  "<int>" },
	{ "max-wall-clock", 'w', POPT_ARG_STRING, NULL, 'w',
	  "maximum wall-clock time in seconds", 
	  "<seconds>" },
	{ "mpimap", 'm', POPT_ARG_STRING, NULL, 'm',
	  "use this specific MPI map of blocks to rank", 
	  "<mpimap_file>" },
	POPT_AUTOHELP
	POPT_TABLEEND
    };

    // Initialization
    //
    // Each process writes to a separate log file,
    // but all processes may write to stdout as well.
    program_return_flag = SUCCESS; // We'll start out optimistically :)
#   ifdef _MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &(G.num_mpi_proc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(G.my_mpi_rank));
    UNUSED_VARIABLE(status);
    G.mpi_parallel = 1;
    master = (G.my_mpi_rank == 0);
    sprintf(log_file_name, "e3mpi.%04d.log", G.my_mpi_rank);
#   else
    G.mpi_parallel = 0;
    G.num_mpi_proc = 0;
    G.my_mpi_rank = 0;
    master = true;
#   ifdef E3RAD
    sprintf(log_file_name, "e3rad.log");
#   else
    sprintf(log_file_name, "e3shared.log");
#   endif
#   endif
    if ((G.logfile = fopen(log_file_name, "w")) == NULL) {
        printf( "\nCould not open %s; BAILING OUT\n", log_file_name );
        exit(FILE_ERROR);
    }

    // Configuration from command-line.
    optCon = poptGetContext(NULL, argc, (const char**)argv, optionsTable, 0);
    if (argc < 2) {
	poptPrintUsage(optCon, stderr, 0);
	goto Quit;
    }
    strcpy(job_name, "");
    strcpy(mpimap_file_name, "");
    // Do options processing as per the popt example.
    while ((c = poptGetNextOpt(optCon)) >= 0) {
	// printf("received option %c\n", c);
	switch (c) {
	case 'f':
	    strcpy(job_name, poptGetOptArg(optCon));
	    dotpy = strstr(job_name, ".py");
	    if ( dotpy ) { *dotpy = '\0'; } // cut off the .py extension
	    break;
	case 'r':
	    do_run_simulation = true;
	    break;
	case 'z':
	    zip_files = true;
	    break;
	case 'a':
	    zip_files = false;
	    break;
	case 'q':
	    with_heat_flux_files = true;
	    break;
	case 's':
	    with_surface_files = true;
	    break;
	case 't':
	    start_tindx = static_cast<size_t>(atoi(poptGetOptArg(optCon)));
	    break;
	case 'v':
	    G.verbosity_level = static_cast<size_t>(atoi(poptGetOptArg(optCon)));
	    break;
	case 'w':
	    strcpy( text_buf, poptGetOptArg(optCon) );
	    sscanf( text_buf, "%d", &max_wall_clock );
	    break;
	case 'm':
	    strcpy(mpimap_file_name, poptGetOptArg(optCon));
	    break;
	}
    }

    if ( G.verbosity_level >= 1 ) {
	// Now, that we have sorted who we are and have some config established,
	// write an introductory message.
#       ifdef _MPI
	if ( master ) {
	    printf("e3main: C++,MPI version.\n");
	}
	MPI_Get_processor_name(node_name, &node_name_len);
	printf("e3main: process %d of %d active on node %s\n", 
	       G.my_mpi_rank, G.num_mpi_proc, node_name );
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#       elif E3RAD
	printf("e3rad: C++ version for radiation transport calculations.\n");
	printf("------------------------------------------------------------------\n"
	       "WARNING: This executable only computes the radiative source\n"
	       "         terms and heat fluxes - no time integration is performed!\n"
	       "------------------------------------------------------------------\n");
#       ifdef _OPENMP
	printf("OpenMP version using %d thread(s).\n", omp_get_max_threads());
#       endif
#       else
	printf("e3main: C++,shared-memory version.\n");
#       ifdef _OPENMP
	printf("OpenMP version using %d thread(s).\n", omp_get_max_threads());
#       error "OpenMP version not functional..."    
#       endif
#       endif
	if ( master ) {
	    cout << "Source code revision string: " << get_revision_string() << endl;
	    printf( "job_name=%s\n", job_name );
	    if ( max_wall_clock > 0 ) {
		printf( "max_wall_clock=%d seconds\n", max_wall_clock );
	    } else {
		printf( "Run will not be limited by wall clock.\n" );
	    }
	}
    } // end if ( G.verbosity_level

    G.base_file_name = string(job_name);
    // Read the static configuration parameters from the INI file
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( SUCCESS != read_config_parameters(G.base_file_name+".config", master, start_tindx) ) goto Quit; 
    // Read the time-step control parameters from a separate file.
    // These will also be read at the start of each time step.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( SUCCESS != read_control_parameters(G.base_file_name+".control", master, true) ) goto Quit;
    // At this point, we know the number of blocks in the calculation.
    // Depending on whether we are running all blocks in the one process
    // or we are running a subset of blocks in this process, talking to
    // the other processes via MPI, we need to decide what blocks belong
    // to the current process.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( SUCCESS != assign_blocks_to_mpi_rank(mpimap_file_name, master) ) goto Quit;

    // Real work is delegated to functions that know what to do...
    //
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    if ( do_run_simulation ) {
	if ( G.verbosity_level >= 1 && master ) {
	    printf("Run simulation...\n");
	}
	try {
	    // The simulation proper.
	    run_status = prepare_to_integrate(start_tindx);
	    if (run_status != SUCCESS) goto Quit;
#           ifndef E3RAD
	    if ( G.sequence_blocks ) {
		run_status = integrate_blocks_in_sequence();
		if (run_status != SUCCESS) goto Quit;
	    } else {
		run_status = integrate_in_time(-1.0);
		if (run_status != SUCCESS) goto Quit;
	    }
#           else
	    run_status = radiation_calculation();
	    if (run_status != SUCCESS) goto Quit;
#           endif
	    finalize_simulation();
	} catch (std::runtime_error &e) {
	    cout << "Well, this *is* embarrassing! The calculation failed." << endl;
	    cout << e.what() << endl;
	    run_status = FAILURE;
	}
    } else {
	printf( "NOTHING DONE -- because you didn't ask...\n" );
    }

    // Finalization.
    //
    Quit: /* nop */;
    fclose(G.logfile);
    eilmer_finalize();
    if ( G.verbosity_level >= 1 ) {
#       ifdef _MPI
	printf("e3main: end of process %d.\n", G.my_mpi_rank);
#       elif E3RAD
	printf("e3rad: done.\n");
#       else
	printf("e3main: done.\n");
#       endif
    }
#   ifdef _MPI
    if (run_status != SUCCESS) 
	MPI_Abort(MPI_COMM_WORLD, program_return_flag);
    else
	MPI_Finalize();
#   endif
    return program_return_flag;
} // end main

//------------------------------------------------------------------------------

void ensure_directory_is_present( string pathname )
{
    string commandstring;
    int cmd_return_status;
    if ( master ) {
	if ( access(pathname.c_str(), F_OK) != 0 ) {
	    commandstring = "mkdir " + pathname;
	    cmd_return_status = system(commandstring.c_str());
	    if ( cmd_return_status != 0 ) {
		cerr << "Problem with system cmd: " << commandstring << endl;
		cerr << "Quitting simulation." << endl;
		exit(FILE_ERROR);
	    }
	}
    }
#   ifdef _MPI
    // We really do need this barrier.
    // Even if we are not the master node, we have to wait for the master
    // to check and (possibly) make the directory.
    MPI_Barrier(MPI_COMM_WORLD);
#   endif
    return;
} // end ensure_directory_is_present()


void do_system_cmd( string commandstring )
{
    int cmd_return_status;
    cmd_return_status = system(commandstring.c_str());
    if ( cmd_return_status != 0 ) {
	cerr << "Problem with system cmd: " << commandstring << endl;
	cerr << "Quitting simulation." << endl;
	exit(FILE_ERROR);
    }
    return;
} // end do_system_cmd()

//-----------------------------------------------------------------------------

int prepare_to_integrate(size_t start_tindx)
{
    global_data &G = *get_global_data_ptr();
    string filename, commandstring, jbstring, tindxstring, ifstring;
    char jbcstr[10], tindxcstr[10], ifcstr[10];

    // Allocate and initialise memory.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    for ( Block *bdp : G.my_blocks ) {
        if ( bdp->array_alloc(G.dimensions) != SUCCESS ) exit( MEMORY_ERROR );
	bdp->bind_interfaces_to_cells(G.dimensions);
    }
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
    if ( allocate_send_and_receive_buffers() != 0 ) exit( MEMORY_ERROR );
#   endif
    // Read block grid and flow data; write history-file headers.
    // Note that the global simulation time is set by the last flow data read.
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
#   endif
    sprintf( tindxcstr, "t%04d", static_cast<int>(start_tindx));
    tindxstring = tindxcstr;
    for ( Block *bdp : G.my_blocks ) {
        if ( G.verbosity_level >= 2 ) printf( "----------------------------------\n" );
	sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) );
	jbstring = jbcstr;
	// Read grid from the specified tindx files.
	if ( G.moving_grid ) {
	    filename = "grid/"+tindxstring+"/"+G.base_file_name+".grid"+jbstring+"."+tindxstring;
	} else {
	    // Read grid from the tindx=0 files, always.
	    filename = "grid/t0000/"+G.base_file_name+".grid"+jbstring+".t0000";
	}
        if (bdp->read_grid(filename, G.dimensions, zip_files) != SUCCESS) {
	    return FAILURE;
	}
	// Read flow data from the specified tindx files.
	filename = "flow/"+tindxstring+"/"+G.base_file_name+".flow"+jbstring+"."+tindxstring;
        if (bdp->read_solution(filename, &(G.sim_time), G.dimensions, zip_files) != SUCCESS) {
	    return FAILURE;
	}
	if ( G.BGK == 2 ) {
	    filename = "flow/"+tindxstring+"/"+G.base_file_name+".BGK"+jbstring+"."+tindxstring;
	    if ( access(filename.c_str(), F_OK) != 0 ) {
		// previous BGK velocity distributions do exist, try to read them in
		if (bdp->read_BGK(filename, &(G.sim_time), G.dimensions, zip_files) != SUCCESS) {
		    return FAILURE;
		}
	    }
	} else if ( G.BGK == 1) {
	    // assume equilibrium velocity distribution, generate from conserved props
	    if (bdp->initialise_BGK_equilibrium() != SUCCESS) {
		return FAILURE;
	    }
	    filename = "flow/"+tindxstring+"/"+G.base_file_name+".BGK"+jbstring+"."+tindxstring;
	    bdp->write_BGK(filename, G.sim_time, G.dimensions, zip_files);
	}
    } // end for *bdp

    // History file header is only written for a fresh start.
    ensure_directory_is_present("hist"); // includes Barrier

    for ( Block *bdp : G.my_blocks ) {
        if ( G.verbosity_level >= 2 ) printf( "----------------------------------\n" );
	sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) );
	jbstring = jbcstr;
	filename = "hist/" + G.base_file_name + ".hist"+jbstring;
	if ( access(filename.c_str(), F_OK) != 0 ) {
	    // History file does not yet exist; write header.
	    bdp->write_history(filename, G.sim_time, true);
	}
	for ( int iface: bdp->transient_profile_faces ) {
	    filename = "./" + G.base_file_name + ".blk"+jbstring+"."+get_face_name(iface)+".profile";
	    if ( access(filename.c_str(), F_OK) != 0 ) {
		// History file does not yet exist; write header.
		bdp->write_profile(filename, iface, G.sim_time, true);
	    }
	}
	// Read in heat-flux vectors if present
	for ( int iface = NORTH; iface <= ((G.dimensions == 3)? BOTTOM : WEST); ++iface ) {
	    sprintf( ifcstr, ".s%04d", iface );
	    ifstring = ifcstr;
	    filename = "heat/"+tindxstring+"/"+G.base_file_name+".heat"+jbstring+ifstring+"."+tindxstring;
	    bdp->bcp[iface]->read_surface_heat_flux(filename, G.dimensions, zip_files);
	}
	// Read in surface vectors if present
	for ( int iface = NORTH; iface <= ((G.dimensions == 3)? BOTTOM : WEST); ++iface ) {
	    sprintf( ifcstr, ".s%04d", iface );
	    ifstring = ifcstr;
	    filename = "surf/"+tindxstring+"/"+G.base_file_name+".surf"+jbstring+ifstring+"."+tindxstring;
	    bdp->bcp[iface]->read_surface_data(filename, G.dimensions, zip_files);
	}
    } // end for *bdp
    output_counter = start_tindx;
    filename = G.base_file_name; filename += ".times";
    if ( master ) {
	if ((G.timestampfile = fopen(filename.c_str(), "a")) == NULL) {
	    cerr << "Could not open " << filename << "; BAILING OUT" << endl;
	    return FILE_ERROR;
	}
	// The zeroth entry in the timestampfile has been written already by e3prep.py.
	// fprintf( G.timestampfile, "# count sim_time dt_global\n" );
	// fprintf( G.timestampfile, "%04d %e %e\n", output_counter, G.sim_time, G.dt_global );
	// fflush( G.timestampfile );
    }

#   ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#   endif

    // Prepare data within the primary finite-volume cells.
    // This includes both geometric data and flow state data.
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    bool first = true;
    // NOTE: these are temporary vectors that are only used
    // in the case that the NORTH boundary is a CONJUGATE_HT_BC.
    // We need to store this information from each MPI rank
    // so that we can assemble it globally and pass it on
    // to some objects for initialisation.
    // We split open the Vector3 objects because it's easier
    // to pass around arrays of doubles using MPI functions.
    vector<double> local_cell_x_pos; // storage of near wall x-pos of local block
    vector<double> local_cell_y_pos; // storage of near wall y-pos of local block
    vector<double> local_iface_x_pos; // storage of interface x-pos of local block
    vector<double> local_iface_y_pos; // storage of interface y-pos of local block
    vector<double> local_iface_nx; // storage of interface n.x of local block
    vector<double> local_iface_ny;
    for ( Block *bdp : G.my_blocks ) {
	bdp->compute_primary_cell_geometric_data(G.dimensions, 0);
	bdp->compute_distance_to_nearest_wall_for_all_cells(G.dimensions, 0);
	bdp->compute_secondary_cell_geometric_data(G.dimensions, 0);
	bdp->set_base_qdot(G, 0);  // this need be done only once per block
	bdp->identify_reaction_zones(G, 0);
	bdp->identify_turbulent_zones(G, 0);
	for ( FV_Cell *cp: bdp->active_cells ) {
	    cp->encode_conserved(0, 0, bdp->omegaz, with_k_omega);
	    // Even though the following call appears redundant at this point,
	    // fills in some gas properties such as Prandtl number that is
	    // needed for both the cfd_check and the BLomax turbulence model.
	    cp->decode_conserved(0, 0, bdp->omegaz, with_k_omega);
	}
	if ( G.conjugate_ht_active ) {
	    if ( bdp->bcp[NORTH]->type_code == CONJUGATE_HT ) {
		// Need to gather some geometric information
		int j = bdp->jmax;
		for ( size_t i = bdp->imin; i <= bdp->imax; ++i ) {
		    FV_Cell *cell = bdp->get_cell(i, j);
		    local_cell_x_pos.push_back(cell->pos[0].x);
		    local_cell_y_pos.push_back(cell->pos[0].y);
		    FV_Interface *iface = bdp->get_ifj(i, j+1);
		    local_iface_x_pos.push_back(iface->pos.x);
		    local_iface_y_pos.push_back(iface->pos.y);
		    local_iface_nx.push_back(iface->n.x);
		    local_iface_ny.push_back(iface->n.y);
		}
	    }
	}
	// update the 'global' data
	if (first) {
	    G.L_min = bdp->L_min;
	    G.bounding_box_min = bdp->bounding_box_min;
	    G.bounding_box_max = bdp->bounding_box_max;
	    first = false;
	} else {
	    if (bdp->L_min < G.L_min) G.L_min = bdp->L_min;

	    if (bdp->bounding_box_max.x > G.bounding_box_max.x) G.bounding_box_max.x = bdp->bounding_box_max.x;
	    if (bdp->bounding_box_max.y > G.bounding_box_max.y) G.bounding_box_max.y = bdp->bounding_box_max.y;
	    if (bdp->bounding_box_max.z > G.bounding_box_max.z) G.bounding_box_max.z = bdp->bounding_box_max.z;

	    if (bdp->bounding_box_min.x < G.bounding_box_min.x) G.bounding_box_min.x = bdp->bounding_box_min.x;
	    if (bdp->bounding_box_min.y < G.bounding_box_min.y) G.bounding_box_min.y = bdp->bounding_box_min.y;
	    if (bdp->bounding_box_min.z < G.bounding_box_min.z) G.bounding_box_min.z = bdp->bounding_box_min.z;
	}
	
    } // end for *bdp
    
#   ifdef _MPI
    // Ensure that all blocks have computed geometry before continuing
    MPI_Barrier(MPI_COMM_WORLD);

    // update the true global minimum cell size
    MPI_Allreduce(MPI_IN_PLACE, &(G.L_min), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // update the global bounding box
    MPI_Allreduce(MPI_IN_PLACE, &(G.bounding_box_min.x), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &(G.bounding_box_min.y), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &(G.bounding_box_min.z), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &(G.bounding_box_max.x), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &(G.bounding_box_max.y), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &(G.bounding_box_max.z), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#   endif

    if (G.MHD && (G.divB_damping_length <= 0.0)) {
		Vector3 L_xyz = 0.5*(G.bounding_box_max - G.bounding_box_min);
		if (L_xyz.z > 0.0) {
		    G.divB_damping_length = 1.0 / (3.141592653589793 * sqrt(1/(L_xyz.x*L_xyz.x) + 1/(L_xyz.y*L_xyz.y) + 1/(L_xyz.z*L_xyz.z)));
		} else {
		    G.divB_damping_length = 1.0 / (3.141592653589793 * sqrt(1/(L_xyz.x*L_xyz.x) + 1/(L_xyz.y*L_xyz.y)));
		}
    }

    if ( G.conjugate_ht_active ) {
	// Need to gather some data for coordination of wall update
	// at global level
	vector<double> cell_xs, cell_ys, wall_xs, wall_ys, wall_nxs, wall_nys;
#       ifdef _MPI
	// For MPI version, we need to gather geometric information
	// from each rank.
	// Deal with x positions of near wall cell centres
	cell_xs.resize(G.T_gas_near_wall.size());
	double *plcxs = &local_cell_x_pos[0];
	double *pcxs = &cell_xs[0];
	int *recvcounts = &(G.recvcounts[0]);
	int *displs = &(G.displs[0]);
	MPI_Gatherv(plcxs, local_cell_x_pos.size(), MPI_DOUBLE, pcxs, recvcounts, displs,
		    MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Deal with y positions of near wall cell centres
	cell_ys.resize(G.T_gas_near_wall.size());
	double *plcys = &local_cell_y_pos[0];
	double *pcys = &cell_ys[0];
	MPI_Gatherv(plcys, local_cell_y_pos.size(), MPI_DOUBLE, pcys, recvcounts, displs,
		    MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Deal with x positions at wall (interface)
	wall_xs.resize(G.T_gas_near_wall.size());
	double *plixs = &local_iface_x_pos[0];
	double *pixs = &wall_xs[0];
	MPI_Gatherv(plixs, local_iface_x_pos.size(), MPI_DOUBLE, pixs, recvcounts, displs,
		    MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Deal with y positions at wall (interface)
	wall_ys.resize(G.T_gas_near_wall.size());
	double *pliys = &local_iface_y_pos[0];
	double *piys = &wall_ys[0];
	MPI_Gatherv(pliys, local_iface_y_pos.size(), MPI_DOUBLE, piys, recvcounts, displs,
		    MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Deal with n.x at wall (interface)
	wall_nxs.resize(G.T_gas_near_wall.size());
	double *plinxs = &local_iface_nx[0];
	double *pinxs = &wall_nxs[0];
	MPI_Gatherv(plinxs, local_iface_nx.size(), MPI_DOUBLE, pinxs, recvcounts, displs,
		    MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Deal with n.y at wall (interface)
	wall_nys.resize(G.T_gas_near_wall.size());
	double *plinys = &local_iface_ny[0];
	double *pinys = &wall_nys[0];
	MPI_Gatherv(plinys, local_iface_ny.size(), MPI_DOUBLE, pinys, recvcounts, displs,
		    MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#       else
        // In shared memory, the local storage vectors should have
	// all the values needed.
	cell_xs = local_cell_x_pos;
	cell_ys = local_cell_y_pos;
	wall_xs = local_iface_x_pos;
	wall_ys = local_iface_y_pos;
	wall_nxs = local_iface_nx;
	wall_nys = local_iface_ny;
#       endif
	if ( master ) {
	    // Also, assemble and initialise the conjugate interface object
	    vector<Vector3> gas_cell_pos, solid_cell_pos, iface_pos, iface_normal;
	    gas_cell_pos.resize(cell_xs.size());
	    iface_pos.resize(cell_xs.size());
	    iface_normal.resize(cell_xs.size());
	    for ( size_t i = 0; i < cell_xs.size(); ++i ) {
		gas_cell_pos[i].x = cell_xs[i];
		gas_cell_pos[i].y = cell_ys[i];
		iface_pos[i].x = wall_xs[i];
		iface_pos[i].y = wall_ys[i];
		iface_normal[i].x = wall_nxs[i];
		iface_normal[i].y = wall_nys[i];
	    }
	    // Ask wall model to give its cell centre positions for
	    // the adjacent cells to the coupled interface
	    get_near_wall_solid_cell_pos(*(G.wm), solid_cell_pos);
	    if ( solid_cell_pos.size() != gas_cell_pos.size() ) {
		cout << "Error matching cells in gas and solid domains\n";
		cout << "for the conjugate heat transfer model.\n";
		cout << "number of solid cells: " << solid_cell_pos.size() << endl;
		cout << "number of gas cells: " << gas_cell_pos.size() << endl;
		cout << "Bailing out!\n";
		exit(FAILURE);
	    }
	    // Initialise interface object
	    G.conj_ht_iface = new Conjugate_HT_Interface(gas_cell_pos, solid_cell_pos, 
							 iface_pos, iface_normal);
	    // Ask wall model to write a solution to file if this is the initial solution
	    if ( start_tindx == 0 )
		write_solution(*(G.wm), G.sim_time, start_tindx);
	    
	}
    }

    // Exchange boundary cell geometry information so that we can
    // next calculate secondary-cell geometries.
    // [TODO] Check that we don't break the usage for integrate_blocks_in_sequence.
#   ifdef _MPI
    // Before we try to exchange data, everyone's internal data should be up-to-date.
    MPI_Barrier( MPI_COMM_WORLD );
    // Now, it's safe to do the exchange for full-face connections.
    mpi_exchange_boundary_data(COPY_CELL_LENGTHS, 0);
    copy_mapped_cell_data_via_mpi(COPY_CELL_LENGTHS, 0);
#   else
    for ( Block *bdp : G.my_blocks ) {
        exchange_shared_boundary_data(bdp->id, COPY_CELL_LENGTHS, 0);
    }
    copy_mapped_cell_data_via_shmem(COPY_CELL_LENGTHS, 0);
#   endif

    // Start up the Lua interpreter and load the external file
    // containing the user-defined functions -- if appropriate.
    if ( G.udf_file.length() > 0 ) {
	L = luaL_newstate();
	luaL_openlibs(L); // load the standard libraries
	// Set up some global variables that might be handy in 
	// the Lua environment.
	// We give the number of blocks on a per rank basis
	lua_pushinteger(L, static_cast<int>(G.my_blocks.size()));
	lua_setglobal(L, "nblks");
	int nsp = get_gas_model_ptr()->get_number_of_species();
	int nmodes = get_gas_model_ptr()->get_number_of_modes();
	lua_pushinteger(L, nsp);
	lua_setglobal(L, "nsp");
	lua_pushinteger(L, nmodes);
	lua_setglobal(L, "nmodes");
	lua_newtable(L); // table for block info at TOS
	for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    lua_pushinteger(L, static_cast<int>(jb));
	    lua_newtable(L);
	    // Push entries into new anoymous table
	    lua_pushstring(L, "id");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->id));
	    lua_settable(L, -3);
	    lua_pushstring(L, "nicells");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->nni));
	    lua_settable(L, -3);
	    lua_pushstring(L, "imin");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->imin));
	    lua_settable(L, -3);
	    lua_pushstring(L, "imax");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->imax));
	    lua_settable(L, -3);
	    lua_pushstring(L, "njcells");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->nnj));
	    lua_settable(L, -3);
	    lua_pushstring(L, "jmin");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->jmin));
	    lua_settable(L, -3);
	    lua_pushstring(L, "jmax");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->jmax));
	    lua_settable(L, -3);
	    lua_pushstring(L, "nkcells");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->nnk));
	    lua_settable(L, -3);
	    lua_pushstring(L, "kmin");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->kmin));
	    lua_settable(L, -3);
	    lua_pushstring(L, "kmax");
	    lua_pushinteger(L, static_cast<int>(G.my_blocks[jb]->imax));
	    lua_settable(L, -3);
	    // Label anonymous table with index
	    lua_settable(L, -3);
	}
	lua_setglobal(L, "blks");
	// Register functions so that they are accessible 
	// from the Lua environment.
	register_luafns(L);
	// Presume that the user-defined functions are in the specified file.
	if ( luaL_loadfile(L, G.udf_file.c_str()) || lua_pcall(L, 0, 0, 0) ) {
	    handle_lua_error(L, "Could not run user file: %s", lua_tostring(L, -1));
	}
	lua_settop(L, 0); // clear the stack
    }
    // Initialise radiation transport if appropriate
    if ( G.radiation ) {
    	if ( get_radiation_transport_model_ptr()->initialise() ) {
	    cerr << "Problem with initialisation of radiation transport data\n";
	    cerr << "Exiting program." << endl;
	    exit(FAILURE);
	}
    }
#ifdef GPU_CHEM
    // Initialise the gpu chemistry update, first count number of cells in simulation
    int ncells = 0;
    for ( Block *bdp : G.my_blocks ) ncells += bdp->active_cells.size();
    gchem = init_gpu_module(ncells);
#endif
    start = time(NULL); // start of wallclock timing
    return SUCCESS;
} // end prepare_to_integrate()

//------------------------------------------------------------------------

int call_udf(double t, size_t step, std::string udf_fn_name)
{
    lua_getglobal(L, udf_fn_name.c_str());  // function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");
    lua_pushinteger(L, static_cast<int>(step)); lua_setfield(L, -2, "step");
    int number_args = 1; // table of {t, step}
    int number_results = 0; // no results returned on the stack.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user-defined function: %s\n", lua_tostring(L, -1));
    }
    lua_settop(L, 0); // clear the stack
    return SUCCESS;
} // end call_udf()

int call_udf2(double t, size_t step, size_t blk_id, std::string udf_fn_name)
{
    lua_getglobal(L, udf_fn_name.c_str());  // function to be called
    if ( udf_fn_name == "before_grid_update" && lua_isnil(L, -1) ) {
	// If the user hasn't defined before_grid_update,
	// we'll just return without doing anything
	return SUCCESS;
    }
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");
    lua_pushinteger(L, static_cast<int>(step)); lua_setfield(L, -2, "step");
    lua_pushinteger(L, static_cast<int>(blk_id)); lua_setfield(L, -2, "block_id");
    int number_args = 1; // table of {t, step}
    int number_results = 0; // no results returned on the stack.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user-defined function: %s\n", lua_tostring(L, -1));
    }
    lua_settop(L, 0); // clear the stack
    return SUCCESS;
} // end call_udf()

/// \brief Add to the components of the source vector, Q, via a Lua udf.
/// This is done (occasionally) just after the for inviscid source vector calculation.
int add_udf_source_vector_for_cell( FV_Cell *cell, size_t gtl, double t )
{
    // Call the user-defined function which returns a table
    // of source term values.
    // These are added to the inviscid source terms 
    // that were computed earlier in the time step.

    global_data &G = *get_global_data_ptr();
    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();

    lua_getglobal(L, "source_vector");  // Lua function to be called
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, t); lua_setfield(L, -2, "t");

    // Pack the interesting cell data in a table with named fields.
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, cell->pos[gtl].x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, cell->pos[gtl].y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, cell->pos[gtl].z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, cell->volume[gtl]); lua_setfield(L, -2, "vol");
    lua_pushnumber(L, cell->fs->gas->p); lua_setfield(L, -2, "p");
    lua_pushnumber(L, cell->fs->gas->rho); lua_setfield(L, -2, "rho"); 
    lua_pushnumber(L, cell->fs->vel.x); lua_setfield(L, -2, "u"); 
    lua_pushnumber(L, cell->fs->vel.y); lua_setfield(L, -2, "v");
    lua_pushnumber(L, cell->fs->vel.z); lua_setfield(L, -2, "w");
    lua_pushnumber(L, cell->fs->gas->a); lua_setfield(L, -2, "a");
    lua_pushnumber(L, cell->fs->gas->mu); lua_setfield(L, -2, "mu");
    lua_newtable(L); // A table for the temperatures
    for ( size_t i = 0; i < nmodes; ++i ) {
	lua_pushinteger(L, i);
	lua_pushnumber(L, cell->fs->gas->T[i]);
	lua_settable(L, -3);
    }
    // At this point, the table of temperatures should be TOS.
    lua_setfield(L, -2, "T");
    lua_newtable(L); // Another table for the mass fractions
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_pushnumber(L, cell->fs->gas->massf[isp]);
	lua_settable(L, -3);
    }
    // At this point, the table of mass fractions should be TOS.
    lua_setfield(L, -2, "massf");
    lua_newtable(L); // A table for the thermal conductivities.
    for ( size_t i = 0; i < nmodes; ++i ) {
	lua_pushinteger(L, i);
	lua_pushnumber(L, cell->fs->gas->k[i]);
	lua_settable(L, -3);
    }
    // At this point, the table of conductivities should be TOS.
    lua_setfield(L, -2, "k");

    // Derivatives needed for Zac Denman's verification of the turbulence model.
    // The following has been mostly lifted from FV_Cell::k_omega_time_derivatives()
    // in cell.cxx.
    double dudx, dudy, dudz;
    double dvdx, dvdy, dvdz;
    double dwdx, dwdy, dwdz;
    double dtkedx, dtkedy, dtkedz;
    double domegadx, domegady, domegadz;
    if ( G.dimensions == 2 ) {
        dudx = 0.25 * (cell->vtx[0]->dudx + cell->vtx[1]->dudx + cell->vtx[2]->dudx + cell->vtx[3]->dudx);
        dudy = 0.25 * (cell->vtx[0]->dudy + cell->vtx[1]->dudy + cell->vtx[2]->dudy + cell->vtx[3]->dudy);
	dudz = 0.0;
        dvdx = 0.25 * (cell->vtx[0]->dvdx + cell->vtx[1]->dvdx + cell->vtx[2]->dvdx + cell->vtx[3]->dvdx);
        dvdy = 0.25 * (cell->vtx[0]->dvdy + cell->vtx[1]->dvdy + cell->vtx[2]->dvdy + cell->vtx[3]->dvdy);
	dvdz = 0.0;
	dwdx = 0.0; dwdy = 0.0; dwdz = 0.0;
        dtkedx = 0.25 * (cell->vtx[0]->dtkedx + cell->vtx[1]->dtkedx + cell->vtx[2]->dtkedx + cell->vtx[3]->dtkedx);
        dtkedy = 0.25 * (cell->vtx[0]->dtkedy + cell->vtx[1]->dtkedy + cell->vtx[2]->dtkedy + cell->vtx[3]->dtkedy);
	dtkedz = 0.0;
        domegadx = 0.25 * (cell->vtx[0]->domegadx + cell->vtx[1]->domegadx + cell->vtx[2]->domegadx + cell->vtx[3]->domegadx);
        domegady = 0.25 * (cell->vtx[0]->domegady + cell->vtx[1]->domegady + cell->vtx[2]->domegady + cell->vtx[3]->domegady);
	domegadz = 0.0;
    } else { // 3D cartesian
        dudx = 0.125 * (cell->vtx[0]->dudx + cell->vtx[1]->dudx + cell->vtx[2]->dudx + cell->vtx[3]->dudx +
                        cell->vtx[4]->dudx + cell->vtx[5]->dudx + cell->vtx[6]->dudx + cell->vtx[7]->dudx);
        dudy = 0.125 * (cell->vtx[0]->dudy + cell->vtx[1]->dudy + cell->vtx[2]->dudy + cell->vtx[3]->dudy +
                        cell->vtx[4]->dudy + cell->vtx[5]->dudy + cell->vtx[6]->dudy + cell->vtx[7]->dudy);
        dudz = 0.125 * (cell->vtx[0]->dudz + cell->vtx[1]->dudz + cell->vtx[2]->dudz + cell->vtx[3]->dudz +
                        cell->vtx[4]->dudz + cell->vtx[5]->dudz + cell->vtx[6]->dudz + cell->vtx[7]->dudz);
        dvdx = 0.125 * (cell->vtx[0]->dvdx + cell->vtx[1]->dvdx + cell->vtx[2]->dvdx + cell->vtx[3]->dvdx +
                        cell->vtx[4]->dvdx + cell->vtx[5]->dvdx + cell->vtx[6]->dvdx + cell->vtx[7]->dvdx);
        dvdy = 0.125 * (cell->vtx[0]->dvdy + cell->vtx[1]->dvdy + cell->vtx[2]->dvdy + cell->vtx[3]->dvdy +
                        cell->vtx[4]->dvdy + cell->vtx[5]->dvdy + cell->vtx[6]->dvdy + cell->vtx[7]->dvdy);
        dvdz = 0.125 * (cell->vtx[0]->dvdz + cell->vtx[1]->dvdz + cell->vtx[2]->dvdz + cell->vtx[3]->dvdz +
                        cell->vtx[4]->dvdz + cell->vtx[5]->dvdz + cell->vtx[6]->dvdz + cell->vtx[7]->dvdz);
        dwdx = 0.125 * (cell->vtx[0]->dwdx + cell->vtx[1]->dwdx + cell->vtx[2]->dwdx + cell->vtx[3]->dwdx +
                        cell->vtx[4]->dwdx + cell->vtx[5]->dwdx + cell->vtx[6]->dwdx + cell->vtx[7]->dwdx);
        dwdy = 0.125 * (cell->vtx[0]->dwdy + cell->vtx[1]->dwdy + cell->vtx[2]->dwdy + cell->vtx[3]->dwdy +
                        cell->vtx[4]->dwdy + cell->vtx[5]->dwdy + cell->vtx[6]->dwdy + cell->vtx[7]->dwdy);
        dwdz = 0.125 * (cell->vtx[0]->dwdz + cell->vtx[1]->dwdz + cell->vtx[2]->dwdz + cell->vtx[3]->dwdz +
                        cell->vtx[4]->dwdz + cell->vtx[5]->dwdz + cell->vtx[6]->dwdz + cell->vtx[7]->dwdz);
        dtkedx = 0.125 * (cell->vtx[0]->dtkedx + cell->vtx[1]->dtkedx + cell->vtx[2]->dtkedx + cell->vtx[3]->dtkedx +
                          cell->vtx[4]->dtkedx + cell->vtx[5]->dtkedx + cell->vtx[6]->dtkedx + cell->vtx[7]->dtkedx);
        dtkedy = 0.125 * (cell->vtx[0]->dtkedy + cell->vtx[1]->dtkedy + cell->vtx[2]->dtkedy + cell->vtx[3]->dtkedy +
                          cell->vtx[4]->dtkedy + cell->vtx[5]->dtkedy + cell->vtx[6]->dtkedy + cell->vtx[7]->dtkedy);
        dtkedz = 0.125 * (cell->vtx[0]->dtkedz + cell->vtx[1]->dtkedz + cell->vtx[2]->dtkedz + cell->vtx[3]->dtkedz +
                          cell->vtx[4]->dtkedz + cell->vtx[5]->dtkedz + cell->vtx[6]->dtkedz + cell->vtx[7]->dtkedz);
        domegadx = 0.125 * (cell->vtx[0]->domegadx + cell->vtx[1]->domegadx + cell->vtx[2]->domegadx + cell->vtx[3]->domegadx +
                            cell->vtx[4]->domegadx + cell->vtx[5]->domegadx + cell->vtx[6]->domegadx + cell->vtx[7]->domegadx);
        domegady = 0.125 * (cell->vtx[0]->domegady + cell->vtx[1]->domegady + cell->vtx[2]->domegady + cell->vtx[3]->domegady +
                            cell->vtx[4]->domegady + cell->vtx[5]->domegady + cell->vtx[6]->domegady + cell->vtx[7]->domegady);
        domegadz = 0.125 * (cell->vtx[0]->domegadz + cell->vtx[1]->domegadz + cell->vtx[2]->domegadz + cell->vtx[3]->domegadz +
                            cell->vtx[4]->domegadz + cell->vtx[5]->domegadz + cell->vtx[6]->domegadz + cell->vtx[7]->domegadz);
    } // end if ( G.dimensions...
    lua_pushnumber(L, dudx); lua_setfield(L, -2, "dudx");
    lua_pushnumber(L, dudy); lua_setfield(L, -2, "dudy");
    lua_pushnumber(L, dudz); lua_setfield(L, -2, "dudz");
    lua_pushnumber(L, dvdx); lua_setfield(L, -2, "dvdx");
    lua_pushnumber(L, dvdy); lua_setfield(L, -2, "dvdy");
    lua_pushnumber(L, dvdz); lua_setfield(L, -2, "dvdz");
    lua_pushnumber(L, dwdx); lua_setfield(L, -2, "dwdx");
    lua_pushnumber(L, dwdy); lua_setfield(L, -2, "dwdy");
    lua_pushnumber(L, dwdz); lua_setfield(L, -2, "dwdz");
    lua_pushnumber(L, dtkedx); lua_setfield(L, -2, "dtkedx");
    lua_pushnumber(L, dtkedy); lua_setfield(L, -2, "dtkedy");
    lua_pushnumber(L, dtkedz); lua_setfield(L, -2, "dtkedz");
    lua_pushnumber(L, domegadx); lua_setfield(L, -2, "domegadx");
    lua_pushnumber(L, domegady); lua_setfield(L, -2, "domegady");
    lua_pushnumber(L, domegadz); lua_setfield(L, -2, "domegadz");
    
    // After all of this we should have ended up with the cell-data table at TOS 
    // with another table (containing t...} and function-name 
    // at pseudo-index locations -2 and -3 respectively.

    int number_args = 2; // table of {t...}, cell-data-table
    int number_results = 1; // one table of results returned on the stack.
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "error running user source_terms function: %s\n",
			 lua_tostring(L, -1));
    }
    // Assume that there is a single table at the TOS
    // cout << "Get mass, momentum and energy source terms." << endl;
    lua_getfield(L, -1, "mass"); cell->Q->mass += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_x"); cell->Q->momentum.x += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_y"); cell->Q->momentum.y += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "momentum_z"); cell->Q->momentum.z += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "total_energy"); cell->Q->total_energy += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "romega"); cell->Q->omega += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "rtke"); cell->Q->tke += lua_tonumber(L, -1); lua_pop(L, 1);
    lua_getfield(L, -1, "radiation"); cell->Q_rE_rad += lua_tonumber(L, -1); lua_pop(L, 1);
    cell->Q->total_energy += cell->Q_rE_rad;

    // cout << "Table of sources for mass fractions" << endl;
    lua_getfield(L, -1, "species"); // put mass-fraction table at TOS 
    for ( size_t isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_gettable(L, -2);
	cell->Q->massf[isp] += lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove species table from top-of-stack

    // cout << "Table of sources for individual energies" << endl;
    lua_getfield(L, -1, "energies"); // put energy table at TOS 
    for ( size_t i = 1; i < nmodes; ++i ) {
	lua_pushinteger(L, i);
	lua_gettable(L, -2);
	cell->Q->energies[i] += lua_tonumber(L, -1);
	lua_pop(L, 1); // remove the number to leave the table at TOS
    }
    lua_pop(L, 1); // remove energies table from top-of-stack

    lua_settop(L, 0); // clear the stack
    // cout << "End of add_udf_source_vector_for_cell()" << endl;
    return SUCCESS;
} // end add_udf_source_vector_for_cell()


/// brief set to the components of the vertex velocity, via a Lua udf.
/// This is done (only for moving grid) for inviscid flux calculation.
int add_udf_velocity_for_vtx( Block *bdp, size_t gtl)
{
    // Call the user-defined function which returns a table
    // of vertex velocity.
    // These are added to the inviscid source terms for moving grid
    
    global_data &G = *get_global_data_ptr();

    size_t krangemax = ( G.dimensions == 2 ) ? bdp->kmax : bdp->kmax+1;
    for ( size_t k = bdp->kmin; k <= krangemax; ++k ) {
	for ( size_t j = bdp->jmin; j <= bdp->jmax+1; ++j ) {
	    for ( size_t i = bdp->imin; i <= bdp->imax+1; ++i ) {
	   
	        FV_Vertex *vtx = bdp->get_vtx(i,j,k);
	        
                lua_getglobal(L, "vtx_velocity");  // Lua function to be called
                lua_newtable(L); // creates a table that is now at the TOS
	    
                lua_pushinteger(L, i); lua_setfield(L, -2, "i");
                lua_pushinteger(L, j); lua_setfield(L, -2, "j");
                lua_pushinteger(L, k); lua_setfield(L, -2, "k");
                lua_pushnumber(L, vtx->pos[0].x); lua_setfield(L, -2, "x");
                lua_pushnumber(L, vtx->pos[0].y); lua_setfield(L, -2, "y");
                lua_pushnumber(L, vtx->pos[0].z); lua_setfield(L, -2, "z");
                lua_pushnumber(L, G.sim_time); lua_setfield(L, -2, "t"); 
                lua_pushnumber(L, G.dt_global); lua_setfield(L, -2, "dt");                                
                lua_pushinteger(L, bdp->id); lua_setfield(L, -2, "bdp_id");
                
                int number_args = 1; // table of {i,j,k,bdp_id}
                int number_results = 1; // one table of results returned on the stack.
                if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	            handle_lua_error(L, "error running user vtx_velocity function: %s\n", 
	                             lua_tostring(L, -1));
                }
                
                lua_getfield(L, -1, "vel_x"); vtx->vel[gtl].x = lua_tonumber(L, -1); lua_pop(L, 1);
                lua_getfield(L, -1, "vel_y"); vtx->vel[gtl].y = lua_tonumber(L, -1); lua_pop(L, 1);
                lua_getfield(L, -1, "vel_z"); vtx->vel[gtl].z = lua_tonumber(L, -1); lua_pop(L, 1);
                
                lua_settop(L, 0); // clear the stack
	    }
	}
    }
            
    return SUCCESS;
} // end add_udf_velocity_for_vtx()


//---------------------------------------------------------------------------

int integrate_blocks_in_sequence(void)
// This procedure integrates the blocks two-at-a-time in sequence.
//
// The idea is to approximate the space-marching approach of sm_3d
// and, hopefully, achieve a significant speed-up over integration
// of the full array of blocks.
//
// It is assumed that we have the restricted case 
// of all the blocks making a single line (west to east) with
// supersonic flow at the west face of block 0 and 
// extrapolate_out on the east face of the final block.
//
// The above assumptions still stands but the code has been modified
// such that the calcluation now does two active blocks at a time.
// The calculation initially does blocks 0 and 1, and then moves over
// one block for the next step calculating blocks 1 and 2 next.
// This process is then repeated for the full calculation

{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    Block *bdp;
    double time_slice = G.max_time / (G.nblock - 1);
    BoundaryCondition *bcp_save;
    int status_flag = SUCCESS;
    //---------------------------------------------------------------------------------------
    // Check compatability...
    // This procedure does not work with MPI jobs, nor in the presence of MappedCellBCs.
    if ( G.mpi_parallel ) {
	throw std::runtime_error("Block-sequence integration not implemented with MPI.");
    }
    bool found_mapped_cell_bc = false;
    for ( Block *bdp : G.my_blocks ) {
	int number_faces = (G.dimensions == 3 ? 6: 4);
	for ( int iface = 0; iface < number_faces; ++iface ) {
	    if ( bdp->bcp[iface]->type_code == MAPPED_CELL ) found_mapped_cell_bc = true;
	}
    }
    if ( found_mapped_cell_bc ) {
	throw std::runtime_error("Block-sequence integration not implemented with mapped-cell exchange.");
    }
    //--------------------------------------------------------------------------------------

    // We should be good to run...
    // Initially deactivate all blocks
    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	G.bd[jb].active = false;
    }

    if ( G.verbosity_level >= 1 ) {
	cout << "Integrate Block 0 and Block 1" << endl;
    }

    // Start by setting up block 0.
    bdp = &(G.bd[0]);
    bdp->active = true;
    // Apply the assumed SupINBC to the west face and propogate across the block.
    bdp->bcp[WEST]->apply_convective(0.0);
    bdp->propagate_data_west_to_east( G.dimensions );
    for ( FV_Cell *cp: bdp->active_cells ) {
	cp->encode_conserved(0, 0, bdp->omegaz, with_k_omega);
	// Even though the following call appears redundant at this point,
	// fills in some gas properties such as Prandtl number that is
	// needed for both the cfd_check and the BLomax turbulence model.
	cp->decode_conserved(0, 0, bdp->omegaz, with_k_omega);
    }

    // Now set up block 1
    bdp = &(G.bd[1]);
    bdp->active = true;
    // Save the original east boundary condition and apply the temporary condition
    // ExtrapolateOutBC for the calculation.
    bcp_save = bdp->bcp[EAST];
    bdp->bcp[EAST] = new ExtrapolateOutBC(bdp, EAST, 0);
    // Read in data from block 0 and propagate across the block.
    exchange_shared_boundary_data(1, COPY_FLOW_STATE, 0);
    bdp->propagate_data_west_to_east( G.dimensions );
    for ( FV_Cell *cp: bdp->active_cells ) {
	cp->encode_conserved(0, 0, bdp->omegaz, with_k_omega);
	cp->decode_conserved(0, 0, bdp->omegaz, with_k_omega);
    }

    // Integrate just the first two blocks in time, hopefully to steady state.
    G.bd[0].active = true;
    G.bd[1].active = true;
    integrate_in_time(time_slice);

    // The rest of the blocks.
    for ( size_t jb = 2; jb < (G.nblock); ++jb ) {
	// jb-2, jb-1, jb
	// jb-2 is the block to be deactivated, jb-1 has been iterated and now
	// becomes the left most block and jb is the new block to be iterated
	if ( G.verbosity_level >= 1 ) {
	    cout << "Integrate Block " << jb << endl;
	}
	// Make the block jb-2 inactive.
	bdp = &(G.bd[jb-2]);
	bdp->active = false;

	// block jb-1 - reinstate the previous boundary condition on east face
	// but leave the block active
	bdp = &(G.bd[jb-1]);
	delete bdp->bcp[EAST];
	bdp->bcp[EAST] = bcp_save;

	// Set up new block jb to be integrated
	bdp = &(G.bd[jb]);
	bdp->active = true;
	if ( jb < G.nblock-1 ) {
	    // Cut off the east boundary of the current block 
	    // from the downstream blocks if there are any.
	    bcp_save = bdp->bcp[EAST];
	    bdp->bcp[EAST] = new ExtrapolateOutBC(bdp, EAST, 0);
	}
	// Now copy the starting data into the WEST ghost cells
	// and propagate it across the current block.
	exchange_shared_boundary_data(jb, COPY_FLOW_STATE, 0);
	bdp->propagate_data_west_to_east(G.dimensions);
	for ( FV_Cell *cp: bdp->active_cells ) {
	    cp->encode_conserved(0, 0, bdp->omegaz, with_k_omega);
	    cp->decode_conserved(0, 0, bdp->omegaz, with_k_omega);
	}
	// Integrate just the two currently active blocks in time,
	// hopefully to steady state.
	G.bd[jb-1].active = true;
	G.bd[jb].active = true;
	integrate_in_time(jb*time_slice);
    }
    // Before leaving, we want all blocks active for output.
    for ( size_t jb = 0; jb < G.nblock; ++jb ) {
	G.bd[jb].active = true;
    }
    return status_flag;
} // end integrate_blocks_in_sequence()

//---------------------------------------------------------------------------

int write_solution_data(std::string tindxstring)
// This function only for use below, in the main time-stepping loop.
{
    global_data &G = *get_global_data_ptr();
    int js, final_s;
    char jbcstr[10], jscstr[10];
    std::string foldername = "flow/"+tindxstring;
    std::string jsstring, jbstring, filename;
    ensure_directory_is_present(foldername); // includes Barrier
    for ( Block *bdp : G.my_blocks ) {
	sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr; 
	filename = foldername+"/"+ G.base_file_name+".flow"+jbstring+"."+tindxstring;
	bdp->write_solution(filename, G.sim_time, G.dimensions, zip_files);
    }

    if ( G.moving_grid ) {
	foldername = "grid/"+tindxstring;
	ensure_directory_is_present(foldername); // includes Barrier
	for ( Block *bdp : G.my_blocks ) {
	    sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr; 
	    filename = foldername+"/"+ G.base_file_name+".grid"+jbstring+"."+tindxstring;
	    bdp->write_grid(filename, G.sim_time, G.dimensions, zip_files);
	}
	if ( G.write_vertex_velocities ) {
	    ensure_directory_is_present("vel");
	    foldername = "vel/"+tindxstring;
	    ensure_directory_is_present(foldername); // includes Barrier
	    // Loop over blocks
	    for ( Block *bdp : G.my_blocks ) {
		sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr;
		final_s = ((G.dimensions == 3)? BOTTOM : WEST);
		// Loop over boundaries/surfaces
		for ( js = NORTH; js <= final_s; ++js ) {
		    sprintf( jscstr, ".s%04d", js ); jsstring = jscstr;
		    filename = foldername+"/"+ G.base_file_name+".vel" \
			+jbstring+jsstring+"."+tindxstring;
		    bdp->bcp[js]->write_vertex_velocities(filename, G.sim_time, G.dimensions);
		    if ( zip_files ) do_system_cmd("gzip -f "+filename);
		}
	    } // end for ( Block *bdp
	} // end if ( G.write_vertex_velocities
    } // end if ( G.moving_grid

    // Compute, store and write heat-flux data, if viscous simulation
    if ( with_heat_flux_files && G.viscous ) {
	foldername = "heat/"+tindxstring;
	ensure_directory_is_present(foldername); // includes Barrier
	for ( Block *bdp : G.my_blocks ) {
	    sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr;
	    final_s = ((G.dimensions == 3)? BOTTOM : WEST);
	    // Loop over boundaries/surfaces
	    for ( js = NORTH; js <= final_s; ++js ) {
		sprintf( jscstr, ".s%04d", js ); jsstring = jscstr;
		filename = foldername+"/"+ G.base_file_name+".heat" \
		    +jbstring+jsstring+"."+tindxstring;
		bdp->bcp[js]->compute_surface_heat_flux();
		bdp->bcp[js]->write_surface_heat_flux(filename,G.sim_time);
		if ( zip_files ) do_system_cmd("gzip -f "+filename);
	    }
	} // for *bdp
    }

    // Compute, store and write surface data, if viscous simulation
    if ( with_surface_files ) {
        ensure_directory_is_present("surf"); // includes Barrier
	foldername = "surf/"+tindxstring;
	ensure_directory_is_present(foldername); // includes Barrier
	for ( Block *bdp : G.my_blocks ) {
	    sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr;
	    final_s = ((G.dimensions == 3)? BOTTOM : WEST);
	    // Loop over boundaries/surfaces
	    for ( js = NORTH; js <= final_s; ++js ) {
		sprintf( jscstr, ".s%04d", js ); jsstring = jscstr;
		filename = foldername+"/"+ G.base_file_name+".surf" \
		    +jbstring+jsstring+"."+tindxstring;
		bdp->bcp[js]->write_surface_data(filename,G.sim_time);
		if ( zip_files ) do_system_cmd("gzip -f "+filename);
	    }
	} // for *bdp
    }
    return SUCCESS;
} // end write_solution_data()

int write_temp_solution_data()
// This function only for flow induced moving grid.
{
    global_data &G = *get_global_data_ptr();
    char jbcstr[10];
    ensure_directory_is_present("temp");
    std::string foldername = "temp/flow";
    std::string jbstring, filename;
    ensure_directory_is_present(foldername); // includes Barrier
    for ( Block *bdp : G.my_blocks ) {
	sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr; 
	filename = foldername+"/"+ G.base_file_name+".flow"+jbstring;
	bdp->write_solution(filename, G.sim_time, G.dimensions, zip_files);
    }

    foldername = "temp/grid";
    ensure_directory_is_present(foldername); // includes Barrier
    for ( Block *bdp : G.my_blocks ) {
	sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr; 
	filename = foldername+"/"+ G.base_file_name+".grid"+jbstring;
	bdp->write_grid(filename, G.sim_time, G.dimensions, zip_files);
    }
    
    if ( master ) {
        FILE *tempfile;
        std::string tempname = "temp/control.dat";    
        tempfile = fopen(tempname.c_str(), "w");
        fprintf(tempfile, "# Control parameters from the simulation\n");
        fprintf(tempfile, "# simulation time, time step\n");    
        fprintf(tempfile, "%20.12e %20.12e\n", G.sim_time, G.dt_global);
        fclose(tempfile);
    }

    return SUCCESS;
} // end write_temp_solution_data()

int integrate_in_time(double target_time)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t n_active_blocks;
    string jbstring, jsstring, tindxstring;
    char jbcstr[10], tindxcstr[10];
    string filename, commandstring, foldername;
    std::vector<double> dt_record;
    double stopping_time;
    bool finished_time_stepping;
    bool viscous_terms_are_on;
    bool diffusion_terms_are_on;
    bool do_cfl_check_now;
    int cfl_result;
#   ifdef _MPI
    int cfl_result_2;
#   endif
    int status_flag = SUCCESS;
    dt_record.resize(G.my_blocks.size()); // Just the blocks local to this process.

    if ( G.verbosity_level >= 1 && master ) {
	printf( "Integrate in time\n" );
	fflush(stdout);
    }
    if ( target_time <= 0.0 ) {
	stopping_time = G.max_time;
    } else {
	stopping_time = target_time;
    }
    // The next time for output...
    G.t_plot = G.sim_time + G.dt_plot;
    G.t_his = G.sim_time + G.dt_his;
    G.t_moving = 0.0;
    // if this is fluid structure interaction simulation
    if ( G.flow_induced_moving ) {
        G.t_moving = G.sim_time + G.dt_moving;
    }
    
    // Flags to indicate that the saved output is fresh.
    // On startup or restart, it is assumed to be so.
    output_just_written = true;
    history_just_written = true;
    av_output_just_written = true;

    G.step = 0; // Global Iteration Count
    G.dt_acc = 0.0; // Initialise to 0.0 to begin acculumating dt increments for conjugate ht problem
    do_cfl_check_now = false;
    if ( G.heat_time_stop == 0.0 ) {
	// We don't want heating at all.
	G.heat_factor = 0.0;
    } else if ( G.sim_time >= G.heat_time_start && G.sim_time < G.heat_time_stop ) {
	// Apply full heat-source effects because both time limits are set
	// and we are within them.
	G.heat_factor = 1.0;
    } else {
	// Looks like we want heating at some period in time but it is not now.
	G.heat_factor = 0.0;
    }
    if ( G.ignition_time_stop == 0.0 ) {
	// We don't want ignition zone applied.
	G.ignition_zone_active = false;
    } else if ( G.sim_time >= G.ignition_time_start && G.sim_time < G.ignition_time_stop ) {
	// Apply ignition zone because both time limits are set
	// and we are within them.
	G.ignition_zone_active = true;
    } else {
	// Looks like we want ignition zone at some period in time but it is not now.
	G.ignition_zone_active = false;
    }
    if ( G.viscous ) {
	// We have requested viscous effects but they may be delayed.
	if ( G.viscous_time_delay > 0.0 && G.sim_time < G.viscous_time_delay ) {
	    // We will initially turn-down viscous effects and
	    // only turn them up when the delay time is exceeded.
	    G.viscous_factor = 0.0;
	    viscous_terms_are_on = false;
	} else {
	    // No delay in applying full viscous effects.
	    G.viscous_factor = 1.0;
	    viscous_terms_are_on = true;
	}
    } else {
	// We haven't requested viscous effects at all.
	G.viscous_factor = 0.0;
	viscous_terms_are_on = false;
    }

    if ( G.diffusion ) {
	// We have requested diffusion effects but they may be delayed.
	if ( G.diffusion_time_delay > 0.0 && G.sim_time < G.diffusion_time_delay ) {
	    // We will initially turn-down diffusion effects and
	    // only turn them up when the delay time is exceeded.
	    G.diffusion_factor = 0.0;
	    diffusion_terms_are_on = false;
	} else {
	    // No delay in applying full diffusion effects.
	    G.diffusion_factor = 1.0;
	    diffusion_terms_are_on = true;
	}
    } else {
	// We haven't requested diffusion effects at all.
	G.diffusion_factor = 0.0;
	diffusion_terms_are_on = false;
    }

    // Spatial filter may be applied occasionally.
    if ( G.filter_flag ) {
	if ( G.sim_time > G.filter_tstart ) {
	    G.filter_next_time = G.sim_time + G.filter_dt;
	    if ( G.filter_next_time > G.filter_tend ) G.filter_flag = false;
	} else {
	    G.filter_next_time = G.filter_tstart;
	}
    }

    // Store the initial temperature at the monitor points.
    for ( Block *bdp : G.my_blocks ) {
	for ( size_t im = 0; im < bdp->mncell; ++ im ) {
	    FV_Cell *cp = bdp->get_cell(bdp->micell[im]+bdp->imin, 
					bdp->mjcell[im]+bdp->jmin, 
					bdp->mkcell[im]+bdp->kmin);
	    bdp->initial_T_value[im] = cp->fs->gas->T[0];
#           if 0
	    // Some Debug...
	    cout << "block id=" << bdp->id << " cell at x=" << cp->pos[0].x << " y=" << cp->pos[0].y
		 << " z=" << cp->pos[0].z << " initial T=" << bdp->initial_T_value[im] << endl;
#           endif
	}
    }

    // Normally, we can terminate upon either reaching 
    // a maximum time or upon reaching a maximum iteration count.
    finished_time_stepping = (G.sim_time >= stopping_time || G.step >= G.max_step);

    //----------------------------------------------------------------
    //                 Top of main time-stepping loop
    //----------------------------------------------------------------
    while ( !finished_time_stepping ) {
#       ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#       endif    
	if ( (G.step/G.control_count)*G.control_count == G.step ) {
	    // Reparse the time-step control parameters as frequently as specified.
	    read_control_parameters(G.base_file_name+".control", master, false);
	}
	// One of the control parameters is G.max_time, so we need to
	// make sure that the stopping_time is updated accordingly.
	if ( target_time <= 0.0 ) {
	    stopping_time = G.max_time;
	} else {
	    stopping_time = target_time;
	}
        // Continue taking steps so long as there is at least one active block.
#       ifdef _MPI
        // FIX-ME -- at the moment we assume all blocks are active
	// in the distributed-memory version.
	n_active_blocks = G.nblock;
#       else
        n_active_blocks = 0;
	for ( Block *bdp : G.my_blocks ) {
            if ( bdp->active ) ++n_active_blocks;
        }
#       endif
        if ( n_active_blocks == 0 ) {
	    printf( "There are no active blocks at step %d\n", static_cast<int>(G.step) );
	    status_flag = FAILURE;
            break;
        }

        // 0. Alter configuration setting if necessary.
	if ( G.viscous && !viscous_terms_are_on && G.sim_time >= G.viscous_time_delay ) {
	    // We want to turn on the viscous effects only once (if requested)
	    // and, when doing so in the middle of a simulation, 
	    // reduce the time step to ensure that these new terms
	    // do not upset the stability of the calculation.
	    if ( G.verbosity_level >= 1 ) printf( "Turning on viscous effects.\n" );
	    viscous_terms_are_on = true;
	    if ( G.step > 0 ) {
		if ( G.verbosity_level >= 1 ) printf( "Start softly with viscous effects.\n" );
		G.viscous_factor = 0.0;  
		G.dt_global *= 0.2;
	    } else {
		if ( G.verbosity_level >= 1 ) printf( "Start with full viscous effects.\n" );
		G.viscous_factor = 1.0;
	    }
	}
	if ( viscous_terms_are_on && G.viscous_factor < 1.0 ) {
	    incr_viscous_factor(G.viscous_factor_increment);
	    if ( G.verbosity_level >= 1 ) printf( "Increment viscous_factor to %f\n", G.viscous_factor );
	    do_cfl_check_now = true;
	}
	if ( G.diffusion && !diffusion_terms_are_on && G.sim_time >= G.diffusion_time_delay ) {
	    // We want to turn on the diffusion effects only once (if requested)
	    // and, when doing so in the middle of a simulation,
	    // reduce the time step to ensure that these new terms
	    // do not upset the stability of the calculation.
	    if ( G.verbosity_level >= 1 ) printf( "Turning on diffusion effects.\n" );
	    diffusion_terms_are_on = true;
	    if ( G.step > 0 ) {
		if ( G.verbosity_level >= 1 ) printf( "Start softly with diffusion effects.\n" );
		G.diffusion_factor = 0.0;
		G.dt_global *= 0.2;
	    } else {
		if ( G.verbosity_level >= 1 ) printf( "Start with full diffusion effects.\n" );
		G.diffusion_factor = 1.0;
	    }
	}
	if ( diffusion_terms_are_on && G.diffusion_factor < 1.0 ) {
	    incr_diffusion_factor( G.diffusion_factor_increment );
	    if ( G.verbosity_level >= 1 ) printf( "Increment diffusion_factor to %f\n", G.diffusion_factor );
	    do_cfl_check_now = true;
	}

	if ( G.heat_time_stop > 0.0 ) {
	    // We want heating at some time.
	    if ( G.sim_time >= G.heat_time_start && G.sim_time < G.heat_time_stop ) {
		// We are within the period of heating but
		// we may be part way through turning up the heat gradually.
		if ( G.heat_factor < 1.0 ) {
		    incr_heat_factor(G.heat_factor_increment);
		    if ( G.verbosity_level >= 1 ) printf("Increment heat_factor to %f\n", G.heat_factor);
		    do_cfl_check_now = true;
		} 
	    } else {
		// We are outside the period of heating.
		G.heat_factor = 0.0;
	    }
	}

	if ( G.ignition_time_stop > 0.0 ) {
            if ( G.sim_time >= G.ignition_time_start && G.sim_time < G.ignition_time_stop ) {
	        G.ignition_zone_active = true;
            } else {
	        G.ignition_zone_active = false;
            }
        }

	// 0.a. call out to user-defined function
	if ( G.udf_file.length() > 0 ) {
	    call_udf( G.sim_time, G.step, "at_timestep_start" );
	}
	// 0.b. We need to remember some flow information for later computing residuals.
	for ( Block *bdp : G.my_blocks ) {
	    bdp->init_residuals( G.dimensions );
	}

	// 1. Set the size of the time step.
	if ( G.step == 0 ) {
	    // When starting a new calculation,
	    // set the global time step to the initial value.
	    do_cfl_check_now = false;
		// if we are using sequence_blocks we don't want to reset the dt_global for each block
		if ( G.sequence_blocks && G.dt_global != 0 ) {
		    /* do nothing i.e. keep dt_global from previous block */ ;
		} else { 
		    G.dt_global = G.dt_init;
		}
	} else if ( !G.fixed_time_step && (G.step/G.cfl_count)*G.cfl_count == G.step ) {
	    // Check occasionally 
	    do_cfl_check_now = true;
	} // end if (G.step == 0 ...

	if ( do_cfl_check_now ) {
	    // Adjust the time step to be the minimum allowed
	    // for any active block. 
	    G.dt_allow = 1.0e6; /* outrageously large so that it gets replace immediately */
	    G.cfl_max = 0.0;
	    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
		Block *bdp = G.my_blocks[jb];
		if ( bdp->active ) {
		    cfl_result = bdp->determine_time_step_size();
		    if ( cfl_result != 0 ) {
			program_return_flag = DT_SEARCH_FAILED;
			status_flag = FAILURE;
			goto conclusion;
		    }
		}
	    } // end for jb loop
	    // If we arrive here, cfl_result will be zero, indicating that all local blocks 
	    // have successfully determined a suitable time step.
#           ifdef _MPI
	    MPI_Allreduce( &cfl_result, &cfl_result_2, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
	    if ( cfl_result_2 != 0 ) {
		program_return_flag = DT_SEARCH_FAILED;
		status_flag = FAILURE;
		goto conclusion;
	    }
#           endif
	    // Get an overview of the allowable timestep.
	    // Note that, for an MPI job, jb may not be the same as bdp->id.
	    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
		dt_record[jb] = 0.0;
	    }
	    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
		Block *bdp = G.my_blocks[jb];
		if ( !bdp->active ) continue;
		if ( bdp->dt_allow < G.dt_allow ) G.dt_allow = bdp->dt_allow;
		dt_record[jb] = bdp->dt_allow;
		if ( bdp->cfl_max > G.cfl_max ) G.cfl_max = bdp->cfl_max;
		if ( bdp->cfl_min < G.cfl_min ) G.cfl_min = bdp->cfl_min;
	    }
#           ifdef _MPI
	    // Finding the minimum allowable time step is very important for stability.
	    MPI_Allreduce(MPI_IN_PLACE, &(G.dt_allow), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, &(G.cfl_max), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, &(G.cfl_min), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#           endif
	    // Change the actual time step, as needed.
	    if ( G.dt_allow <= G.dt_global ) {
		// If we need to reduce the time step, do it immediately.
		G.dt_global = G.dt_allow;
	    } else {
		// Make the transitions to larger time steps gentle.
		G.dt_global = min(G.dt_global*1.5, G.dt_allow);
		// The user may supply, explicitly, a maximum time-step size.
		G.dt_global = min(G.dt_global, G.dt_max);
	    }
	    do_cfl_check_now = 0;  // we have done our check for now
	} // end if do_cfl_check_now 

        if ( (G.step > 5) && (G.cfl_max < G.cfl_tiny) ) {
            // Record the worst CFL
            G.cfl_tiny = G.cfl_max;
            G.time_tiny = G.sim_time;
        }
        if ( G.cfl_max > 100.0 ) {
            // If the CFL has gone crazy, bail out.
            printf( "\nCFL = %e: breaking main loop\n", G.cfl_max );
	    status_flag = FAILURE;
            break;
        }

	if (G.MHD) {
	    // update the divergence cleaning speed
	    update_MHD_c_h();
	}

        // 2. Attempt a time step.
	// 2aa. Compute wall conduction if conjugate heat transfer available
	if ( G.wall_update_count == 0 ) {
	    cout << "ERROR: The 'wall_update_count' should be an integer greater than or equal to 1.\n";
	    cout << "ERROR: Check value set in .control file.\n";
	    cout << "BAILING OUT!\n";
	    exit(BAD_INPUT_ERROR);
	}
#       ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#       endif
	if ( G.conjugate_ht_active && ((G.step/G.wall_update_count)*G.wall_update_count == G.step) ) {
	    int flag = gather_near_wall_gas_data(G);
	    if ( flag != SUCCESS ) {
		cout << "Error gathering temperatures and thermal conductivities near the wall from all ranks.\n";
		cout << "Bailing out!\n";
		exit(FAILURE);
	    }
	    // Now only master has properly updated T and q vector of near wall temperatures in gas
	    // so we'll only do our work on the master.
	    if ( master ) {
		// Compute fluxes at interface between gas and solid.
		// 1. Gather temperatures and thermal conductivities from solid
		flag = gather_near_wall_solid_data(*(G.wm), G.T_solid_near_wall, G.k_solid_near_wall);
		if ( flag != SUCCESS ) {
		    cout << "Error retrieving near wall temperatures from wall conduction solver.\n";
		    cout << "Bailing out!\n";
		    exit(FAILURE);
		}
		// 2. Given the temperatures, find the interface flux
		flag = G.conj_ht_iface->compute_flux(G.T_gas_near_wall, G.k_gas_near_wall,
						     G.T_solid_near_wall, G.k_solid_near_wall, 
						     G.q_wall, G.T_wall);	//Changed by JB 02-08-14
		if ( flag != SUCCESS ) {
		    cout << "Error computing interface fluxes at gas/solid interface.\n";
		    cout << "Bailing out!\n";
		    exit(FAILURE);
		}
		// 3. Now do the wall update
		if ( G.dt_acc == 0.0 ) {
		    // Just started but need to do update.
		    G.dt_acc = G.dt_global;
		}
		flag = update_state(*(G.wm), G.dt_acc, G.q_wall, G.T_wall); //Changed by JB 02-08-14
		// Reset dt_acc to 0.0
		G.dt_acc = 0.0;
	    }
#           ifdef _MPI
	    // For MPI version:
	    // Now broadcast the fluxes so that every rank
	    // has an update copy of the vector.
	    MPI_Barrier(MPI_COMM_WORLD);
	    flag = broadcast_wall_values(G);
	    if ( flag != SUCCESS ) {
		cout << "Error broadcasting the updated wall values to all ranks.\n";
		cout << "Bailing out!\n";
		exit(FAILURE);
	    }
#           endif
	    // For shared memory version:
	    // In shared memory version, G.q_wall is up-to-date globally
	    // and available for all blocks.
	}
	
	// 2a.
	// Assgin moving vertex velocity
	if ( G.moving_grid ) {
	    // The user might want to do some customization before the grid update step.
	    // If so, we provide a user-defined function as the customization point.
	    
//	    if ( G.udf_file.length() > 0 ) {
//		for ( Block *bdp : G.my_blocks ) {
//		    call_udf2( G.sim_time, G.step, bdp->id, "before_grid_update" );
//		}
//	    }
//      code has been moved to after G.sim_time >= G.t_moving.

	    // Place a barrier so that each rank has completed the call before moving on
#           ifdef _MPI	        
	    MPI_Barrier( MPI_COMM_WORLD );
#           endif	      
	    // Grid-movement is done after a specified point in time.
	    // As Paul Petrie-repar suggested
            // The edge length (2D) or interface area (3D), interface velocity, normal vector
            // should be remained as the same as level 0 for both stages of the pc
            // The default gird movement strategy is set by equations, the more complex way
            // is flow induced moving grid, which vertex moving velocity will be defined
            // the external code 
            if ( G.flow_induced_moving ) {
                for ( Block *bdp : G.my_blocks ) {
                    bdp->clear_vertex_velocities(G.dimensions);
                }
            } // end if ( G.flow_induced_moving )   
	    if ( G.sim_time >= G.t_moving ) { 
	        // The user might want to do some customization before the grid update step.
	        // If so, we provide a user-defined function as the customization point.
 
	        //if ( G.udf_vtx_velocity_flag == 1 ) {  
	        if ( G.udf_file.length() > 0 ) {
		        for ( Block *bdp : G.my_blocks ) {
		            call_udf2( G.sim_time, G.step, bdp->id, "before_grid_update" );
		        }
	        }
	        // Place a barrier so that each rank has completed the call before moving on
#               ifdef _MPI	        
	        MPI_Barrier( MPI_COMM_WORLD );
#               endif	  

	        if ( G.flow_induced_moving ) { // flow induced grid movement
	            write_temp_solution_data();
#                   ifdef _MPI	        
	            MPI_Barrier( MPI_COMM_WORLD );
#                   endif	        
	            if ( master ) {
	                if ( system("./flow_induced_moving.sh") != 0 ) {
	                    printf( "Error: check the external code for flow induced moving grid.\n");
	                }
	            }
#                   ifdef _MPI	        
	            MPI_Barrier( MPI_COMM_WORLD );
#                   endif
	            for ( Block *bdp : G.my_blocks ) {
	                sprintf( jbcstr, "b%04d", static_cast<int>(bdp->id) );
	                string jbstring = jbcstr;
	                string filename = "temp/"+jbstring;
                        if ( bdp->read_vertex_velocities(filename, G.dimensions) != SUCCESS ) {
                            printf( "Error: check the output vertex velocities.\n");
	                }
                    } // end for *bdp
#                   ifdef _MPI	        
	            MPI_Barrier( MPI_COMM_WORLD );
#                   endif
                    // check CFL condition if flow induced moving and
                    // grid has been adjusted at this time step                    
                    do_cfl_check_now = true;
                    // if we are using fixed time step, no cfl check
                    if ( G.fixed_time_step ) {
                        do_cfl_check_now = false;
                    }
	        } // end if ( G.flow_induced_moving )
	        else { // user-defined grid movement
	            for ( Block *bdp : G.my_blocks ) {
	                if ( G.udf_vtx_velocity_flag == 1 ) { // typical way for moving grid 
	                    add_udf_velocity_for_vtx(bdp, 0);
	                } else { // mainly for shock fitting
		            bdp->set_geometry_velocities(G.dimensions, 0);
		        }
	            }  // end for
                } // end else
                G.t_moving = G.sim_time + G.dt_moving;
	    } // end if ( G.sim_time >= G.t_moving )
	} // end if ( G.moving_grid )
		
	// explicit or implicit update of the inviscid terms.
	int break_loop2 = 0;
	switch ( G.implicit_mode ) {
	case 0: // explicit update of convective terms and, maybe, the viscous terms
	    if ( G.moving_grid )     
		break_loop2 = gasdynamic_increment_with_moving_grid(G.dt_global);
	    else
		break_loop2 = gasdynamic_explicit_increment_with_fixed_grid(G.dt_global);
	    break;
	case 1:
	    break_loop2 = gasdynamic_point_implicit_inviscid_increment(G.dt_global);
	    break;
	case 2:
	    break_loop2 = gasdynamic_fully_implicit_inviscid_increment(G.dt_global);
	} // end switch ( G.implicit_mode )
	if ( break_loop2 ) {
	    printf("Breaking main loop:\n");
	    printf("    time step failed at inviscid increment.\n");
	    status_flag = FAILURE;
	    break;
	}
        
	// 2b. Piston step.
	// Code removed 24-Mar-2013.

	// 2b. Recalculate all geometry if moving grid.
	if ( G.moving_grid ) {
	    for ( Block *bdp : G.my_blocks ) {
		bdp->compute_primary_cell_geometric_data(G.dimensions, 0);
		bdp->compute_distance_to_nearest_wall_for_all_cells(G.dimensions, 0);
		bdp->compute_secondary_cell_geometric_data(G.dimensions, 0);
	    }
	 }
	
	// 2c. Increment because of viscous effects may be done
	//     separately to the convective terms.
	if ( G.viscous && G.separate_update_for_viscous_terms ) {
	    // We now have the option of explicit or point implicit update
	    // of the viscous terms, thanks to Ojas Joshi, EPFL.
	    int break_loop = 0;
	    switch ( G.implicit_mode ) {
	    case 0:
		break_loop = gasdynamic_separate_explicit_viscous_increment();
		break;
	    case 1:
		break_loop = gasdynamic_point_implicit_viscous_increment();
		break;
	    case 2:
	    	break_loop = gasdynamic_fully_implicit_viscous_increment();
	    } // end switch ( G.implicit_mode )
	    if ( break_loop ) {
		printf("Breaking main loop:\n");
		printf("    time step failed at viscous increment.\n");
		status_flag = FAILURE;
		break;
	    }
	    if ( G.turbulence_model == TM_K_OMEGA && G.separate_update_for_k_omega_source ) {
		for ( Block *bdp : G.my_blocks ) {
		    if ( !bdp->active ) continue;
		    for ( FV_Cell *cp: bdp->active_cells )
			cp->update_k_omega_properties(G.dt_global);
		}
	    }
	} // end if ( G.viscous )

        // 2d. Chemistry step. 
	//     Allow finite-rate evolution of species due
        //     to chemical reactions
#ifdef GPU_CHEM
        if ( G.reacting && G.sim_time >= G.reaction_time_start ) {
	    // Gather all cells across all active blocks
	    vector<FV_Cell*> cells;
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		for ( FV_Cell *cp: bdp->active_cells ) {
		    cells.push_back(cp);
		}
	    }
	    update_chemistry(*gchem, G.dt_global, cells);
	}
#else
        if ( G.reacting && G.sim_time >= G.reaction_time_start ) {
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		for ( FV_Cell *cp: bdp->active_cells ) {
#ifdef GPU_CHEM_ALGO
		    if ( cp->chemical_increment(G.dt_global) != SUCCESS ) {
			cout << "Chemistry problem using simplified stepping algorithm.\n";
			cout << "Bailing out!\n";
			exit(NUMERICAL_ERROR);
		    }
#else    
		    if ( cp->chemical_increment(G.dt_global, G.T_frozen) != SUCCESS ) {
			cout << "In block: " << bdp->id << " the chemical increment failed on cell:\n";
			vector<size_t> ijk(bdp->to_ijk_indices(cp->id));
			cout << "[i,j,k]= [" << ijk[0] << "," << ijk[1] << "," << ijk[2] << "]\n";
			cout << "The global timestep was: " << G.dt_global << endl;
			cout << "The chemistry timestep was: " << cp->dt_chem << endl;
			cout << "Bailing out at this point!\n";
			exit(NUMERICAL_ERROR);
		    }
#endif
		}
	    }
	}
#endif
	// 2e. Thermal step.
	//     Allow finite-rate evolution of thermal energy
	//     due to transfer between thermal energy modes.
	if ( G.thermal_energy_exchange && G.sim_time >= G.reaction_time_start  ) {
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		for ( FV_Cell *cp: bdp->active_cells ) {
		    if ( cp->thermal_increment(G.dt_global, G.T_frozen_energy) != SUCCESS ) {
			cout << "In block: " << bdp->id << " the thermal increment failed on cell:\n";
			vector<size_t> ijk(bdp->to_ijk_indices(cp->id));
			cout << "[i,j,k]= [" << ijk[0] << "," << ijk[1] << "," << ijk[2] << "]\n";
			cout << "The global timestep was: " << G.dt_global << endl;
			cout << "The thermal timestep was: " << cp->dt_therm << endl;
			cout << "Bailing out at this point!\n";
			exit(NUMERICAL_ERROR);
		    }
		}
	    }
	}

        // 3. Update the time record and (occasionally) print status.
        ++G.step;
        output_just_written = false;
        history_just_written = false;
	av_output_just_written = false;
	G.dt_acc += G.dt_global;
	// G.sim_time += G.dt_global; 2013-04-07 have moved increment of sim_time
	// into the inviscid gasdynamic update

        if ( ((G.step / G.print_count) * G.print_count == G.step) && master ) {
            // Print the current time-stepping status.
            now = time(NULL);
            if ( G.verbosity_level >= 1 ) {
		printf("Step=%7d t=%10.3e dt=%10.3e %s\n",
		       static_cast<int>(G.step), G.sim_time, G.dt_global,
		       time_to_go(start, now, G.step, G.max_step,
				  G.dt_global, G.sim_time, G.max_time) );
	    }
            fprintf(G.logfile, "Step=%7d t=%10.3e dt=%10.3e %s\n",
		    static_cast<int>(G.step), G.sim_time, G.dt_global,
		    time_to_go(start, now, G.step, G.max_step, G.dt_global, G.sim_time, G.max_time) );
            fprintf(G.logfile, "CFL_min = %e, CFL_max = %e, dt_allow = %e\n",
		    G.cfl_min, G.cfl_max, G.dt_allow );
            fprintf(G.logfile, "Smallest CFL_max so far = %e at t = %e\n",
		    G.cfl_tiny, G.time_tiny );
	    for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
		Block *bdp = G.my_blocks[jb];
                if ( !bdp->active ) continue;
                fprintf(G.logfile, " dt[%d]=%e", static_cast<int>(bdp->id), dt_record[jb] );
            }
	    if ( n_active_blocks == 1 ) {
		fprintf(G.logfile, "\nThere is %d active block.\n",
			static_cast<int>(n_active_blocks) );
	    } else {
		fprintf(G.logfile, "\nThere are %d active blocks.\n",
			static_cast<int>(n_active_blocks) );
	    }
	    fflush( stdout );
        } // end if


        // 4. (Occasionally) Write out an intermediate solution
        if ( (G.sim_time >= G.t_plot) && !output_just_written ) {
	    ++output_counter;
	    if ( master ) {
	        fprintf( G.timestampfile, "%04d %e %e\n", static_cast<int>(output_counter),
			 G.sim_time, G.dt_global );
		fflush( G.timestampfile );
	    }
	    sprintf( tindxcstr, "t%04d", static_cast<int>(output_counter) ); // C string
	    write_solution_data(tindxcstr);
	    if ( G.conjugate_ht_active && master ) {
		write_solution(*(G.wm), G.sim_time, static_cast<int>(output_counter));
	    }
	    output_just_written = true;
            G.t_plot += G.dt_plot;
        }
	if ( (G.write_at_step > 0) && (G.step == G.write_at_step) && !write_at_step_has_been_done ) {
	    // Write the solution once-off, most likely for debug.
	    write_solution_data("txxxx");
	    write_at_step_has_been_done = true;
	}
        if ( (G.sim_time >= G.t_his) && !history_just_written ) {
	    for ( Block *bdp : G.my_blocks ) {
		sprintf(jbcstr, ".b%04d", static_cast<int>(bdp->id)); jbstring = jbcstr;
		filename = "hist/"+G.base_file_name+".hist"+jbstring;
                bdp->write_history(filename, G.sim_time);
		bdp->print_forces(G.logfile, G.sim_time, G.dimensions);
		for ( int iface: bdp->transient_profile_faces ) {
		    filename = "./" + G.base_file_name + ".blk"+jbstring+"."+get_face_name(iface)+".profile";
		    bdp->write_profile(filename, iface, G.sim_time, false);
		}
	    }
            history_just_written = true;
            G.t_his += G.dt_his;
        }

	// Velocity profile recording (Andrew Denman  19-June-2003)
	// Code removed 24-Mar-2013
	    
        // 5. For steady-state approach, check the residuals for mass and energy.
        if ( (G.step / G.print_count) * G.print_count == G.step ) {
            G.mass_residual = 0.0;
            G.energy_residual = 0.0;
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		bdp->compute_residuals(G.dimensions, 0);
		fprintf( G.logfile, "RESIDUAL mass block %d max: %e at (%g,%g,%g)\n",
			 static_cast<int>(bdp->id), bdp->mass_residual, bdp->mass_residual_loc.x,
			 bdp->mass_residual_loc.y, bdp->mass_residual_loc.z );
		fprintf( G.logfile, "RESIDUAL energy block %d max: %e at (%g,%g,%g)\n",
			 static_cast<int>(bdp->id), bdp->energy_residual, bdp->energy_residual_loc.x,
			 bdp->energy_residual_loc.y, bdp->energy_residual_loc.z );
            } // end for *bdp
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		if ( bdp->mass_residual > G.mass_residual ) G.mass_residual = bdp->mass_residual; 
		if ( bdp->energy_residual > G.energy_residual ) G.energy_residual = bdp->energy_residual; 
            }
#           ifdef _MPI
	    MPI_Allreduce(MPI_IN_PLACE, &(G.mass_residual), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	    MPI_Allreduce(MPI_IN_PLACE, &(G.energy_residual), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#           endif
            fprintf( G.logfile, "RESIDUAL mass global max: %e step %d time %g\n",
		     G.mass_residual, static_cast<int>(G.step), G.sim_time );
            fprintf( G.logfile, "RESIDUAL energy global max: %e step %d time %g\n",
		     G.energy_residual, static_cast<int>(G.step), G.sim_time );
	    fflush( G.logfile );
        }

	// 6. Spatial filter may be applied occasionally.
	if ( G.filter_flag && G.sim_time > G.filter_next_time ) {
	    // [TODO] Note that this ins't implemented for MappedCellBCs...
	    for ( size_t ipass = 0; ipass < G.filter_npass; ++ipass ) {
#               ifdef _MPI
		MPI_Barrier(MPI_COMM_WORLD);
		mpi_exchange_boundary_data(COPY_FLOW_STATE, 0);
#               else
		for ( Block *bdp : G.my_blocks ) {
		    if ( !bdp->active ) continue;
		    exchange_shared_boundary_data(bdp->id, COPY_FLOW_STATE, 0);
		}
#               endif
		for ( Block *bdp : G.my_blocks ) {
		    if ( !bdp->active ) continue;
		    bdp->apply_spatial_filter_diffusion( G.filter_mu, G.filter_npass, G.dimensions );
		    for ( FV_Cell *cp: bdp->active_cells )
			cp->encode_conserved(0, 0, bdp->omegaz, with_k_omega);
		}
		for ( Block *bdp : G.my_blocks ) {
		    if ( !bdp->active ) continue;
		    apply_convective_bc(*bdp, G.sim_time, G.dimensions);
		    if ( G.viscous ) apply_viscous_bc(*bdp, G.sim_time, G.dimensions); 
		}
#               ifdef _MPI
		MPI_Barrier(MPI_COMM_WORLD);
		mpi_exchange_boundary_data(COPY_FLOW_STATE, 0);
#               else
		for ( Block *bdp : G.my_blocks ) {
		    if ( !bdp->active ) continue;
		    exchange_shared_boundary_data(bdp->id, COPY_FLOW_STATE, 0);
		}
#               endif
		for ( Block *bdp : G.my_blocks ) {
		    if ( !bdp->active ) continue;
		    bdp->apply_spatial_filter_anti_diffusion(G.filter_mu, G.filter_npass, G.dimensions);
		    for ( FV_Cell *cp: bdp->active_cells )
			cp->encode_conserved(0, 0, bdp->omegaz, with_k_omega);
		}
		for ( Block *bdp : G.my_blocks ) {
		    if ( !bdp->active ) continue;
		    apply_convective_bc( *bdp, G.sim_time, G.dimensions );
		    if ( G.viscous ) apply_viscous_bc(*bdp, G.sim_time, G.dimensions); 
		}
	    } // end for ipass
	    G.filter_next_time = G.sim_time + G.filter_dt;
	    if ( G.filter_next_time > G.filter_tend ) G.filter_flag = false;
	} // end if ( G.filter_flag )

	// 7. Call out to user-defined function.
	if ( G.udf_file.length() > 0 ) {
	    call_udf( G.sim_time, G.step, "at_timestep_end" );
	}


        // 8. Loop termination criteria:
        //    (1) reaching a maximum simulation time or target time
        //    (2) reaching a maximum number of steps
        //    (3) finding that the "halt_now" parameter has been set 
	//        in the control-parameter file.
        //        This provides a semi-interactive way to terminate the 
        //        simulation and save the data.
	//    (4) Exceeding a maximum number of wall-clock seconds.
	//    (5) Having the temperature at one of the control points exceed 
	//        the preset tolerance.  
	//        This is mainly for the radiation-coupled simulations.
	//    (-) Exceeding an allowable delta(f_rad) / f_rad_org factor
	//
	//    Note that the max_time and max_step control parameters can also
	//    be found in the control-parameter file (which may be edited
	//    while the code is running).
        if ( G.sim_time >= stopping_time ) {
            finished_time_stepping = true;
            if ( G.verbosity_level >= 1 && master )
		printf( "Integration stopped: reached maximum simulation time.\n" );
        }
        if ( G.step >= G.max_step ) {
            finished_time_stepping = true;
            if ( G.verbosity_level >= 1 && master )
		printf( "Integration stopped: reached maximum number of steps.\n" );
        }
        if ( G.halt_now == 1 ) {
            finished_time_stepping = true;
            if ( G.verbosity_level >= 1 && master )
		printf( "Integration stopped: Halt set in control file.\n" );
        }
	now = time(NULL);
	if ( max_wall_clock > 0 && ( static_cast<int>(now - start) > max_wall_clock ) ) {
            finished_time_stepping = true;
            if ( G.verbosity_level >= 1 && master )
		printf( "Integration stopped: reached maximum wall-clock time.\n" );
	}
	if ( G.halt_on_large_flow_change ) {
	    // Test the monitor points and see if any has experienced a large change in temperature.
	    int large_T_change = 0;
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		for ( size_t im = 0; im < bdp->mncell; ++ im ) {
		    FV_Cell *cp = bdp->get_cell(bdp->micell[im]+bdp->imin,
						bdp->mjcell[im]+bdp->jmin,
						bdp->mkcell[im]+bdp->kmin);
		    if ( fabs(bdp->initial_T_value[im] - cp->fs->gas->T[0]) > G.tolerance_in_T ) {
			large_T_change = 1;
			if ( G.verbosity_level >= 1 ) {
			    cout << "block=" << bdp->id 
				 << " cell at x=" << cp->pos[0].x << " y=" << cp->pos[0].y 
				 << " z=" << cp->pos[0].z << ","
				 << " initial T=" << bdp->initial_T_value[im] 
				 << " current T=" << cp->fs->gas->T[0] << endl;
			    cout << "Temperature change exceeded tolerance of " 
				 << G.tolerance_in_T << endl;
			}
		    }
		} // for im
	    } // for Block
#           ifdef _MPI
	    MPI_Allreduce(MPI_IN_PLACE, &(large_T_change), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#           endif
	    if ( large_T_change == 1 ) {
		finished_time_stepping = true;
		if ( G.verbosity_level >= 1 && master )
		    printf( "Integration stopped: maximum temperature change exceeded.\n" );
	    }
	} // if G.halt_on_large_flow_change

	// [todo] Dan, I removed the remote control from this section of code
	//        bacause it seemed simply to turn off the following code.
	//        I guess that it's unfinished work.
#	if 0
	if ( G.radiation ) {
	    if ( check_radiation_scaling() ) {
	    	finished_time_stepping = true;
	    	if ( G.verbosity_level >= 1 && master )
		    printf( "Integration stopped: radiation source term needs updating.\n" );
	    }
	}
#       endif
	    	
    } // end while
    //----------------------------------------------------------------
    //                Bottom of main time-stepping loop
    //----------------------------------------------------------------

 conclusion:
    dt_record.clear();
    return status_flag;
} // end integrate_in_time()


int finalize_simulation( void )
{
    global_data &G = *get_global_data_ptr();
    string filename, commandstring, foldername, jbstring, jsstring;
    char tindxcstr[10], jbcstr[10];

    // Write out the final solution only if it has NOT just been written
    // as part of the main time-stepping loop.
    if ( !output_just_written ) {
	++output_counter;
	if ( master ) {
	    fprintf( G.timestampfile, "%04d %e %e\n", static_cast<int>(output_counter),
		     G.sim_time, G.dt_global );
	    fflush( G.timestampfile );
	}
	sprintf( tindxcstr, "t%04d", static_cast<int>(output_counter) ); // C string
	write_solution_data(tindxcstr);

	if ( G.conjugate_ht_active && master ) {
	    write_solution(*(G.wm), G.sim_time, static_cast<int>(output_counter));
	}

	output_just_written = true;
	G.t_plot += G.dt_plot;
    }
    // For the history files, we don't want to double-up on solution data.
    if ( !history_just_written ) {
	for ( Block *bdp : G.my_blocks ) {
	    sprintf( jbcstr, ".b%04d", static_cast<int>(bdp->id) ); jbstring = jbcstr;
	    filename = "hist/"+G.base_file_name+".hist"+jbstring;
            bdp->write_history( filename, G.sim_time );
	    for ( int iface: bdp->transient_profile_faces ) {
		filename = "./" + G.base_file_name + ".blk"+jbstring+"."+get_face_name(iface)+".profile";
		bdp->write_profile(filename, iface, G.sim_time, false);
	    }
	}
        history_just_written = true;
    }
    if ( G.verbosity_level >= 1 && master )
	printf( "\nTotal number of steps = %d\n", static_cast<int>(G.step) );

    filename = G.base_file_name; filename += ".finish";
    if ( master ) {
	write_finishing_data( &G, filename );
	fclose( G.timestampfile );
    }
#   ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#   endif

    for ( Block *bdp : G.my_blocks ) bdp->array_cleanup(G.dimensions); 
#   ifdef _MPI
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD); // just to reduce the jumble in stdout
    if ( delete_send_and_receive_buffers() != SUCCESS ) exit( MEMORY_ERROR );
#   endif
    return SUCCESS;
} // end finalize_simulation()

//------------------------------------------------------------------------

int gasdynamic_explicit_increment_with_fixed_grid(double dt)
// Time level of grid stays at 0.
// 2013-04-07 also updated G.sim_time
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int step_status_flag;
    using std::swap;
    double t0 = G.sim_time;
    // Set the time-step coefficients for the stages of the update scheme.
    double c2 = 1.0; // default for PC_UPDATE
    double c3 = 1.0; // default for PC_UPDATE
    switch ( get_gasdynamic_update_scheme() ) {
    case EULER_UPDATE:
    case PC_UPDATE: c2 = 1.0; c3 = 1.0; break;
    case MIDPOINT_UPDATE: c2 = 0.5; c3 = 1.0; break;
    case CLASSIC_RK3_UPDATE: c2 = 0.5; c3 = 1.0; break;
    case TVD_RK3_UPDATE: c2 = 1.0; c3 = 0.5; break;
    case DENMAN_RK3_UPDATE: c2 = 1.0; c3 = 0.5; break; 
    // FIX-ME: Check Andrew Denman's coefficients.
    default: 
	throw std::runtime_error("gasdynamic_inviscid_increment_with_fixed_grid(): "
				 "unknown update scheme.");
    }
    int attempt_number = 0;
    do {
	//  Preparation for the predictor-stage of inviscid gas-dynamic flow update.
	++attempt_number;
	step_status_flag = 0;
	for ( Block *bdp : G.my_blocks ) {
	    if ( !bdp->active ) continue;
	    bdp->clear_fluxes_of_conserved_quantities(G.dimensions);
	    for ( FV_Cell *cp: bdp->active_cells ) cp->clear_source_vector();
	}
#       ifdef _MPI
        // Before we try to exchange data, everyone's internal data should be up-to-date.
	MPI_Barrier( MPI_COMM_WORLD );
	// Now, it's safe to do the exchange for full-face connections.
	mpi_exchange_boundary_data(COPY_FLOW_STATE, 0);
	copy_mapped_cell_data_via_mpi(COPY_FLOW_STATE, 0);
#       else
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) exchange_shared_boundary_data(bdp->id, COPY_FLOW_STATE, 0);
	}
	copy_mapped_cell_data_via_shmem(COPY_FLOW_STATE, 0);
#       endif
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) {
		apply_convective_bc( *bdp, G.sim_time, G.dimensions );
		// We've put this detector step here because it needs the ghost-cell data
		// to be current, as it should be just after a call to apply_convective_bc().
		if ( get_flux_calculator() == FLUX_ADAPTIVE )
		    bdp->detect_shock_points(G.dimensions);
	    }
	}
	// Non-local radiation transport needs to be performed a-priori for parallelization.
	// Note that Q_rE_rad is not re-evaluated for subsequent stages of the update.
	if ( G.radiation ) perform_radiation_transport();
	// Make a copy of Q_rE_rad so that we can reinstate it at each later stage of the update.
	for ( Block *bdp : G.my_blocks ) {
	    if ( !bdp->active ) continue;
	    for ( FV_Cell *cp: bdp->active_cells ) cp->Q_rE_rad_save = cp->Q_rE_rad;
	}	

	// First-stage of gas-dynamic update.
	for ( Block *bdp : G.my_blocks ) {
	    G.t_level = 0;
	    if ( !bdp->active ) continue;
	    bdp->inviscid_flux(G.dimensions);
	    if ( G.viscous && !G.separate_update_for_viscous_terms ) {	    
		apply_viscous_bc(*bdp, G.sim_time, G.dimensions);
		// if ( G.turbulence_model == TM_K_OMEGA ) apply_menter_boundary_correction(*bdp, 0);
		if ( G.dimensions == 2 ) viscous_derivatives_2D(bdp, 0); else viscous_derivatives_3D(bdp, 0); 
		estimate_turbulence_viscosity(&G, bdp);
		if ( G.dimensions == 2 ) viscous_flux_2D(bdp); else viscous_flux_3D(bdp); 
	    } // end if ( G.viscous )
	    for ( FV_Cell *cp: bdp->active_cells ) {
		cp->add_inviscid_source_vector(0, bdp->omegaz);
		if ( G.udf_source_vector_flag == 1 )
		    add_udf_source_vector_for_cell(cp, 0, G.sim_time);
		if ( G.viscous && !G.separate_update_for_viscous_terms ) {
		    cp->add_viscous_source_vector(with_k_omega && !G.separate_update_for_k_omega_source);
		}
		cp->time_derivatives(0, 0, G.dimensions, with_k_omega);
		bool force_euler = false;
		cp->stage_1_update_for_flow_on_fixed_grid(G.dt_global, force_euler, with_k_omega);
		cp->decode_conserved(0, 1, bdp->omegaz, with_k_omega);
	    } // end for *cp
	} // end of for jb...

	if ( number_of_stages_for_update_scheme(get_gasdynamic_update_scheme()) >= 2 ) {
	    // Preparation for second-stage of gas-dynamic update.
	    G.sim_time = t0 + c2*dt;
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		bdp->clear_fluxes_of_conserved_quantities(G.dimensions);
		for ( FV_Cell *cp: bdp->active_cells ) cp->clear_source_vector();
	    }
#           ifdef _MPI
	    // Before we try to exchange data, everyone's internal data should be up-to-date.
	    MPI_Barrier( MPI_COMM_WORLD );
	    // Now, it's safe to do the exchange for full-face connections.
	    mpi_exchange_boundary_data(COPY_FLOW_STATE, 0);
	    copy_mapped_cell_data_via_mpi(COPY_FLOW_STATE, 0);
#           else
	    for ( Block *bdp : G.my_blocks ) {
		if ( bdp->active )
		    exchange_shared_boundary_data(bdp->id, COPY_FLOW_STATE, 0);
	    }
	    copy_mapped_cell_data_via_shmem(COPY_FLOW_STATE, 0);
#           endif
	    // Second stage of gas-dynamic update.
	    for ( Block *bdp : G.my_blocks ) {
		G.t_level = 1;
		if ( !bdp->active ) continue;
		apply_convective_bc(*bdp, G.sim_time, G.dimensions);
		bdp->inviscid_flux(G.dimensions);
		if ( G.viscous && !G.separate_update_for_viscous_terms ) {
		    apply_viscous_bc(*bdp, G.sim_time, G.dimensions);
		    // if ( G.turbulence_model == TM_K_OMEGA ) apply_menter_boundary_correction(*bdp, 1);
		    if ( G.dimensions == 2 ) viscous_derivatives_2D(bdp, 0); else viscous_derivatives_3D(bdp, 0); 
		    estimate_turbulence_viscosity(&G, bdp);
		    if ( G.dimensions == 2 ) viscous_flux_2D(bdp); else viscous_flux_3D(bdp); 
		} // end if ( G.viscous )
		for ( FV_Cell *cp: bdp->active_cells ) {
		    // Radiation transport was calculated once before staged update; recover saved value.
		    cp->Q_rE_rad = cp->Q_rE_rad_save; 
		    cp->add_inviscid_source_vector(0, bdp->omegaz);
		    if ( G.udf_source_vector_flag == 1 )
			add_udf_source_vector_for_cell(cp, 0, G.sim_time);
		    if ( G.viscous && !G.separate_update_for_viscous_terms ) {
			cp->add_viscous_source_vector(with_k_omega && !G.separate_update_for_k_omega_source);
		    }
		    cp->time_derivatives(0, 1, G.dimensions, with_k_omega);
		    cp->stage_2_update_for_flow_on_fixed_grid(G.dt_global, with_k_omega);
		    cp->decode_conserved(0, 2, bdp->omegaz, with_k_omega);
		} // end for ( *cp
	    } // end for jb loop
	} // end if ( number_of_stages_for_update_scheme() >= 2 

	if ( number_of_stages_for_update_scheme(get_gasdynamic_update_scheme()) >= 3 ) {
	    // Preparation for third stage of gasdynamic update.
	    G.sim_time = t0 + c3*dt;
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		bdp->clear_fluxes_of_conserved_quantities(G.dimensions);
		for ( FV_Cell *cp: bdp->active_cells ) cp->clear_source_vector();
	    }
#           ifdef _MPI
	    // Before we try to exchange data, everyone's internal data should be up-to-date.
	    MPI_Barrier( MPI_COMM_WORLD );
	    // Now, it's safe to do the exchange for full-face connections.
	    mpi_exchange_boundary_data(COPY_FLOW_STATE, 0);
	    copy_mapped_cell_data_via_mpi(COPY_FLOW_STATE, 0);
#           else
	    for ( Block *bdp : G.my_blocks ) {
		if ( bdp->active )
		    exchange_shared_boundary_data(bdp->id, COPY_FLOW_STATE, 0);
	    }
	    copy_mapped_cell_data_via_shmem(COPY_FLOW_STATE, 0);
#           endif
	    for ( Block *bdp : G.my_blocks ) {
		G.t_level = 2;
		if ( !bdp->active ) continue;
		apply_convective_bc(*bdp, G.sim_time, G.dimensions);
		bdp->inviscid_flux( G.dimensions );
		if ( G.viscous && !G.separate_update_for_viscous_terms ) {
		    apply_viscous_bc(*bdp, G.sim_time, G.dimensions);
		    // if ( G.turbulence_model == TM_K_OMEGA ) apply_menter_boundary_correction(*bdp, 2);
		    if ( G.dimensions == 2 ) viscous_derivatives_2D(bdp, 0); else viscous_derivatives_3D(bdp, 0); 
		    estimate_turbulence_viscosity(&G, bdp);
		    if ( G.dimensions == 2 ) viscous_flux_2D(bdp); else viscous_flux_3D(bdp); 
		} // end if ( G.viscous )
		for ( FV_Cell *cp: bdp->active_cells ) {
		    // Radiation transport was calculated once before staged update; recover saved value.
		    cp->Q_rE_rad = cp->Q_rE_rad_save; 
		    cp->add_inviscid_source_vector(0, bdp->omegaz);
		    if ( G.udf_source_vector_flag == 1 )
			add_udf_source_vector_for_cell(cp, 0, G.sim_time);
		    if ( G.viscous && !G.separate_update_for_viscous_terms ) {
			cp->add_viscous_source_vector(with_k_omega && !G.separate_update_for_k_omega_source);
		    }
		    cp->time_derivatives(0, 2, G.dimensions, with_k_omega);
		    cp->stage_3_update_for_flow_on_fixed_grid(G.dt_global, with_k_omega);
		    cp->decode_conserved(0, 3, bdp->omegaz, with_k_omega);
		} // for *cp
	    } // end for *bdp
	} // end if ( number_of_stages_for_update_scheme() >= 3
   
	// 2d. Check the record of bad cells and if any cells are bad, 
	//     fail this attempt at taking a step,
	//     set everything back to the initial state and
	//     reduce the time step for the next attempt
	int most_bad_cells = do_bad_cell_count(0);
	if ( !G.adjust_invalid_cell_data && most_bad_cells > 0 ) {
	    step_status_flag = 1;
	}
#       ifdef _MPI
	MPI_Allreduce(MPI_IN_PLACE, &(step_status_flag), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#       endif
	if ( step_status_flag != 0 ) {
	    G.dt_global = G.dt_reduction_factor * G.dt_global;
	    printf("Attempt %d failed: reducing dt to %e.\n", attempt_number, G.dt_global);
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		for ( FV_Cell *cp: bdp->active_cells ) {
		    cp->decode_conserved(0, 0, bdp->omegaz, with_k_omega);
		}
	    } // end for *bdp
	} // end if step_status_flag

    } while (attempt_number < 3 && step_status_flag == 1);

    // Get the end conserved data into U[0] for next step.
    size_t end_indx = 2;
    switch ( get_gasdynamic_update_scheme() ) {
    case EULER_UPDATE: end_indx = 1; break;
    case PC_UPDATE: 
    case MIDPOINT_UPDATE: end_indx = 2; break;
    case TVD_RK3_UPDATE:
    case CLASSIC_RK3_UPDATE:
    case DENMAN_RK3_UPDATE: end_indx = 3; break;
    default:
	throw std::runtime_error("gasdynamic_inviscid_increment_with_fixed_grid(): "
				 "unknown update scheme.");
    }
    for ( Block *bdp : G.my_blocks ) {
	if ( !bdp->active ) continue;
	for ( FV_Cell *cp: bdp->active_cells ) {
	    swap(cp->U[0], cp->U[end_indx]);
	}
    } // end for *bdp
    G.sim_time = t0 + dt;
    return step_status_flag;
} // end gasdynamic_inviscid_increment_with_fixed_grid()


int gasdynamic_increment_with_moving_grid(double dt)
// We have implemented only the simplest consistent two-stage update scheme. 
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int step_status_flag;
    using std::swap;

    // FIX-ME moving grid: this is a work/refactoring in progress... PJ
    // 25-Mar-2103 except for superficial changes, it is the same as 
    // Andrew's implementation.
    // 31-Mar-2013 let's get serious about refactoring.
    // 02-Apr-2013 have finally arrived at working on this section
    // but would like Andrew to check and fix my indexing of the
    // grid and flow time-levels in the following code.

    double t0 = G.sim_time;
    int attempt_number = 0;
    do {
	//  Preparation for the first-stage of inviscid gas-dynamic flow update.
	++attempt_number;
	step_status_flag = 0;		
        
	// predict new vertex positions
	for ( Block *bdp : G.my_blocks ) {
            bdp->predict_vertex_positions(G.dimensions, G.dt_global);
	    bdp->compute_primary_cell_geometric_data(G.dimensions, 1); // for start of next stage
	    bdp->set_gcl_interface_properties(G.dimensions, 1, G.dt_global);		
	}	
	for ( Block *bdp : G.my_blocks ) {
	    if ( !bdp->active ) continue;
	    bdp->clear_fluxes_of_conserved_quantities(G.dimensions);
	    for ( FV_Cell *cp: bdp->active_cells ) cp->clear_source_vector();
	}
#       ifdef _MPI
        // Before we try to exchange data, everyone's data should be up-to-date.
	MPI_Barrier( MPI_COMM_WORLD );
	mpi_exchange_boundary_data(COPY_FLOW_STATE, 1);
	// Separate exchange of interface data to avoid message collisions.
#       else
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) {
		exchange_shared_boundary_data(bdp->id, COPY_FLOW_STATE, 1);
	    }
	}
#       endif
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) {
		apply_convective_bc( *bdp, G.sim_time, G.dimensions );
		// We've put this detector step here because it needs the ghost-cell data
		// to be current, as it should be just after a call to apply_convective_bc().
		if ( get_flux_calculator() == FLUX_ADAPTIVE )
		    bdp->detect_shock_points( G.dimensions );
	    }
	}
#       ifdef _MPI
        // Before we try to exchange data, everyone's data should be up-to-date.
	// Separate exchange of interface data to avoid message collisions.
	MPI_Barrier( MPI_COMM_WORLD );
	mpi_exchange_boundary_data(COPY_INTERFACE_DATA, 1);
#       else
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) {
		exchange_shared_boundary_data(bdp->id, COPY_INTERFACE_DATA, 1);
	    }
	}
#       endif
	// Non-local radiation transport needs to be performed a-priori for parallelization.
	// Note that Q_rE_rad is not re-evaluated for subsequent stages of the update.
	if ( G.radiation ) perform_radiation_transport();
	// Make a copy of Q_rE_rad so that we can reinstate it at each later stage of the update.
	for ( Block *bdp : G.my_blocks ) {
	    if ( !bdp->active ) continue;
	    for ( FV_Cell *cp: bdp->active_cells ) cp->Q_rE_rad_save = cp->Q_rE_rad;
	}					
	// First-stage of gas-dynamic update.
	for ( Block *bdp : G.my_blocks ) {
	    if ( !bdp->active ) continue;
	    bdp->inviscid_flux( G.dimensions );
            if ( G.viscous && !G.separate_update_for_viscous_terms ) {
		apply_viscous_bc(*bdp, G.sim_time, G.dimensions);
		// if ( G.turbulence_model == TM_K_OMEGA ) apply_menter_boundary_correction(*bdp, 0);
		if ( G.dimensions == 2 ) viscous_derivatives_2D(bdp, 1); else viscous_derivatives_3D(bdp, 1); 
		    estimate_turbulence_viscosity(&G, bdp);
		if ( G.dimensions == 2 ) viscous_flux_2D(bdp); else viscous_flux_3D(bdp); 
	    } // end if ( G.viscous )	    
	    for ( FV_Cell *cp: bdp->active_cells ) {
		// Radiation transport was calculated once before staged update; recover saved value.
		cp->Q_rE_rad = cp->Q_rE_rad_save; 	    
		cp->add_inviscid_source_vector(1, bdp->omegaz);
		if ( G.udf_source_vector_flag == 1 ) add_udf_source_vector_for_cell(cp, 1, G.sim_time);
		if ( G.viscous && !G.separate_update_for_viscous_terms ) {
		    cp->add_viscous_source_vector(with_k_omega && !G.separate_update_for_k_omega_source);
		}		
		cp->time_derivatives(1, 0, G.dimensions, with_k_omega);
		cp->stage_1_update_for_flow_on_moving_grid(G.dt_global, with_k_omega);
		cp->decode_conserved(1, 1, bdp->omegaz, with_k_omega);
	    } // end for *cp
	} // end of for *bdp...

	// Preparation for second-stage of gas-dynamic update.
	for ( Block *bdp : G.my_blocks ) {
	    bdp->correct_vertex_positions(G.dimensions, G.dt_global);
	    bdp->compute_primary_cell_geometric_data(G.dimensions, 2);
	    bdp->set_gcl_interface_properties(G.dimensions, 2, G.dt_global);	
	}
	for ( Block *bdp : G.my_blocks ) {
	    if ( !bdp->active ) continue;
	    bdp->clear_fluxes_of_conserved_quantities(G.dimensions);
	    for ( FV_Cell *cp: bdp->active_cells ) cp->clear_source_vector();
	}
#       ifdef _MPI
	MPI_Barrier( MPI_COMM_WORLD );
	mpi_exchange_boundary_data(COPY_FLOW_STATE, 2);
	// Separate exchange of interface data to avoid message collisions.
#       else
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) {
		exchange_shared_boundary_data(bdp->id, COPY_FLOW_STATE, 2);
	    }
	}
#       endif
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) apply_convective_bc( *bdp, G.sim_time, G.dimensions );
	}
#       ifdef _MPI
        // Before we try to exchange data, everyone's data should be up-to-date.
	// Separate exchange of interface data to avoid message collisions.
	MPI_Barrier( MPI_COMM_WORLD );
	mpi_exchange_boundary_data(COPY_INTERFACE_DATA, 2);
#       else
	for ( Block *bdp : G.my_blocks ) {
	    if ( bdp->active ) {
		exchange_shared_boundary_data(bdp->id, COPY_INTERFACE_DATA, 2);
	    }
	}
#       endif
	// Second-stage of gas-dynamic update.
	for ( Block *bdp : G.my_blocks ) {
	    if ( !bdp->active ) continue;
	    bdp->inviscid_flux( G.dimensions );
            if ( G.viscous && !G.separate_update_for_viscous_terms ) {
		apply_viscous_bc(*bdp, G.sim_time, G.dimensions);	
		// if ( G.turbulence_model == TM_K_OMEGA ) apply_menter_boundary_correction(*bdp, 1);
		if ( G.dimensions == 2 ) viscous_derivatives_2D(bdp, 2); else viscous_derivatives_3D(bdp, 2); 
		    estimate_turbulence_viscosity(&G, bdp);
		if ( G.dimensions == 2 ) viscous_flux_2D(bdp); else viscous_flux_3D(bdp); 
	    } // end if ( G.viscous )	    
	    for ( FV_Cell *cp: bdp->active_cells ) {
		// Radiation transport was calculated once before staged update; recover saved value.
		cp->Q_rE_rad = cp->Q_rE_rad_save; 
		cp->add_inviscid_source_vector(1, bdp->omegaz);
		if ( G.udf_source_vector_flag == 1 ) add_udf_source_vector_for_cell(cp, 2, G.sim_time);
		if ( G.viscous && !G.separate_update_for_viscous_terms ) {
		    cp->add_viscous_source_vector(with_k_omega && !G.separate_update_for_k_omega_source);
		}		
		cp->time_derivatives(2, 1, G.dimensions, with_k_omega);
		cp->stage_2_update_for_flow_on_moving_grid(G.dt_global, with_k_omega);
		cp->decode_conserved(2, 2, bdp->omegaz, with_k_omega);
	    } // end for *cp
	} // end for jb loop
	   
	// 2d. Check the record of bad cells and if any cells are bad, 
	//     fail this attempt at taking a step,
	//     set everything back to the initial state and
	//     reduce the time step for the next attempt
	int most_bad_cells = do_bad_cell_count(2);
	if ( !G.adjust_invalid_cell_data && most_bad_cells > 0 ) {
	    step_status_flag = 1;
	}
#       ifdef _MPI
	MPI_Allreduce(MPI_IN_PLACE, &(step_status_flag), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#       endif
	if ( step_status_flag != 0 ) {
	    G.dt_global = G.dt_reduction_factor * G.dt_global;
	    printf("Attempt %d failed: reducing dt to %e.\n", attempt_number, G.dt_global);
	    for ( Block *bdp : G.my_blocks ) {
		if ( !bdp->active ) continue;
		for ( FV_Cell *cp: bdp->active_cells ) {
		    cp->decode_conserved(0, 0, bdp->omegaz, with_k_omega);
		}
	    } // end for *bdp
	} // end if step_status_flag
	
    } while (attempt_number < 3 && step_status_flag == 1);

    for ( Block *bdp : G.my_blocks ) {
	if ( !bdp->active ) continue;
	for ( FV_Cell *cp: bdp->active_cells ) {
	    swap(cp->U[0], cp->U[2]);
	    cp->copy_grid_level_to_level(2, 0);
	}
    }
    
    G.sim_time = t0 + dt;
    return step_status_flag;
} // end gasdynamic_inviscid_increment_with_moving_grid()


int gasdynamic_separate_explicit_viscous_increment()
// A simple first-order Euler step, of just the viscous-terms contribution.
// Note that it starts from scratch, following the inviscid update.
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    bool force_euler = true;
    using std::swap;
    for ( Block *bdp : G.my_blocks ) {
        if ( !bdp->active ) continue;
	bdp->clear_fluxes_of_conserved_quantities(G.dimensions);
	for ( FV_Cell *cp: bdp->active_cells ) cp->clear_source_vector();
	apply_viscous_bc(*bdp, G.sim_time, G.dimensions);
	// if ( G.turbulence_model == TM_K_OMEGA ) apply_menter_boundary_correction(*bdp, 0);
	if ( G.dimensions == 2 ) viscous_derivatives_2D(bdp, 0); else viscous_derivatives_3D(bdp, 0); 
	estimate_turbulence_viscosity(&G, bdp);
	if ( G.dimensions == 2 ) viscous_flux_2D(bdp); else viscous_flux_3D(bdp); 
	for ( FV_Cell *cp: bdp->active_cells ) {
	    cp->add_viscous_source_vector(with_k_omega);
	    cp->time_derivatives(0, 0, G.dimensions, with_k_omega);
	    cp->stage_1_update_for_flow_on_fixed_grid(G.dt_global, force_euler, with_k_omega);
	    swap(cp->U[0], cp->U[1]);
	    cp->decode_conserved(0, 0, bdp->omegaz, with_k_omega);
	} // end for *cp
    } // end for *bdp...
    return SUCCESS;
} // int gasdynamic_separate_explicit_viscous_increment()


int do_bad_cell_count(size_t gtl)
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    int most_bad_cells = 0;
    for ( Block *bdp : G.my_blocks ) {
	size_t bad_cell_count = bdp->count_invalid_cells(G.dimensions, gtl);
	if ( static_cast<int>(bad_cell_count) > most_bad_cells ) 
	    most_bad_cells = static_cast<int>(bad_cell_count); 
	if ( bad_cell_count > G.max_invalid_cells ) {
	    printf( "   Too many bad cells (i.e. %d > %d) in block[%d].\n", 
		    static_cast<int>(bad_cell_count), 
		    static_cast<int>(G.max_invalid_cells), 
		    static_cast<int>(bdp->id) );
	}
    } // end for *bdp
#   ifdef _MPI
    MPI_Allreduce(MPI_IN_PLACE, &most_bad_cells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#   endif
    return most_bad_cells;
} // end do_bad_cell_count()


/// \brief Write out simulation data at end of simulation, such as: no. steps, final time.
int write_finishing_data( global_data *G, std::string filename )
{
    FILE *fp;
    if ((fp = fopen(filename.c_str(), "a")) == NULL) {
	cerr << "Could not open " << filename << "; BAILING OUT" << endl;
	exit( FILE_ERROR );
    }
    fprintf(fp, "[simulation_end]\n");
    fprintf(fp, "final_time = %.12e \n", G->sim_time);
    fprintf(fp, "dt = %.12e \n", G->dt_allow);
    fprintf(fp, "no_steps = %d\n", static_cast<int>(G->step));
    fclose( fp );
    return SUCCESS;

}

/// \brief Check that the maximum delta(f_rad) / f_rad_org value is not exceedingly large
int check_radiation_scaling( void )
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    double global_f_max = 0.0;
    for ( Block *bdp : G.my_blocks ) {
    	double block_f_max = 0.0;
	for ( FV_Cell *cellp : bdp->active_cells ) {
	    double f_max = cellp->rad_scaling_ratio();
	    if ( f_max > block_f_max) block_f_max = f_max;
	}
	if ( block_f_max > global_f_max ) global_f_max = block_f_max;
    }
#   ifdef _MPI
    MPI_Allreduce(MPI_IN_PLACE, &global_f_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#   endif
    // printf("Global radiation f_max = %f\n", global_f_max);
    if ( global_f_max > 1000000.0 )
    	return FAILURE;
    else
    	return SUCCESS;
}

//---------------------------------------------------------------------------

int radiation_calculation()
{
    global_data &G = *get_global_data_ptr();  // set up a reference
    // Apply viscous THEN inviscid boundary conditions to match environment for
    // radiation calculation in gasdynamic_inviscid_increment()
    for ( Block *bdp : G.my_blocks ) {
	if ( !bdp->active ) continue;
	if ( G.viscous ) apply_viscous_bc( *bdp, G.sim_time, G.dimensions ); 
	apply_convective_bc( *bdp, G.sim_time, G.dimensions );
    }
    if ( G.radiation ) perform_radiation_transport();
    return SUCCESS;
} // end radiation_calculation()

void perform_radiation_transport()
{
    RadiationTransportModel * rtm = get_radiation_transport_model_ptr();
    global_data &G = *get_global_data_ptr();
    Block * bdp;

#   ifdef E3RAD
    // Always recompute radiation with e3rad - parallelised via openmp
    rtm->compute_Q_rad_for_flowfield();
    // store the radiation scaling parameters for each cell
    size_t jb;
#   ifdef _OPENMP
#   pragma omp parallel for private(jb) schedule(runtime)
#   endif
    for ( jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	if ( !bdp->active ) continue;
	for ( FV_Cell *cp: bdp->active_cells ) cp->store_rad_scaling_params();
    }
#   else
    // e3shared.exe and e3mpi.exe version - currently no parallelisation
    // Determine if a scaled or complete radiation call is required
    int ruf = G.radiation_update_frequency;
    if ( (ruf == 0) || ((G.step/ruf)*ruf != G.step) ) {
	// rescale
	for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( !bdp->active ) continue;
//	    for ( FV_Cell *cp: bdp->active_cells ) cp->rescale_Q_rE_rad();
	}
    } else if ( G.radiation_scaling ){
	// recompute
	rtm->compute_Q_rad_for_flowfield();
	// store the radiation scaling parameters for each cell
	for ( size_t jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( !bdp->active ) continue;
	    for ( FV_Cell *cp: bdp->active_cells ) cp->store_rad_scaling_params();
	}
    }
#   endif

    return;
} // end perform_radiation_transport()
