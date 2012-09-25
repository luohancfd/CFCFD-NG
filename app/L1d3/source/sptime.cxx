/** \file sptime.cxx
 * \ingroup l1d3
 * \brief Postprocessor for l1d.cxx: this one produces data for X-T diagrams.
 *
 * This program is used to read the flow
 * solution data for the One-Dimensional Lagrangian
 * code and then generate a generic plotting file
 * containing space-time data for a particular variable.
 *
 * \author PA Jacobs
 * \author Jan Martinez-Schramm (TECPLOT output)
 *
 * \version 1.0 --  30-Dec-91
 * \version 6.0 --  05-Jun-00, Allowed for adaptivity.
 * \version 6.1 --  17-Nov-02, Jan Martinez-Schramm added TECPLOT output
 * \version 24-Jul-06, C++ port.
 */

/*-----------------------------------------------------------------*/

#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <valarray>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l_kernel.hh"
#include "l1d.hh"
#include "l_io.hh"
#include "l_misc.hh"

/*-----------------------------------------------------------------*/

int main(int argc, char **argv)
{
    struct simulation_data SD;
    int js, jp, jd;
    vector<int> nnx_initial;
    std::vector<slug_data> A;               /* several gas slugs        */
    std::vector<piston_data> Pist;          /* room for several pistons */
    std::vector<diaphragm_data> Diaph;      /* diaphragms            */

    double tstart, tstop;
    int i, max_sol, tecplot_format;
    char oname[40], pname[40], iname[40];
    char base_file_name[32], name_tag[10];
    FILE *infile,            /* flow solution file           */
        *outfile;           /* plot file                    */
    char var_name[132],var_name_tec[132];

#   define  BLOCKS  10
    double **xarray[BLOCKS], **tarray[BLOCKS], **varray[BLOCKS];
    int nt, nt_read, nt_write, ix, nx, nxtot, option;
    int takelog;
    int found, ixmin, ixmax;
    double dx, xloc, xL, xR;
    int command_line_error;
    int echo_input;

    /*
     * ----------
     * INITIALIZE
     * ----------
     */
    printf("\n---------------------------");
    printf("\nLagrangian 1D Postprocessor");
    printf("\nSpace - Time Plots         ");
    printf("\n---------------------------");
    printf("\n\n");

    strcpy(base_file_name, "default");
    max_sol = 500;
    tstart = 0.0;
    tstop = 0.0;
    takelog = 0;
    tecplot_format = 0;
#   define  SELECT_RHO  1
#   define  SELECT_U    2
#   define  SELECT_E    3
#   define  SELECT_P    4
#   define  SELECT_A    5
#   define  SELECT_T    6
#   define  SELECT_TAU  7
#   define  SELECT_Q    8
#   define  SELECT_S    9
    option = SELECT_P;
    echo_input = 0;

    /*
     * Decode command line arguments.
     */
    command_line_error = 0;
    if (argc < 2) {
        /* no arguments supplied */
        command_line_error = 1;
        goto usage;
    }   /* end if */

    i = 1;
    while (i < argc) {
        /* process the next command-line argument */

        if (strcmp(argv[i], "-f") == 0) {
            /* Set the base file name. */
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
        } else if (strcmp(argv[i], "-tstop") == 0) {
            /* Set the end solution time. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%lf", &tstop);
            i++;
            printf("Setting tstop = %e\n", tstop);
        } else if (strcmp(argv[i], "-tstart") == 0) {
            /* Set the start solution time. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%lf", &tstart);
            i++;
            printf("Setting tstart = %e\n", tstart);
        } else if (strcmp(argv[i], "-maxsol") == 0) {
            /* Set the maximum number of solutions. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%d", &max_sol);
            i++;
            printf("Setting max_sol = %d\n", max_sol);
        } else if (strcmp(argv[i], "-tecplot") == 0) {
            tecplot_format = 1;
            i++;
            printf("Will write output in TECPLOT format.\n");
        } else if (strcmp(argv[i], "-log") == 0) {
            takelog = 1;
            i++;
            printf("Take logarithm.\n");
        } else if (strcmp(argv[i], "-p") == 0) {
            option = SELECT_P;
            i++;
            printf("Pressure, option = %d.\n", option);
        } else if (strcmp(argv[i], "-rho") == 0) {
            option = SELECT_RHO;
            i++;
            printf("Density, option = %d.\n", option);
        } else if (strcmp(argv[i], "-tau") == 0) {
            option = SELECT_TAU;
            i++;
            printf("Density, option = %d.\n", option);
        } else if (strcmp(argv[i], "-u") == 0) {
            option = SELECT_U;
            i++;
            printf("Velocity, option = %d.\n", option);
        } else if (strcmp(argv[i], "-e") == 0) {
            option = SELECT_E;
            i++;
            printf("Internal energy, option = %d.\n", option);
        } else if (strcmp(argv[i], "-T") == 0) {
            option = SELECT_T;
            i++;
            printf("Temperature, option = %d.\n", option);
        } else if (strcmp(argv[i], "-a") == 0) {
            option = SELECT_A;
            i++;
            printf("Sound speed, option = %d.\n", option);
        } else if (strcmp(argv[i], "-q") == 0) {
            option = SELECT_Q;
            i++;
            printf("Heat transfer, option = %d.\n", option);
        } else if (strcmp(argv[i], "-S") == 0) {
            option = SELECT_S;
            i++;
            printf("Entropy, option = %d.\n", option);
        } else if (strcmp(argv[i], "-help") == 0) {
            command_line_error = 1;
            printf("Printing usage message...\n\n");
            goto usage;
        } else {
            printf("Unknown option: %s\n", argv[i]);
            command_line_error = 1;
            goto usage;
        }   /* end if */
    }   /* end while */

    /*
     * If the command line arguments are incorrect,
     * write some advice to the console then finish.
     */
    usage:
    if (command_line_error == 1) {
        printf("Purpose:\n");
        printf("Pick up a solution file <base_file_name>.Ls, select one\n");
        printf("of the flow variables and then either save the x-t data\n");
        printf
            ("for that variable in a format suitable for contour plotting.\n");
        printf
            ("When reformatting, a logarithmic scale can be used, however,\n");
        printf
            ("taking the logarithm of velocity or shear stress is probably\n");
        printf("a bad idea.\n");
        printf("\n");

        printf("Command-line options: (defaults are shown in parentheses)\n");
        printf("-f <base_file_name>       (default)\n");
        printf("-echo                     (no echo)\n");
        printf("-tstart <time>            (0.0)\n");
        printf("-tstop <time>             (0.0)\n");
        printf("-maxsol <n>               (%d)\n", max_sol);
        printf("-log                      (linear)\n");
        printf("-tecplot                  (generic)\n");
        printf("-help                     (print this message)\n");
        printf("To select a variable, pick one of:\n");
        printf("-p -rho -u -e -T -a -q -tau -S   (default is -p)\n");
        printf("\n");
        exit(1);
    }   /* end if command_line_error */

    /*
     * * Read the input parameter file.
     */
    strcpy(pname, base_file_name);
    strcat(pname, ".Lp");
    printf("parameterfile: %s\n", pname);
    ConfigParser parameterdict = ConfigParser(pname);
    L_set_case_parameters(&SD, parameterdict, echo_input);
    A.resize(SD.nslug);
    nnx_initial.resize(SD.nslug);
    Pist.resize(SD.npiston);
    Diaph.resize(SD.ndiaphragm);
    Gas_model *gmodel = get_gas_model_ptr();
    for (jp = 0; jp < SD.npiston; ++jp) {
        set_piston_parameters(&(Pist[jp]), jp, parameterdict, SD.dt_init, echo_input);
        Pist[jp].sim_time = 0.0;
    }   /* end for */
    for (jd = 0; jd < SD.ndiaphragm; ++jd) {
        set_diaphragm_parameters(&(Diaph[jd]), jd, parameterdict, echo_input);
        Diaph[jd].sim_time = 0.0;
    }   /* end for */
    for (js = 0; js < SD.nslug; ++js) {
        L_set_slug_parameters(&(A[js]), js, &SD, parameterdict, echo_input);
        A[js].sim_time = 0.0;
        nnx_initial[js] = A[js].nnx;
    }   /* end for */

    /*
     * Allocate enough memory to accumulate the data.
     */
    for (js = 0; js < SD.nslug; ++js) {
        if ((xarray[js] = (double **) calloc(max_sol,sizeof(double *)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
        if ((tarray[js] = (double **) calloc(max_sol,sizeof(double *)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
        if ((varray[js] = (double **) calloc(max_sol,sizeof(double *)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
        for (nt = 0; nt < max_sol; ++nt) {
            nxtot = nnx_initial[js];
            if ((xarray[js][nt] = (double *) calloc(nxtot,sizeof(double)))
                == NULL) {
                printf("Memory alloc problem.\n");
                exit(-1);
            }   /* end if */
            if ((tarray[js][nt] = (double *) calloc(nxtot,sizeof(double)))
                == NULL) {
                printf("Memory alloc problem.\n");
                exit(-1);
            }   /* end if */
            if ((varray[js][nt] = (double *) calloc(nxtot,sizeof(double)))
                == NULL) {
                printf("Memory alloc problem.\n");
                exit(-1);
            }   /* end if */
        }   /* end for (nt = 0;...  */

    }   /* end for (js = 0;...  */

    /*
     * Read all of the solutions and save the requested data.
     */
    strcpy(iname, base_file_name);
    strcat(iname, ".Ls");
    printf("infile       : %s\n", iname);
    if ((infile = fopen(iname, "r")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", iname);
        exit(1);
    }   /* end if */

    nt_write = 0;
    nt_read  = 0;
    for (nt = 0; nt < max_sol; ++nt) {
        for (jp = 0; jp < SD.npiston; ++jp)
            read_piston_solution(&(Pist[jp]), infile);
        for (jd = 0; jd < SD.ndiaphragm; ++jd)
            read_diaphragm_solution(&(Diaph[jd]), infile);
        for (js = 0; js < SD.nslug; ++js)
            L_read_solution(&(A[js]), infile);
	++nt_read;

	if ( A[0].sim_time < tstart ) {
	    /* Don't add this solution to our collection, continue reading. */
	    printf(".");
	    fflush(stdout);
	    continue;
        }

        for (js = 0; js < SD.nslug; ++js) {
            if (A[js].adaptive == 0) {
                /*
                 * For a slug with constant number of cells, 
                 * use individual cell data.
                 */
                nx = 0;
                for (ix = A[js].ixmin; ix <= A[js].ixmax; ++ix) {
                    xarray[js][nt_write][nx] = 0.5 *
                        (A[js].Cell[ix - 1].x + A[js].Cell[ix].x);
                    tarray[js][nt_write][nx] = A[js].sim_time;
                    if (option == SELECT_RHO) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].gas->rho;
                    } else if (option == SELECT_U) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].u;
                    } else if (option == SELECT_E) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].gas->e[0];
                    } else if (option == SELECT_P) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].gas->p;
                    } else if (option == SELECT_A) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].gas->a;
                    } else if (option == SELECT_T) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].gas->T[0];
                    } else if (option == SELECT_TAU) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].shear_stress;
                    } else if (option == SELECT_Q) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].heat_flux;
                    } else if (option == SELECT_S) {
                        varray[js][nt_write][nx] = A[js].Cell[ix].entropy;
                    } else {
                        printf("Invalid option: start again.\n");
                        exit(-1);
                    }   /* end if */

                    if (takelog == 1) {
                        varray[js][nt_write][nx] = log10(fabs(varray[js][nt_write][nx]));
                    }
                    /* end if */
                    ++nx;
                }   /* end for ix */
            } else {
                /*
                 * For a slug with a variable number of cells,
                 * interpolate the data from the actual cells onto
                 * an evenly distributed number of points.
                 */
		struct L_cell* icell = new(struct L_cell);
		icell->gas = new Gas_data(gmodel);
                ixmin = A[js].ixmin;
                ixmax = A[js].ixmax;
                xL = 0.5 * (A[js].Cell[ixmin - 1].x + A[js].Cell[ixmin].x);
                xR = 0.5 * (A[js].Cell[ixmax - 1].x + A[js].Cell[ixmax].x);
                dx = (xR - xL) / (nnx_initial[js] - 1);
                for (nx = 0; nx < nnx_initial[js]; ++nx) {
                    xloc = xL + dx * nx;
                    found = L_interpolate_cell_data(&(A[js]), xloc, *icell);
                    xarray[js][nt_write][nx] = xloc;
                    tarray[js][nt_write][nx] = A[js].sim_time;
                    if (option == SELECT_RHO) {
                        varray[js][nt_write][nx] = icell->gas->rho;
                    } else if (option == SELECT_U) {
                        varray[js][nt_write][nx] = icell->u;
                    } else if (option == SELECT_E) {
                        varray[js][nt_write][nx] = icell->gas->e[0];
                    } else if (option == SELECT_P) {
                        varray[js][nt_write][nx] = icell->gas->p;
                    } else if (option == SELECT_A) {
                        varray[js][nt_write][nx] = icell->gas->a;
                    } else if (option == SELECT_T) {
                        varray[js][nt_write][nx] = icell->gas->T[0];
                    } else if (option == SELECT_TAU) {
                        varray[js][nt_write][nx] = icell->shear_stress;
                    } else if (option == SELECT_Q) {
                        varray[js][nt_write][nx] = icell->heat_flux;
                    } else if (option == SELECT_S) {
                        varray[js][nt_write][nx] = icell->entropy;
                    } else {
                        printf("Invalid option: start again.\n");
                        exit(-1);
                    }   /* end if */

                    if (takelog == 1) {
                        varray[js][nt_write][nx] = log10(fabs(varray[js][nt_write][nx]));
                    }   /* end if */
                }   /* end for */
            }   /* end if */
        }   /* end for js */

	printf("%%");
	fflush(stdout);
        ++nt_write;

        if (A[0].sim_time >= tstop || nt_write >= max_sol)
            break;
    }   /* end for nt... */

    printf("\n");
    if (infile != NULL)
        fclose(infile);
    printf("Number of solutions read: %d\n", nt_read);
    printf("Final time = %g\n", A[0].sim_time);


    printf("Now, write out the contouring data.\n");

    strcpy(oname, base_file_name);
    if (takelog == 1) {
        strcpy(name_tag, "_log");
        strcat(oname, name_tag);
	strcpy(var_name, "log-");
	strcpy(var_name_tec, "log-");
    } else {
	strcpy(var_name, "");
	strcpy(var_name_tec, "");
    }   /* end if */
    if (option == SELECT_P) {
        strcpy(name_tag, "_p");
	strcat(var_name, "p,Pa");
	strcat(var_name_tec, "p:[Pa]");
    } else if (option == SELECT_RHO) {
        strcpy(name_tag, "_rho");
	strcat(var_name, "rho,kg/m^3");
	strcat(var_name_tec, "rho:[kg/m^3]");
    } else if (option == SELECT_TAU) {
        strcpy(name_tag, "_tau");
	strcat(var_name, "tau0,Pa");
	strcat(var_name_tec, "tau0:[Pa]");
    } else if (option == SELECT_T) {
        strcpy(name_tag, "_T");
	strcat(var_name, "T,K");
	strcat(var_name_tec, "T:[K]");
    } else if (option == SELECT_U) {
        strcpy(name_tag, "_u");
	strcat(var_name, "u,m/s");
	strcat(var_name_tec, "u:[m/s]");
    } else if (option == SELECT_E) {
        strcpy(name_tag, "_e");
	strcat(var_name, "e,J/kg");
	strcat(var_name_tec, "e:[J/kg]");
    } else if (option == SELECT_Q) {
        strcpy(name_tag, "_q");
	strcat(var_name, "q,W/m^2");
	strcat(var_name_tec, "q:[W/m^2]");
    } else if (option == SELECT_A) {
        strcpy(name_tag, "_a");
	strcat(var_name, "a,m/s");
	strcat(var_name_tec, "a:[m/s]");
    } else if (option == SELECT_S) {
        strcpy(name_tag, "_S");
	strcat(var_name, "S2,J/kg/K");
	strcat(var_name_tec, "S2:[J/kg/K]");
    }   /* end if */
    strcat(oname, name_tag);
    if ( tecplot_format ) {
	strcat(oname, ".plt");
    } else {
	strcat(oname, ".gen");
    }
    printf("Output file: %s\n", oname);

    if ((outfile = fopen(oname, "w")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", oname);
        exit(-1);
    }   /* end if */
    if ( tecplot_format ) {
	fprintf(outfile, "TITLE=\"L1D:Space-time plot\"\n");
	fprintf(outfile, "VARIABLES= \"x:[m]\", \"t:[s]\", \"%s\"\n", 
		var_name_tec );
    } else {
	fprintf(outfile, "L1D:space-time-plot\n");
	fprintf(outfile, "%d   <== number of columns\n", 3);
	fprintf(outfile, "x,m\n");
	fprintf(outfile, "t,s\n");
	fprintf(outfile, "%s\n", var_name );
	fprintf(outfile, "%d   <== number of slugs\n", SD.nslug);
    }

    for (js = 0; js < SD.nslug; ++js) {
	if ( tecplot_format ) {
	    fprintf(outfile, "ZONE T=\"SLUG-%0d\"\n",js);
	    fprintf(outfile, "I=%d, J=%d, K=1, F=POINT\n", 
		    nt_write, nnx_initial[js]);
	} else {
	    fprintf(outfile, "%d %d  <== ny, nx\n", 
		    nt_write, nnx_initial[js]);
	}
        for (nx = 0; nx < nnx_initial[js]; ++nx) {
            for (nt = 0; nt < nt_write; ++nt) {
                fprintf(outfile, "%e %e %e\n", xarray[js][nt][nx],
                        tarray[js][nt][nx], varray[js][nt][nx]);
            }   /* end for nt */
        }   /* end for nx */
    }   /* end for (js = 0;... */

    if (outfile != NULL) fclose(outfile);
    printf("Number of solutions written: %d\n", nt_write);

    return 0;
}   /* end function main */

/*================ end of sptime.c =================*/
