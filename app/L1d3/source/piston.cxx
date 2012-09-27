/** \file piston.cxx
 * \ingroup l1d2
 * \brief Postprocessor for l1d.cxx -- extract piston data.
 *
 * This program is used to read the flow
 * solution data for the One-Dimensional Lagrangian
 * code and then generate a generic plotting file
 * containing space-time data for a particular piston.
 *
 * \author PA Jacobs
 *
 * \version 1.0 -- 03-Feb-92, adapted from sptime.c
 * \version 24-Jul-06, C++ port.
 */

/*-----------------------------------------------------------------*/

#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../../../lib/util/source/config_parser.hh"
#include "l1d.hh"
#include "l_io.hh"
#include "l_diaph.hh"
#include "l_piston.hh"

/*-----------------------------------------------------------------*/

int main(int argc, char **argv)
{
    int js, jp, jd;
    std::vector<GasSlug> A;               /* several gas slugs        */
    std::vector<PistonData> Pist;          /* room for several pistons */
    std::vector<DiaphragmData> Diaph;      /* diaphragms            */

    double tstop;
    int i, max_sol;
    char oname[40], base_file_name[32], name_tag[10];
    char pname[40], iname[40];
    FILE *infile,            /* flow solution file           */
        *outfile;           /* plot file                    */

#   define  BLOCKS  10
    double *xarray[BLOCKS], *tarray[BLOCKS], *varray[BLOCKS];
    double *aarray[BLOCKS], *marray[BLOCKS];
    int nt, nttot, ip;
    int command_line_error;
    int echo_input;

    /*
     * INITIALIZE
     */
    printf("\n---------------------------");
    printf("\nLagrangian 1D Postprocessor");
    printf("\nPiston Position-Time Plots ");
    printf("\n---------------------------");
    printf("\n\n");

    strcpy(base_file_name, "default");
    max_sol = 5000;
    tstop = 0.0;
    ip = 0;
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
            /* Set the target solution time. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%lf", &tstop);
            i++;
            printf("Setting tstop = %e\n", tstop);
        } else if (strcmp(argv[i], "-p") == 0) {
            /* Set the piston index. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%d", &ip);
            i++;
            printf("Setting ip = %d\n", ip);
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
        printf("Pick up a solution file <base_file_name>.Ls, extract\n");
        printf("the data for a particular piston and save that data\n");
        printf("in a format suitable for GNU-Plot.\n");
        printf("This file will be named <base_file_name>_p<ip>.dat.\n");
        printf("\n");

        printf("Command-line options: (defaults are shown in parentheses)\n");
        printf("-f <base_file_name>       (default)\n");
        printf("-echo                     (no echo)\n");
        printf("-tstop <time>             (0.0)\n");
        printf("-p <ip>                   (0)\n");
        printf("-help          (print this message)\n");
        exit(1);
    } // end if command_line_error
    
    // Read the input parameter file.
    strcpy(pname, base_file_name);
    strcat(pname, ".Lp");
    SimulationData SD = SimulationData(pname, echo_input);
    for (jp = 0; jp < SD.npiston; ++jp) {
        Pist.push_back(PistonData(jp, SD.dt_init, pname, echo_input));
        Pist[jp].sim_time = 0.0;
    }
    for (jd = 0; jd < SD.ndiaphragm; ++jd) {
        Diaph.push_back(DiaphragmData(jd, pname, echo_input));
        Diaph[jd].sim_time = 0.0;
    }
    ConfigParser parameterdict = ConfigParser(pname);
    for (js = 0; js < SD.nslug; ++js) {
        A.push_back(GasSlug(js, SD, pname, echo_input));
        A[js].sim_time = 0.0;
    }

    /*
     * Allocate enough memory to accumulate the data.
     */
    for (jp = 0; jp < SD.npiston; ++jp) {
        if ((xarray[jp] = (double *) calloc(max_sol,sizeof(double)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
        if ((tarray[jp] = (double *) calloc(max_sol,sizeof(double)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
        if ((varray[jp] = (double *) calloc(max_sol,sizeof(double)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
        if ((aarray[jp] = (double *) calloc(max_sol,sizeof(double)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
        if ((marray[jp] = (double *) calloc(max_sol,sizeof(double)))
            == NULL) {
            printf("Memory alloc problem.\n");
            exit(-1);
        }   /* end if */
    }   /* for (jp = 0;...  */

    /*
     * Check selected piston index.
     */
    if (ip >= SD.npiston) {
        printf("Invalid piston selected: ip = %d\n", ip);
        printf("Valid range: 0..%d\n", SD.npiston - 1);
        exit(-1);
    }   /* end if */

    /*
     * * Read all of the solutions and save the requested data.
     */
    strcpy(iname, base_file_name);
    strcat(iname, ".Ls");
    printf("infile: %s\n", iname);
    if ((infile = fopen(iname, "r")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", iname);
        exit(1);
    }   /* end if */

    nttot = 0;
    for (nt = 0; nt < max_sol; ++nt) {
        printf(".");
        fflush(stdout);
        for (jp = 0; jp < SD.npiston; ++jp) Pist[jp].read_state(infile);
        for (jd = 0; jd < SD.ndiaphragm; ++jd) Diaph[jd].read_state(infile);
        for (js = 0; js < SD.nslug; ++js)
            A[js].read_state(infile);
        ++nttot;
        for (jp = 0; jp < SD.npiston; ++jp) {
            xarray[jp][nt] = Pist[jp].x;
            varray[jp][nt] = Pist[jp].V;
            aarray[jp][nt] = Pist[jp].DVDt[0];
            marray[jp][nt] = Pist[jp].mass;
            tarray[jp][nt] = A[0].sim_time;
        }   /* end for jp */
        if (A[0].sim_time >= tstop)
            break;
    }   /* end for nt... */
    printf("\n");
    if (infile != NULL)
        fclose(infile);
    printf("Number of solutions read: %d\n", nttot);
    printf("Final time = %g\n", A[0].sim_time);

    /*
     * Write the history data for the selected piston
     * in a format suitable for GNU-Plot.
     */
    strcpy(oname, base_file_name);
    sprintf(name_tag, "_p%d", ip);
    strcat(oname, name_tag);
    strcat(oname, ".dat");
    printf("Output file: %s\n", oname);
    if ((outfile = fopen(oname, "w")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", oname);
        exit(1);
    }   /* end if */

    fprintf(outfile, "# Piston-Data\n");
    fprintf(outfile, "# %d         <== number of columns\n", 4);
    fprintf(outfile, "# time,s    <== col_1\n");
    fprintf(outfile, "# x,m       <== col_2\n");
    fprintf(outfile, "# V,m/s     <== col_3\n");
    fprintf(outfile, "# a,m/s**2  <== col_4\n");
    fprintf(outfile, "# mass,kg      <== col_5\n");
    fprintf(outfile, "# %d         <== number of plotting blocks\n", 1);
    fprintf(outfile, "# %d      <== number of samples in this block\n", nttot);
    for (nt = 0; nt < nttot; ++nt) {
        fprintf(outfile, "%e %e %e %e %e\n", tarray[ip][nt],
                xarray[ip][nt], varray[ip][nt], aarray[ip][nt], marray[ip][nt]);
    }   /* end for nt */

    if (outfile != NULL)
        fclose(outfile);

    return 0;
}   /* end function main */
