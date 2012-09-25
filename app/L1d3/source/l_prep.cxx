/** \file l_prep.cxx
 * \ingroup l1d3
 * \brief Preprocessor for l1d.cxx.
 *
 * This program is used to initialize the flow
 * solution data for the One-Dimensional Lagrangian
 * code.
 * The parameter file is read to obtain the number of gas
 * slugs, pistons and diaphragms, and their initial states.
 * A tube-area specification file is written and an initial
 * solution file is written.
 *
 * \author PA Jacobs
 *
 * \version 1.0 --  16-Dec-91 : 
 * \version 10.21   02-Apr-97 : Up to this version, changes are described in
 *                     l_prep_old.c
 * \version 11.0    05-Apr-97 : New (more general) input format.
 * \version 11.1    13-Apr-98 : command-line options
 *                     option to cluster cells
 * \version 11.2    22-Apr-98 : option to echo input parameters
 * \version 12.0    19-Jan-99 : specify tube diameters rather than areas
 * \version 24-Jul-06, C++ port.
 */

/*-----------------------------------------------------------------*/

#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l_kernel.hh"
#include "l1d.hh"
#include "l_io.hh"
#include "l_tstep.hh"
#include "l_misc.hh"
extern "C" {
#   include "../../../lib/nm/source/roberts.h"
}

/*-----------------------------------------------------------------*/

int main(int argc, char **argv)
{
    struct simulation_data SD;
    int js, jp, jd;
    vector <slug_data> slug;    /* the slug data structure */
    vector <piston_data> Pist; /* room for several pistons */
    vector <diaphragm_data> Diaph;           /* diaphragms */
    tube_data tube;             /* the tube area spec */

    FILE *infile,                    /* beginning flow state         */
        *areafile,                  /* specification of tube area   */
        *outfile;                   /* computed solution            */
    char pname[40], aname[40], oname[40], dname[40], base_file_name[32];
    FILE *dumpf;                    /* dump files for area and K/L  */

    int ix;
    double xx;
    int iseg, found_segment;
    double x[11000]; 
    // Not that we have made a guess for the required size; 
    // we should use a std::vector but x is used in the call 
    // to distribute_points().
    double beta1, beta2;
    double real_x, alpha2, alpha3;
    double K_over_L, T_Wall;
    int i, command_line_error;
    int echo_input;
    const double myPI = 4.0*atan(1.0);

    UNUSED_VARIABLE(infile);

    /*
     * ----------
     * INITIALIZE
     * ----------
     */

    printf("\n");
    printf("--------------------------\n");
    printf("Lagrangian 1D Preprocessor\n");
    printf("--------------------------\n");
    printf("\n");

    strcpy(base_file_name, "default");
    echo_input = 0; /* by default, don't echo input */

    /*
     * Decode command line arguments.
     */
    command_line_error = 0;
    if (argc < 2) {
	command_line_error = 1;
	goto usage;
    }

    i = 1;
    while (i < argc) {
        /* process the next command-line argument */

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
            printf("Will echo input parameters...\n\n");
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
        printf("Command-line options: (defaults are shown in parentheses)\n");
        printf("-f <base_file_name>       (default)\n");
        printf("-echo                     (no echo)\n");
        printf("-help          (print this message)\n");
        exit(1);
    }   /* end if command_line_error */

    /*
     * * Build the file names.
     */
    strcpy(pname, base_file_name);
    strcat(pname, ".Lp");
    strcpy(aname, base_file_name);
    strcat(aname, ".La");
    strcpy(oname, base_file_name);
    strcat(oname, ".L0");
    strcpy(dname, base_file_name);
    strcat(dname, ".dump");
    printf("\n");
    printf("parameterfile: %s\n", pname);
    printf("areafile     : %s\n", aname);
    printf("outfile      : %s\n", oname);
    printf("dump file    : %s\n", dname);

    ConfigParser parameterdict = ConfigParser(pname);
    L_set_case_parameters(&SD, &tube, parameterdict, echo_input);
    slug.resize(SD.nslug);
    Pist.resize(SD.npiston);
    Diaph.resize(SD.ndiaphragm);
    for (jp = 0; jp < SD.npiston; ++jp) {
        set_piston_parameters(&(Pist[jp]), jp, parameterdict, SD.dt_init, echo_input);
        Pist[jp].sim_time = 0.0;
    }
    for (jd = 0; jd < SD.ndiaphragm; ++jd) {
        set_diaphragm_parameters(&(Diaph[jd]), jd, parameterdict, echo_input);
        Diaph[jd].sim_time = 0.0;
    }
    for (js = 0; js < SD.nslug; ++js) {
        L_set_slug_parameters(&(slug[js]), js, &SD, parameterdict, echo_input);
        slug[js].sim_time = 0.0;
    }
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();

    printf("Set tube area specification\n");
    tube.diam.resize(tube.n);
    tube.area.resize(tube.n);
    tube.T_Wall.resize(tube.n);
    tube.K_over_L.resize(tube.n);
    tube.x1 = tube.xb[0];
    tube.dx = (tube.xb[tube.nseg] - tube.xb[0]) / tube.n;
    for (ix = 0; ix < tube.n; ++ix) {
        real_x = tube.x1 + tube.dx * ix;

        /*
         * Locate the appropriate tube segment and interpolate area.
         * Start the search from the left and stop when the xb[iseg]
         * exceeds the given x-location.
         */
        found_segment = 0;
        for (iseg = 0; iseg <= tube.nseg; ++iseg) {
            if (real_x < tube.xb[iseg]) {
                found_segment = 1;
                /* on leaving this loop, iseg indicates the segment */
                break;
            }   /* end if */
        }   /* end for iseg... */

        if (iseg == 0) {
            /* We are upstream of the tube. */
            tube.diam[ix] = tube.Diamb[iseg];
        } else if (found_segment == 0) {
            /* We are downstream of the tube. */
            tube.diam[ix] = tube.Diamb[tube.nseg];
        } else {
            /* We are between xb[iseg-1] and xb[iseg]. */
            if (tube.linear[iseg - 1] == 1) {
                /* Linear interpolation */
                xx = (real_x - tube.xb[iseg - 1]) /
                    (tube.xb[iseg] - tube.xb[iseg - 1]);
                tube.diam[ix] = tube.Diamb[iseg - 1] * (1.0 - xx) +
                    tube.Diamb[iseg] * xx;
            } else {
                /* Cubic interpolation to give dArea/dx == 0 at ends */
                xx = (real_x - tube.xb[iseg - 1]) /
                    (tube.xb[iseg] - tube.xb[iseg - 1]);
                alpha2 = 3.0 * (tube.Diamb[iseg] - tube.Diamb[iseg - 1]);
                alpha3 = -2.0 / 3.0 * alpha2;
                tube.diam[ix] = tube.Diamb[iseg - 1] +
                    (alpha2 + alpha3 * xx) * xx * xx;
            }   /* end if */
        }   /* end if */
        tube.area[ix] = myPI * 0.25 * tube.diam[ix] * tube.diam[ix];

        /*
         * Pipe-fitting loss coefficients:
         * Most of the tube is "smooth", so assume zero, then
         * search the loss patches to see if we are within one.
         */
        K_over_L = 0.0;
        for (iseg = 0; iseg < tube.nKL; ++iseg) {
            if (real_x >= tube.xbeginK[iseg] && real_x <= tube.xendK[iseg]) {
                K_over_L = tube.K[iseg] /
                    (tube.xendK[iseg] - tube.xbeginK[iseg]);
            }   /* end if */
        }   /* end for iseg */
        tube.K_over_L[ix] = K_over_L;

        /*
         * Local variations of wall temperature:
         * Assume a nominal temperature then
         * search the loss patches to see if we are within one.
         */
        T_Wall = tube.Tnominal;
        for (iseg = 0; iseg < tube.nT; ++iseg) {
            if (real_x >= tube.xbeginT[iseg] && real_x <= tube.xendT[iseg]) {
                T_Wall = tube.Tlocal[iseg];
            }   /* end if */
        }   /* end for iseg */
        tube.T_Wall[ix] = T_Wall;

    }   /* end for ix... */

    printf("Set up gas slugs.\n");

    for (js = 0; js < SD.nslug; ++js) {
        /*
         * Distribute the gas cells along the gas slug.
         * Modified to suit sm_3d stretching functions.
         */
        if ( slug[js].cluster_to_end_1 == 1 ) {
            beta1 = slug[js].cluster_strength;
        } else {
            beta1 = 0.0;
        }
        if ( slug[js].cluster_to_end_2 == 1 ) {
            beta2 = slug[js].cluster_strength;
        } else {
            beta2 = 0.0;
        }
	// FIX-ME: one day, it will be good to use a vector for x
	// such that its size can be adjested as needed.
        distribute_points(slug[js].xbegin, slug[js].xend,
                          slug[js].nnx, x, beta1, beta2 ); 
        i = 0;
        for (ix = slug[js].ixmin - 1; ix <= slug[js].ixmax; ++ix) {
            slug[js].Cell[ix].x = x[i];
            ++i;
        }   /* end for */
        /* Compute the areas at the cell interfaces. */
        L_compute_areas(&(slug[js]), &tube);
        /* Fill each slug with uniform initial conditions. */
        L_fill_data(&(slug[js]));
        L_encode_conserved(&(slug[js]));
        /* 
         * The following line should fill in all of the 
         * extra variables.
         */
        if ( L_decode_conserved(&(slug[js])) != 0 ) {
	    printf( "Failure decoding conserved quantities for slug[%d].\n", js );
	    exit( -1 );
        }
    }   /* end for js... */

    printf("Write area file.\n");

    if ((areafile = fopen(aname, "w")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", aname);
        exit(1);
    }   /* end if */
    L_write_area(&tube, areafile);
    if (areafile != NULL)
        fclose(areafile);

    printf("Write starting solution file.\n");

    if ((outfile = fopen(oname, "w")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", oname);
        exit(1);
    }   /* end if */
    for (jp = 0; jp < SD.npiston; ++jp)
        write_piston_solution(&(Pist[jp]), outfile);
    for (jd = 0; jd < SD.ndiaphragm; ++jd)
        write_diaphragm_solution(&(Diaph[jd]), outfile);
    for (js = 0; js < SD.nslug; ++js)
        L_write_solution(&(slug[js]), outfile);
    if (outfile != NULL)
        fclose(outfile);

    printf("Writing dump file for ");
    printf("diameter, area, L_over_L, T_wall...\n");

    if ((dumpf = fopen(dname, "w")) == NULL) {
        printf("Could not open %s.\n", dname);
        exit(0);
    }   /* end if */
    fprintf(dumpf, "# x  diam  area  K_over_L  T_Wall\n");
    for (ix = 0; ix < tube.n; ++ix) {
        xx = tube.dx * ix + tube.x1;
        fprintf(dumpf, "%e %e %e %e %e\n", xx, tube.diam[ix],
                tube.area[ix], tube.K_over_L[ix], tube.T_Wall[ix]);
    }   /* end for ix... */
    fclose(dumpf);

    return 0;
}   /* end of main() */

/*===================== end of l_prep.c ======================*/
