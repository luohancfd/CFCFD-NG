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
    FILE *infile,                    /* beginning flow state         */
        *areafile,                  /* specification of tube area   */
        *outfile;                   /* computed solution            */
    char pname[40], aname[40], oname[40], dname[40], base_file_name[32];
    FILE *dumpf;                    /* dump files for area and K/L  */

    int ix;
    double xx;
    double x[11000]; 
    // Not that we have made a guess for the required size; 
    // we should use a std::vector but x is used in the call 
    // to distribute_points().
    double beta1, beta2;
    int i, command_line_error;
    int echo_input;

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
    L_set_case_parameters(&SD, parameterdict, echo_input);
    TubeModel tube = TubeModel(pname, echo_input);
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
    tube.write_area(aname);
    tube.write_dump_file(dname);
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
    if (outfile != NULL) fclose(outfile);

    return SUCCESS;
}   /* end of main() */

/*===================== end of l_prep.c ======================*/
