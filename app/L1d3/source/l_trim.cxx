/** \file l_trim.cxx
 * \ingroup l1d2
 * \brief Part of the postprocessing set for l1d.c.
 *
 * This program is used to read the data for a single gas slug
 * and then trim that slug according to user-specified limits.
 *
 * \author PA Jacobs
 *
 * \version 1.0 --  21-Aug-99
 * \version 24-Jul-06, C++ port
 */

/*-----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/util/source/config_parser.hh"
#include "l_kernel.hh"
#include "l1d.hh"
#include "l_io.hh"
#include "l_misc.hh"

int L_read_solution_no_check(struct slug_data *A, FILE * infile);
int L_trim_slug(struct slug_data *A, double xmin, double xmax);

/*-----------------------------------------------------------------*/

int main(int argc, char **argv)
{
    struct simulation_data SimData;
    struct slug_data sd;
    tube_data tube;
    char file_name[40], base_file_name[32], pname[40];
    FILE *fp;
    int i, command_line_error, result_flag;
    double xmin, xmax;

    /*
     * INITIALIZE
     */
    printf("\n---------------------------------------------------");
    printf("\nLagrangian 1D Solution Postprocessor: Trim Gas Slug");
    printf("\n---------------------------------------------------");
    printf("\n\n");

    /*
     * Defaults.
     */
    strcpy(base_file_name, "default");
    strcpy(file_name, "slug");
    xmin = -9.9e10;
    xmax = 9.9e10;

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
            printf("Setting base file_name = %s\n", base_file_name);
        } else if (strcmp(argv[i], "-fs") == 0) {
            /* Set the slug file name. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            strcpy(file_name, argv[i]);
            i++;
            printf("Setting slug file_name = %s\n", file_name);
        } else if (strcmp(argv[i], "-xmin") == 0) {
            /* Set the minimum x-location. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%lf", &xmin);
            i++;
            printf("Setting xmin = %e\n", xmin);
        } else if (strcmp(argv[i], "-xmax") == 0) {
            /* Set the maximum x-location. */
            i++;
            if (i >= argc) {
                command_line_error = 1;
                goto usage;
            }
            sscanf(argv[i], "%lf", &xmax);
            i++;
            printf("Setting xmax = %e\n", xmax);
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
        printf("Pick up a gas slug and trim it such that all cells\n");
        printf("lie completely in the range xmin <= x <= xmax.\n");
        printf("\n");

        printf("Command-line options: (defaults are shown in parentheses)\n");
        printf("-f <base_file_name>       (default)\n");
        printf("-fs <slug_file_name>         (slug)\n");
        printf("-xmin <x>                 (-9.9e10)\n");
        printf("-xmax <x>                 (+9.9e10)\n");
        printf("-help          (print this message)\n");
        exit(1);
    }
    /* end if command_line_error */

    /*
     ** Read the input parameter file.
     */
    strcpy(pname, base_file_name);
    strcat(pname, ".Lp");
    printf("parameterfile: %s\n", pname);
    ConfigParser parameterdict = ConfigParser(pname);
    L_set_case_parameters(&SimData, &tube, parameterdict, 0);
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    printf("After reading case parameters, nsp=%d\n", nsp);

    printf("Read the current gas slug file.\n");
    if ((fp = fopen(file_name, "r")) == NULL) {
        printf("\nCould not open %s; BAILING OUT\n", file_name);
        exit(1);
    }   /* end if */
    nsp = L_read_solution_no_check(&sd, fp);
    if (fp != NULL) fclose(fp);

    result_flag = L_trim_slug(&sd, xmin, xmax);

    if (result_flag == 0) {
        printf("Write the new slug out.\n");
        if ((fp = fopen(file_name, "w")) == NULL) {
            printf("\nCould not open %s; BAILING OUT\n", file_name);
            exit(-1);
        }   /* end if */
        L_write_solution(&sd, fp, nsp);
        if (fp != NULL) fclose(fp);
    } else {
	printf( "Nothing left to write out.\n" );
    }   /* end if */

    return 0;
}   /* end function main */

/*----------------------------------------------------*/

int L_read_solution_no_check(struct slug_data *A, FILE * infile)
{
    /* 
     * Read the flow solution (i.e. the primary variables at the 
     * cell centers) from a disc file.
     * Note that it relys entirely on information in that file; 
     * does not access the .Lp file.
     * 
     * Returns nsp as read from the top of the slug data.
     */
#   define NCHAR 320
    char line[NCHAR], token[NCHAR];
    int ix, nsp, isp, nread;
    double xx, f_sum;
    struct L_cell *c;

    printf("Reading the slug data...\n");

    if (fgets(line, NCHAR, infile) == NULL) {
        printf("Empty flow field file.\n");
        exit(0);
    }   /* end if */
    nread = sscanf(line, "%lf", &(A->sim_time));
    printf("Time = %e\n", A->sim_time);

    /*
     * From here on, we assume that the solution file is ok.
     * i.e. NO checks on the data.
     */

    if ( fgets(line, NCHAR, infile) == NULL ) {
	printf("Problem reading file.\n");
	exit(BAD_INPUT_ERROR);
    }
    nread = sscanf(line, "%d %d", &(A->nnx), &nsp);
    printf("nnx = %d, nsp = %d\n", A->nnx, nsp);
    A->nghost = 2;
    A->nxdim = A->nnx + 2 * A->nghost;
    A->ixmin = 2;
    A->ixmax = A->ixmin + A->nnx - 1;
    if (L_alloc(A) != 0) {
        printf("Could not allocate memory for slug data.\n");
        exit(1);
    }   /* end if */

    for (ix = A->ixmin - 1; ix <= A->ixmax; ++ix) {
	c = &( A->Cell[ix] );
	if ( fgets(line, NCHAR, infile) == NULL ) {
	    printf("Problem reading file.\n");
	    exit(BAD_INPUT_ERROR);
	}
        nread = sscanf(line, "%lf %lf", &(c->x), &(c->area));
	if ( nread != 2 ) {
	    printf( "Cell interface[%d] failed to read x and area:\n", ix );
	    printf( "%s\n", line );
	    exit( -1 );
	}
    }   /* end for */

    for (ix = A->ixmin; ix <= A->ixmax; ++ix) {
	c = &( A->Cell[ix] );

	if ( fgets(line, NCHAR, infile) == NULL ) {
	    printf("Problem reading file.\n");
	    exit(BAD_INPUT_ERROR);
	}
        nread = sscanf(line, "%lf %lf %lf %lf %lf", 
		       &xx, &(c->gas->rho), &(c->u), &(c->gas->e[0]), &(c->L_bar) );
	if ( nread != 5 ) {
	    printf( "Cell[%d] failed to read x, rho, u, e, L_bar from line:\n", ix );
	    printf( "%s\n", line );
	    exit( -1 );
	}

	if ( fgets(line, NCHAR, infile) == NULL ) {
	    printf("Problem reading file.\n");
	    exit(BAD_INPUT_ERROR);
	}
        nread = sscanf(line, "%lf %lf %lf %lf %lf %lf", 
		       &(c->gas->p), &(c->gas->a), &(c->gas->T[0]), 
                       &(c->shear_stress), &(c->heat_flux), &(c->entropy));
	if ( nread != 6 ) {
	    printf( "Cell[%d] failed to read p, a, T, tau0, heatflux, entropy from line:\n", ix );
	    printf( "%s\n", line );
	    exit( -1 );
	}

	/*
         * Now get the species mass fractions from one line of text
	 * and check their sum.
         */
	for ( isp = 0; isp < nsp; ++isp ) {
	    c->gas->massf[isp] = 0.0;
	}
	if ( fgets(line, NCHAR, infile) == NULL ) {
	    printf("Problem reading file.\n");
	    exit(BAD_INPUT_ERROR);
	}
	strcpy( token, strtok( line, " " ) );
	nread = sscanf( token, "%lf", &(c->gas->massf[0]) );
	if ( nread != 1 ) {
	    printf("Cell[%d] problem reading species[%d] from token:%s\n",
		   ix, 0, token);
	    exit(-1);
	}
	f_sum = c->gas->massf[0];
	for ( isp = 1; isp < nsp; ++isp ) {
	    strcpy( token, strtok( NULL, " " ) );
	    nread = sscanf( token, "%lf", &(c->gas->massf[isp]) );
	    if ( nread != 1 ) {
		printf("Cell[%d] problem reading species[%d] from token:%s\n",
		       ix, isp, token);
		exit(-1);
	    }
	    f_sum += c->gas->massf[isp];
	}
	if ( fabs(f_sum - 1.0) > 0.001 ) {
	    printf("Species don't sum correctly %g\n", f_sum );
	}
    }   /* end for */

    printf("Range of x values : %g %g\n",
           A->Cell[A->ixmin - 1].x, A->Cell[A->ixmax].x);

    return nsp;
#   undef NCHAR
}   /* end function L_read_solution_no_check() */

/*----------------------------------------------------*/

int L_trim_slug(struct slug_data *A, double xmin, double xmax)
{
    /* 
     * Purpose...
     * -------
     * Remove cells from the slug until all remaining cells
     * lie within the range xmin <= x <= xmax.
     *
     * Input...
     * -----
     * A          : pointer to the slug data structure
     * xmin, xmax : user-specified limits for the cell positions
     */

    printf("Trimming the slug data...\n");

    /*
     * Check that we have some work to do.
     */
    if (xmin >= xmax) {
        printf("Invalid x-range specified.\n");
        return -1;
    }   /* end if */
    if (xmin > A->Cell[A->ixmax].x) {
        printf("Specified range lies completely to right of slug.\n");
        return -1;
    }   /* end if */
    if (xmax < A->Cell[A->ixmin - 1].x) {
        printf("Specified range lies completely to left of slug.\n");
        return -1;
    }   /* end if */

    /*
     * * Start looking at the cell interfaces from the left
     * * Stop when we are within left end of the specified range.
     */
    while (A->Cell[A->ixmin - 1].x < xmin) {
        ++(A->ixmin);
    };

    /*
     * Do the equivalent from the right.
     */
    while (A->Cell[A->ixmax].x > xmax) {
        --(A->ixmax);
    };

    /* 
     * Recompute the number of cells.
     */
    A->nnx = A->ixmax - A->ixmin + 1;
    printf("New number of cells : %d\n", A->nnx);

    printf("New range of x values : %g %g\n",
           A->Cell[A->ixmin - 1].x, A->Cell[A->ixmax].x);

    return 0;
#   undef NCHAR
}   /* end function L_trim_slug() */

/*=============== end of l_trim.c ==================*/
