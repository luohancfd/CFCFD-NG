/** \file exch.cxx
 * \ingroup eilmer3
 * \brief Boundary data exchange for the synchronous multiblock solver.
 *
 * \author PA Jacobs
 *
 * \version 03-Nov-96: Updated multiple-species code.
 * \version 06-Feb-01: Updated exchange functions for boundary-cell geometry.
 *            Actually delegated the copying to the helper function 
 *            copy_cell_to_cell().
 * \version 02-Mar-08 Elmer3 port
 *
 */

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "exch2d.hh"
#include "exch3d.hh"


/*------------------------------------------------------------------*/

int exchange_shared_boundary_data(int jb, int type_of_copy, size_t gtl)
{
    global_data *G = get_global_data_ptr();
    if ( G->dimensions == 2 ) {
	copy_boundary_data_2D(jb, type_of_copy, gtl);
    } else {
	copy_boundary_data_3D(jb, type_of_copy, gtl);
    }
    return SUCCESS;
} // end copy_bc_data()



/** \brief Copy the boundary data from adjacent blocks to the block
 *         presently being updated, taking account of a diaphragm on
 *         the east/west boundary if it is present.
 *
 * \param jb   : index of the present block
 * \param type_of_copy : see mb_cns.h for the allowed values.
 * \param diaphragm_block : identity of the diaphragm which has a diaphragm 
 *                          at its east boundary
 * \param diaphragm_rupture_time :
 * \param diaphragm_rupture_diameter :
 *
 */
int copy_boundary_data_2D(int jb, int type_of_copy, size_t gtl)
{
    global_data *G = get_global_data_ptr();
    Block *bdp = get_block_data_ptr(jb);
    Block *other_bdp;
    int other_block, other_bndry;
    other_block = bdp->bcp[NORTH]->neighbour_block;
    other_bdp = get_block_data_ptr(other_block);
    if (other_block >= 0) {
        other_bndry = bdp->bcp[NORTH]->neighbour_face;
        copy_to_north_boundary_2D(bdp, other_bdp, other_bndry, type_of_copy, gtl);   
    }
    /* note: the variable "diaphragm_block" defaults to -1 and so will never refer
     *       to an actual block, unless the case id for a diaphragm rupture simulation
     *       is set and the variable is assigned to a real block in mb_special_init.c
     */
 
    /* if you are the block at the upstream side of the diaphragm
       and you are looking at the diaphragm edge, go to the special
       diaphragm function */
    other_block = bdp->bcp[EAST]->neighbour_block;
    other_bdp = get_block_data_ptr(other_block);
    if (other_block >= 0) {
        other_bndry = bdp->bcp[EAST]->neighbour_face;
        if (jb == G->diaphragm_block) {
            copy_to_east_boundary_diaphragm_2D(bdp, other_bdp, other_bndry, 
					       type_of_copy, G->diaphragm_rupture_time, 
					       G->diaphragm_rupture_diameter, G->sim_time,
					       gtl);
        } else {
            copy_to_east_boundary_2D(bdp, other_bdp, other_bndry, type_of_copy, gtl);
        }
    }

    other_block = bdp->bcp[SOUTH]->neighbour_block;
    other_bdp = get_block_data_ptr(other_block);
    if (other_block >= 0) {
        other_bndry = bdp->bcp[SOUTH]->neighbour_face;
        copy_to_south_boundary_2D(bdp, other_bdp, other_bndry, type_of_copy, gtl);
    }

    /* if you are the block at the downstream side of the diaphragm
       and you are looking at the diaphragm edge, go to the special
       diaphragm function */
    other_block = bdp->bcp[WEST]->neighbour_block;
    other_bdp = get_block_data_ptr(other_block);
    if (other_block >= 0) {
        other_bndry = bdp->bcp[WEST]->neighbour_face;
        if (other_block == G->diaphragm_block) {
            copy_to_west_boundary_diaphragm_2D(bdp, other_bdp, other_bndry,
					       type_of_copy, G->diaphragm_rupture_time,
					       G->diaphragm_rupture_diameter, G->sim_time,
					       gtl);
        } else {
            copy_to_west_boundary_2D(bdp, other_bdp, other_bndry, type_of_copy, gtl);
        }
    }

    return SUCCESS;
}   /* end copy_bc_data_2D() */
        

/*-----------------------------------------------------------------*/

double calculate_diaphragm_radius( double sim_time,
				   double diaphragm_rupture_time,
				   double diaphragm_rupture_diameter,
				   double minimum_radius)
{
    double diaphragm_time_fraction, diaphragm_area, diaphragm_radius;

    if ( diaphragm_rupture_time == 0.0 ) {
	diaphragm_time_fraction = 0.0;
    } else {
	diaphragm_time_fraction = sim_time/diaphragm_rupture_time;
    } 

    /*
     * Linearly opening iris based rupture 
     * (experimental opening vs time curve of Rothkopf and Low (1974))
     */  
    if (diaphragm_time_fraction > 1.0) diaphragm_time_fraction = 1.0;
    diaphragm_area = diaphragm_time_fraction * 
	(3.14159265 * diaphragm_rupture_diameter * diaphragm_rupture_diameter/4.0);
    diaphragm_radius = sqrt(diaphragm_area/3.14159265);

    if ( diaphragm_radius < minimum_radius )
	diaphragm_radius = minimum_radius + 1.0e-4;

#   if 0
    /* this is for the new trial diaphragm with the multiple openings */
    /* BUT it won't work in the current function arrangement!! */
    if (diaphragm_time_fraction > 1.0) diaphragm_time_fraction = 1.0;
    if (cell_radius < 0.0 + diaphragm_time_fraction*5.18e-3 ||
	cell_radius > 12.96e-3 - diaphragm_time_fraction*4.32e-3 && 
	cell_radius < 12.96e-3 + diaphragm_time_fraction*4.32e-3 ||
	cell_radius > 24.19e-3 - diaphragm_time_fraction*3.46e-3 && 
	cell_radius < 24.19e-3 + diaphragm_time_fraction*3.46e-3 ) {
	    diaphragm_cell = 0;
	} else {
	    diaphragm_cell = 1;
	}
 #   endif

    return diaphragm_radius;
}   /* end calculate_diaphragm_radius() */


/** \brief Copy data from the B_bndry (WEST only!)
 *         of (source) block B to the EAST boundary of (target) block A,
 *         taking into account the amount that the diaphragm has opened..
 *
 *  \param  A       : pointer to the target block
 *  \param  B       : pointer to the source block
 *  \param B_bndry  : value specifying the source block boundary
 */
int copy_to_east_boundary_diaphragm_2D(Block *A,
				       Block *B,
				       int B_bndry,
				       int type_of_copy,
				       double diaphragm_rupture_time,
				       double diaphragm_rupture_diameter,
				       double sim_time, size_t gtl)
{
    int i_A, j_A, i_B, j_B;
    int jfirst, jlast;
    FV_Cell *src, *dest;
    double cell_radius, diaphragm_radius;

    if (B_bndry == WEST) {
        /*
         * Copy data from cells adjacent to the WEST boundary of
         * block B to the EAST boundary ghost cells of block A.
         */
	diaphragm_radius = 
	    calculate_diaphragm_radius(sim_time, diaphragm_rupture_time,
				       diaphragm_rupture_diameter,
				       B->get_ifj(B->imin,B->jmin)->Ybar);

        i_A = A->imax + 1;
        i_B = B->imin;
        jfirst = A->jmin;
        jlast = A->jmax;

        for (j_A = jfirst; j_A <= jlast; ++j_A) {
	    j_B = B->jmin + (j_A - A->jmin);
	    cell_radius = B->get_ifj(A->imin,j_A)->Ybar;

	    if (cell_radius < diaphragm_radius) {
		/* 
		 * This cell y index does not have diaphragm material next to it.
		 * It gets data from the neighbouring block like a normal boundary.
		 */
		/* Fill the first line of ghost cells. */
		dest = A->get_cell(i_A,j_A);
		src = B->get_cell(i_B,j_B);
		dest->copy_values_from(*src, type_of_copy, gtl);
		/* Fill the second line of ghost cells. */
		dest = A->get_cell(i_A+1,j_A);
		src = B->get_cell(i_B+1,j_B);
		dest->copy_values_from(*src, type_of_copy, gtl);
	    } else {
		/*
		 * This cell y index has diaphragm material next to it.
		 * It gets data from itself as a reflected boundary.
		 */
		/* Fill the first line of ghost cells. */
		dest = A->get_cell(i_A,j_A);
		src = A->get_cell(A->imax,j_A);
		dest->copy_values_from(*src, type_of_copy, gtl);
		/* Fill the second line of ghost cells. */
		dest = A->get_cell(i_A+1,j_A);
		src = A->get_cell(A->imax-1,j_A);
		dest->copy_values_from(*src, type_of_copy, gtl);
	    }   /* end if */
        }   /* end for */

    } else {
        printf("\ncopy_to_east_boundary_diaphragm: invalid boundary\n");
        exit(VALUE_ERROR);
    }   /* end if */

    return SUCCESS;
}   /* end copy_to_east_boundary_diaphragm() */


/** \brief Copy data from the B_bndry (EAST only)
 *         of (source) block B to the WEST boundary of (target) block A.
 *
 *  \param  A       : pointer to the target block
 *  \param  B       : pointer to the source block
 *  \param B_bndry  : value specifying the source block boundary
 */
int copy_to_west_boundary_diaphragm_2D(Block *A,
				       Block *B,
				       int B_bndry,
				       int type_of_copy,
				       double diaphragm_rupture_time,
				       double diaphragm_rupture_diameter,
				       double sim_time, size_t gtl)
{
    int i_A, j_A, i_B, j_B;
    int jfirst, jlast;
    FV_Cell *src, *dest;
    double cell_radius, diaphragm_radius;

    if (B_bndry == EAST) {
        /*
         * Copy data from cells adjacent to the EAST boundary of
         * block B to the WEST boundary ghost cells of block A.
         */
	diaphragm_radius = 
	    calculate_diaphragm_radius(sim_time, diaphragm_rupture_time,
				       diaphragm_rupture_diameter,
				       A->get_ifj(A->imin,A->jmin)->Ybar);

        i_A = A->imin - 1;
        i_B = B->imax;
        jfirst = A->jmin;
        jlast = A->jmax;

        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            j_B = B->jmin + (j_A - A->jmin);
            cell_radius = A->get_ifj(A->imin,j_A)->Ybar;

	    if (cell_radius < diaphragm_radius) {
		/*
		 * This cell y index does not have diaphragm material next to it.
		 * It gets data from the neighbouring block like a normal boundary.
		 */
		/* Fill the first line of ghost cells. */
		dest = A->get_cell(i_A,j_A);
		src = B->get_cell(i_B,j_B);
		dest->copy_values_from(*src, type_of_copy, gtl);
		/* Fill the second line of ghost cells. */
		dest = A->get_cell(i_A-1,j_A);
		src = B->get_cell(i_B-1,j_B);
		dest->copy_values_from(*src, type_of_copy, gtl);
	    } else {
		/* 
		 * This cell y index has diaphragm material next to it and so it gets
		 * data from itself as a reflected boundary
		 */
		/* Fill the first line of ghost cells. */
		dest = A->get_cell(i_A,j_A);
		src = A->get_cell(A->imin,j_A);
		dest->copy_values_from(*src, type_of_copy, gtl);
		/* Fill the second line of ghost cells. */
		dest = A->get_cell(i_A-1,j_A);
		src = A->get_cell(A->imin+1,j_A);
		dest->copy_values_from(*src, type_of_copy, gtl);
	    }   /* end if */
        }   /* end for */

    } else {
        printf("\ncopy_to_west_boundary_diaphragm: invalid boundary\n");
        exit(VALUE_ERROR);
    }   /* end if */

    return SUCCESS;
}   /* end copy_to_west_boundary_diaphragm() */


/** \brief Copy data from the B_bndry (NORTH, EAST, SOUTH, WEST)
 *         of (source) block B to the NORTH boundary of (target) block A.
 *
 *  \param  A       : pointer to the target block
 *  \param  B       : pointer to the source block
 *  \param B_bndry  : value specifying the source block boundary
 *
 */
int copy_to_north_boundary_2D(Block *A, Block *B, int B_bndry,
			      int type_of_copy, size_t gtl)
{
    int i_A, j_A, i_B, j_B;
    int ifirst, ilast;
    FV_Cell *src, *dest;
    if (B_bndry == SOUTH) {
        /*
         * Copy data from cells adjacent to the SOUTH boundary of
         * block B to the NORTH boundary ghost cells of block A.
         */
        j_A = A->jmax + 1;
        j_B = B->jmin;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            i_B = B->imin + (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A+1);
            src = B->get_cell(i_B,j_B+1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == EAST) {
        /*
         * Copy data from cells adjacent to the EAST boundary of
         * block B to the NORTH boundary ghost cells of block A.
         */
        j_A = A->jmax + 1;
        i_B = B->imax;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            j_B = B->jmin + (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A+1);
            src = B->get_cell(i_B-1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == NORTH) {
        /*
         * Copy data from cells adjacent to the NORTH boundary of
         * block B to the NORTH boundary ghost cells of block A.
         */
        j_A = A->jmax + 1;
        j_B = B->jmax;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            i_B = B->imax - (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A+1);
            src = B->get_cell(i_B,j_B-1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == WEST) {
        /*
         * Copy data from cells adjacent to the WEST boundary of
         * block B to the NORTH boundary ghost cells of block A.
         */
        j_A = A->jmax + 1;
        i_B = B->imin;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            j_B = B->jmax - (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A+1);
            src = B->get_cell(i_B+1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else {
        printf("\ncopy_to_north_boundary: invalid boundary\n");
        exit(VALUE_ERROR);
    }   /* end if */

    return SUCCESS;
}   /* end copy_to_north_boundary_2D() */


/** \brief Copy data from the B_bndry (NORTH, EAST, SOUTH, WEST)
 *         of (source) block B to the EAST boundary of (target) block A.
 *
 *  \param  A       : pointer to the target block
 *  \param  B       : pointer to the source block
 *  \param B_bndry  : value specifying the source block boundary
 */
int copy_to_east_boundary_2D(Block *A, Block *B, int B_bndry,
			     int type_of_copy, size_t gtl)
{
    int i_A, j_A, i_B, j_B;
    int jfirst, jlast;
    FV_Cell *src, *dest;

    if (B_bndry == SOUTH) {
        /*
         * Copy data from cells adjacent to the SOUTH boundary of
         * block B to the EAST boundary ghost cells of block A.
         */
        i_A = A->imax + 1;
        j_B = B->jmin;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            i_B = B->imax - (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A+1,j_A);
            src = B->get_cell(i_B,j_B+1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == EAST) {
        /*
         * Copy data from cells adjacent to the EAST boundary of
         * block B to the EAST boundary ghost cells of block A.
         */
        i_A = A->imax + 1;
        i_B = B->imax;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            j_B = B->jmax - (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A+1,j_A);
            src = B->get_cell(i_B-1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == NORTH) {
        /*
         * Copy data from cells adjacent to the NORTH boundary of
         * block B to the EAST boundary ghost cells of block A.
         */
        i_A = A->imax + 1;
        j_B = B->jmax;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            i_B = B->imin + (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A+1,j_A);
            src = B->get_cell(i_B,j_B-1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == WEST) {
        /*
         * Copy data from cells adjacent to the WEST boundary of
         * block B to the EAST boundary ghost cells of block A.
         */
        i_A = A->imax + 1;
        i_B = B->imin;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            j_B = B->jmin + (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A+1,j_A);
            src = B->get_cell(i_B+1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else {
        printf("\ncopy_to_east_boundary: invalid boundary\n");
        exit(VALUE_ERROR);
    }   /* end if */

    return SUCCESS;
}   /* end copy_to_east_boundary_2D() */


/** \brief Copy data from the B_bndry (NORTH, EAST, SOUTH, WEST)
 *         of (source) block B to the SOUTH boundary of (target) block A.
 *
 *  \param  A       : pointer to the target block
 *  \param  B       : pointer to the source block
 *  \param B_bndry  : value specifying the source block boundary
 */
int copy_to_south_boundary_2D(Block *A, Block *B, int B_bndry,
			      int type_of_copy, size_t gtl)
{
    int i_A, j_A, i_B, j_B;
    int ifirst, ilast;
    FV_Cell *src, *dest;

    if (B_bndry == SOUTH) {
        /*
         * Copy data from cells adjacent to the SOUTH boundary of
         * block B to the SOUTH boundary ghost cells of block A.
         */
        j_A = A->jmin - 1;
        j_B = B->jmin;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            i_B = B->imax - (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A-1);
            src = B->get_cell(i_B,j_B+1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == EAST) {
        /*
         * Copy data from cells adjacent to the EAST boundary of
         * block B to the SOUTH boundary ghost cells of block A.
         */
        j_A = A->jmin - 1;
        i_B = B->imax;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            j_B = B->jmax - (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A-1);
            src = B->get_cell(i_B-1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == NORTH) {
        /*
         * Copy data from cells adjacent to the NORTH boundary of
         * block B to the SOUTH boundary ghost cells of block A.
         */
        j_A = A->jmin - 1;
        j_B = B->jmax;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            i_B = B->imin + (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A-1);
            src = B->get_cell(i_B,j_B-1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == WEST) {
        /*
         * Copy data from cells adjacent to the WEST boundary of
         * block B to the SOUTH boundary ghost cells of block A.
         */
        j_A = A->jmin - 1;
        i_B = B->imin;

        ifirst = A->imin;
        ilast = A->imax;
        for (i_A = ifirst; i_A <= ilast; ++i_A) {
            j_B = B->jmin + (i_A - A->imin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A,j_A - 1);
            src = B->get_cell(i_B+1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else {
        printf("\ncopy_to_south_boundary: invalid boundary\n");
        exit(VALUE_ERROR);
    }   /* end if */

    return SUCCESS;
}   /* end copy_to_south_boundary_2D() */


/** \brief Copy data from the B_bndry (NORTH, EAST, SOUTH, WEST)
 *         of (source) block B to the WEST boundary of (target) block A.
 *
 *  \param  A       : pointer to the target block
 *  \param  B       : pointer to the source block
 *  \param B_bndry  : value specifying the source block boundary
 */
int copy_to_west_boundary_2D(Block *A, Block *B, int B_bndry,
			     int type_of_copy, size_t gtl)
{
    int i_A, j_A, i_B, j_B;
    int jfirst, jlast;
    FV_Cell *src, *dest;

    if (B_bndry == SOUTH) {
        /*
         * Copy data from cells adjacent to the SOUTH boundary of
         * block B to the WEST boundary ghost cells of block A.
         */
        i_A = A->imin - 1;
        j_B = B->jmin;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            i_B = B->imin + (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A-1,j_A);
            src = B->get_cell(i_B,j_B+1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == EAST) {
        /*
         * Copy data from cells adjacent to the EAST boundary of
         * block B to the WEST boundary ghost cells of block A.
         */
        i_A = A->imin - 1;
        i_B = B->imax;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            j_B = B->jmin + (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A-1,j_A);
            src = B->get_cell(i_B-1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == NORTH) {
        /*
         * Copy data from cells adjacent to the NORTH boundary of
         * block B to the WEST boundary ghost cells of block A.
         */
        i_A = A->imin - 1;
        j_B = B->jmax;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            i_B = B->imax - (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A-1,j_A);
            src = B->get_cell(i_B,j_B-1);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else if (B_bndry == WEST) {
        /*
         * Copy data from cells adjacent to the WEST boundary of
         * block B to the WEST boundary ghost cells of block A.
         */
        i_A = A->imin - 1;
        i_B = B->imin;

        jfirst = A->jmin;
        jlast = A->jmax;
        for (j_A = jfirst; j_A <= jlast; ++j_A) {
            j_B = B->jmax - (j_A - A->jmin);

            /* Fill the first line of ghost cells. */
            dest = A->get_cell(i_A,j_A);
            src = B->get_cell(i_B,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);

            /* Fill the second line of ghost cells. */
            dest = A->get_cell(i_A-1,j_A);
            src = B->get_cell(i_B+1,j_B);
            dest->copy_values_from(*src, type_of_copy, gtl);
        }   /* end for */

    } else {
        printf("\ncopy_to_west_boundary: invalid boundary\n");
        exit(VALUE_ERROR);
    }   /* end if */

    return SUCCESS;
}   /* end copy_to_west_boundary_2D() */


/*--------------------------------------------------------------------*/
/*                  Copy-via-buffer exchange functions.               */ 
/*--------------------------------------------------------------------*/

/** \brief Copy data into the send buffer from the appropriate
 *         boundary of the current block.
 *
 * See workbook page 33, 18-Jun-02 for details of order.
 *
 * \param bd           : pointer to the target block
 * \param bndry        : value specifying the source block boundary
 * \param type_of_copy : sometimes we want to copy cell geometry,
 *                       other times the flow data
 */
int copy_into_send_buffer_2D(Block *bd, int bndry, int type_of_copy,
			     double *send_buffer, size_t gtl)
{
    int i, j, ib, ii, nx, ny, nv;
    FV_Cell *cell;
    double *bufp;
    nv = number_of_values_in_cell_copy(type_of_copy);

    if (bndry == NORTH) {
        j = bd->jmax;
        nx = bd->nni;
        for (ii = 0; ii < nx; ++ii) {
            /* copy counter-clockwise */
            i = bd->imax - ii;
            /* copy the outer-most line of active cells. */
            ib = ii;
            cell = bd->get_cell(i,j);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
            /* copy the second line of active cells. */
            ib = ii + nx;
            cell = bd->get_cell(i,j-1);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
        }
    } else if (bndry == EAST) {
        i = bd->imax;
        ny = bd->nnj;
        for (ii = 0; ii < ny; ++ii) {
            /* copy counter-clockwise */
            j = bd->jmin + ii;
            /* copy the outer-most line of active cells. */
            ib = ii;
            cell = bd->get_cell(i,j);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
            /* copy the second line of active cells. */
            ib = ii + ny;
            cell = bd->get_cell(i-1,j);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	}
    } else if (bndry == SOUTH) {
        j = bd->jmin;
        nx = bd->nni;
        for (ii = 0; ii < nx; ++ii) {
            /* copy counter-clockwise */
            i = bd->imin + ii;
            /* copy the outer-most line of active cells. */
            ib = ii;
            cell = bd->get_cell(i,j);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
            /* copy the second line of active cells. */
            i = bd->imin + ii;
            ib = ii + nx;
            cell = bd->get_cell(i,j+1);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
        }
    } else if (bndry == WEST) {
        i = bd->imin;
        ny = bd->nnj;
        for (ii = 0; ii < ny; ++ii) {
            /* copy counter-clockwise */
            j = bd->jmax - ii;
            /* copy the outer-most line of active cells. */
            ib = ii;
            cell = bd->get_cell(i,j);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
            /* copy the second line of active cells. */
            ib = ii + ny;
            cell = bd->get_cell(i+1,j);
            bufp = &(send_buffer[ib * nv]);
            cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	}
    } else {
        printf("\ncopy_into_send_buffer_2D(): invalid boundary\n");
        return VALUE_ERROR;
    } // end if

    return SUCCESS;
} // end copy_into_send_buffer_2D()


/** \brief Copy data from the receive buffer into the appropriate
 *         boundary of the current block.
 *
 * See workbook page 33, 18-Jun-02 for details of order.
 */
int copy_from_receive_buffer_2D(Block *bd, int bndry, int type_of_copy,
				double *receive_buffer, size_t gtl)
{
    int i, j, ib, ii, nx, ny, nv;
    FV_Cell *cell;
    double *bufp;
    nv = number_of_values_in_cell_copy(type_of_copy);
    if (bndry == NORTH) {
        j = bd->jmax;
        nx = bd->nni;
        for (ii = 0; ii < nx; ++ii) {
            /* Fill clockwise */
            i = bd->imin + ii;
            /* Fill the first line of ghost cells. */
            ib = ii;
            cell = bd->get_cell(i,j+1);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
            /* Fill the second line of ghost cells. */
            ib = ii + nx;
            cell = bd->get_cell(i,j+2);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
        }
    } else if (bndry == EAST) {
        i = bd->imax;
        ny = bd->nnj;
        for (ii = 0; ii < ny; ++ii) {
            /* Fill clockwise */
            j = bd->jmax - ii;
            /* Fill the first line of ghost cells. */
            ib = ii;
            cell = bd->get_cell(i+1,j);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
            /* Fill the second line of ghost cells. */
            ib = ii + ny;
            cell = bd->get_cell(i+2,j);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
	}
    } else if (bndry == SOUTH) {
        j = bd->jmin;
        nx = bd->nni;
        for (ii = 0; ii < nx; ++ii) {
            /* Fill clockwise */
            i = bd->imax - ii;
            /* Fill the first line of ghost cells. */
            ib = ii;
            cell = bd->get_cell(i,j-1);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
            /* Fill the second line of ghost cells. */
            ib = ii + nx;
            cell = bd->get_cell(i,j-2);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
        }
    } else if (bndry == WEST) {
        i = bd->imin;
        ny = bd->nnj;
        for (ii = 0; ii < ny; ++ii) {
            /* Fill clockwise */
            j = bd->jmin + ii;
            /* Fill the first line of ghost cells. */
            ib = ii;
            cell = bd->get_cell(i-1,j);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
            /* Fill the second line of ghost cells. */
            ib = ii + ny;
            cell = bd->get_cell(i-2,j);
            bufp = &(receive_buffer[ib * nv]);
            cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
	}
    } else {
        printf("\ncopy_from_receive_buffer: invalid boundary\n");
        return VALUE_ERROR;
    } // end if

    return SUCCESS;
} // end copy_from_receive_buffer_2D()

