/** \file exch3d.cxx
 * \ingroup eilmer3
 * \brief Functions to copy boundary data from one 3D block to another.
 *
 * \author PJ
 * \version July 2008, ported from eilmer2.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cell.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "exch3d.hh"

/*--------------------------------------------------------------------*/
/*                  Direct-copy exchange functions.                   */ 
/*--------------------------------------------------------------------*/

/** \brief Coordinates the filling of ghost-cell data from ADJACENT blocks
 *         in a serial calculation where all blocks are in the one data space.
 */
int copy_boundary_data_3D(size_t jb, int type_of_copy, size_t gtl)
{
    Block *bdp = get_block_data_ptr(jb);
    int other_block;

    other_block = bdp->bcp[NORTH]->neighbour_block;
    if (other_block >= 0) {
	copy_into_north_boundary_3D(bdp, get_block_data_ptr(other_block), type_of_copy, gtl);   
    }
    other_block = bdp->bcp[EAST]->neighbour_block;
    if (other_block >= 0) {
	copy_into_east_boundary_3D(bdp, get_block_data_ptr(other_block), type_of_copy, gtl);   
    }
    other_block = bdp->bcp[SOUTH]->neighbour_block;
    if (other_block >= 0) {
	copy_into_south_boundary_3D(bdp, get_block_data_ptr(other_block), type_of_copy, gtl);   
    }
    other_block = bdp->bcp[WEST]->neighbour_block;
    if (other_block >= 0) {
	copy_into_west_boundary_3D(bdp, get_block_data_ptr(other_block), type_of_copy, gtl);   
    }
    other_block = bdp->bcp[TOP]->neighbour_block;
    if (other_block >= 0) {
	copy_into_top_boundary_3D(bdp, get_block_data_ptr(other_block), type_of_copy, gtl);   
    }
    other_block = bdp->bcp[BOTTOM]->neighbour_block;
    if (other_block >= 0) {
	copy_into_bottom_boundary_3D(bdp, get_block_data_ptr(other_block), type_of_copy, gtl);   
    }
    return SUCCESS;
}


/** \brief Do a straight copy of boundary data into the east-boundary
 *         ghost cells.
 */
int copy_into_east_boundary_3D(Block *bp, Block *bp_src, int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t i_dest, i_src, j_dest, j_src, k_dest, k_src, j, k;
    FV_Cell *src, *dest;
    int neighbour_faceId, orientation;

    i_dest = bp->imax;  /* index of the east-most plane of cells */
    neighbour_faceId = bp->bcp[EAST]->neighbour_face;
    orientation = bp->bcp[EAST]->neighbour_orientation;

    if ( neighbour_faceId == WEST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    j_src = k; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 2 ) {
		    j_src = bp_src->nnj - j - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = bp_src->nnj - k - 1; k_src = j;
		} else {
		    printf("copy_into_east_boundary_3D(): ");
		    printf("neighbour face %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imin; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest+1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src+1,j_src,k_src);
		dest = bp->get_cell(i_dest+2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == EAST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = bp_src->nnj - j - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    j_src = bp_src->nnj - k - 1; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 2 ) {
		    j_src = bp_src->nnj - j - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = j; k_src = bp_src->nnk - k - 1;
		} else {
		    printf("copy_into_east_boundary_3D(): ");
		    printf("neighbour face %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imax; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest+1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src-1,j_src,k_src);
		dest = bp->get_cell(i_dest+2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - j - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - k - 1; k_src = j;
		} else {
		    printf("copy_into_east_boundary_3D(): ");
		    printf("neighbour face %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmax; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest+1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src-1,k_src);
		dest = bp->get_cell(i_dest+2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - j - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 2 ) {
		    i_src = j; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; k_src = j;
		} else {
		    printf("copy_into_east_boundary_3D(): ");
		    printf("neighbour face %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmin; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest+1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src);
		dest = bp->get_cell(i_dest+2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == TOP ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - j - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; j_src = bp_src->nnj - j - 1;
		} else if ( orientation == 2 ) {
		    i_src = j; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; j_src = j;
		} else {
		    printf("copy_into_east_boundary_3D(): ");
		    printf("neighbour face %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmax; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest+1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src-1);
		dest = bp->get_cell(i_dest+2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = j; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; j_src = bp_src->nnj - j - 1;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - j - 1; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; j_src = j;
		} else {
		    printf("copy_into_east_boundary_3D(): ");
		    printf("neighbour face %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmin; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest+1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src+1);
		dest = bp->get_cell(i_dest+2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else {
	printf("copy_into_east_boundary_3D(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }

    return SUCCESS;
}   /* end copy_into_east_boundary_3D() */


/** \brief Do a straight copy of boundary data into the west-boundary
 *         ghost cells.
 */
int copy_into_west_boundary_3D(Block *bp, Block *bp_src, int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t i_dest, i_src, j_dest, j_src, k_dest, k_src, j, k;
    FV_Cell *src, *dest;
    int neighbour_faceId, orientation;

    i_dest = bp->imin;  /* index of the west-most plane of cells */
    neighbour_faceId = bp->bcp[WEST]->neighbour_face;
    orientation = bp->bcp[WEST]->neighbour_orientation;

    if ( neighbour_faceId == EAST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    j_src = bp_src->nnj - k - 1; k_src = j;
		} else if ( orientation == 2 ) {
		    j_src = bp_src->nnj - j - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = k; k_src = bp_src->nnk - j - 1;
		} else {
		    printf("copy_into_west_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imax; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest-1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src-1,j_src,k_src);
		dest = bp->get_cell(i_dest-2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == WEST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = bp_src->nnj - j - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 1 ) {
		    j_src = k; k_src = j;
		} else if ( orientation == 2 ) {
		    j_src = j; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = bp_src->nnj - k - 1; k_src = bp_src->nnk - j - 1;
		} else {
		    printf("copy_into_west_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imin; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest-1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src+1,j_src,k_src);
		dest = bp->get_cell(i_dest-2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - j - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; k_src = j;
		} else if ( orientation == 2 ) {
		    i_src = j; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - k - 1; k_src = bp_src->nnk - j - 1;
		} else {
		    printf("copy_into_west_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmax; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest-1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src-1,k_src);
		dest = bp->get_cell(i_dest-2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; k_src = j;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - j - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; k_src = bp_src->nnk - j - 1;
		} else {
		    printf("copy_into_west_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmin; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest-1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src);
		dest = bp->get_cell(i_dest-2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == TOP ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = j; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; j_src = j;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - j - 1; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; j_src = bp_src->nnj - j - 1;
		} else {
		    printf("copy_into_west_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmax; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest-1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src-1);
		dest = bp->get_cell(i_dest-2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - j - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; j_src = j;
		} else if ( orientation == 2 ) {
		    i_src = j; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - k - 1; j_src = bp_src->nnj - j - 1;
		} else {
		    printf("copy_into_west_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmin; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest-1,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src+1);
		dest = bp->get_cell(i_dest-2,j_dest,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* j loop */
    } else {
	printf("copy_into_west_boundary_3D(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }

    return SUCCESS;
}   /* end copy_into_west_boundary_3D() */


/** \brief Do a straight copy of boundary data into the north-boundary
 *         ghost cells.
 */
int copy_into_north_boundary_3D(Block *bp, Block *bp_src, int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t i_dest, i_src, j_dest, j_src, k_dest, k_src, i, k;
    FV_Cell *src, *dest;
    int neighbour_faceId, orientation;

    j_dest = bp->jmax;  /* index of the north-most plane of cells */
    neighbour_faceId = bp->bcp[NORTH]->neighbour_face;
    orientation = bp->bcp[NORTH]->neighbour_orientation;

    if ( neighbour_faceId == SOUTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - k - 1; k_src = i;
		} else {
		    printf("copy_into_north_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmin; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src);
		dest = bp->get_cell(i_dest,j_dest+2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; k_src = i;
		} else if ( orientation == 2 ) {
		    i_src = i; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - k - 1; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_north_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmax; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src-1,k_src);
		dest = bp->get_cell(i_dest,j_dest+2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == EAST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    j_src = bp_src->nnj - k - 1; k_src = i;
		} else if ( orientation == 2 ) {
		    j_src = bp_src->nnj - i - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = k; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_north_boundary_3D(): ");
		    printf("neighbour face %d invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imax; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src-1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == WEST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = bp_src->nnj - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    j_src = k; k_src = i;
		} else if ( orientation == 2 ) {
		    j_src = i; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = bp_src->nnj - k - 1; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_north_boundary_3D(): ");
		    printf("neighbour face %d invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imin; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src+1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == TOP ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = i; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; j_src = i;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; j_src = bp_src->nnj - i - 1;
		} else {
		    printf("copy_into_north_boundary_3D(): ");
		    printf("neighbour face %d invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmax; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src-1);
		dest = bp->get_cell(i_dest,j_dest+2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; j_src = i;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - k - 1; j_src = bp_src->nnj - i - 1;
		} else {
		    printf("copy_into_north_boundary_3D(): ");
		    printf("neighbour face %d invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmin; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest+1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src+1);
		dest = bp->get_cell(i_dest,j_dest+2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else {
	printf("copy_into_north_boundary_3D(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }

    return SUCCESS;
}   /* end copy_into_north_boundary_3D() */


/** \brief Do a straight copy of boundary data into the south-boundary
 *         ghost cells.
 */
int copy_into_south_boundary_3D(Block *bp, Block *bp_src, int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t i_dest, i_src, j_dest, j_src, k_dest, k_src, i, k;
    FV_Cell *src, *dest;
    int neighbour_faceId, orientation;

    j_dest = bp->jmin;  /* index of the south-most plane of cells */
    neighbour_faceId = bp->bcp[SOUTH]->neighbour_face;
    orientation = bp->bcp[SOUTH]->neighbour_orientation;

    if ( neighbour_faceId == NORTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; k_src = i;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_south_boundary_3D(): ");
		    printf(" neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmax; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src-1,k_src);
		dest = bp->get_cell(i_dest,j_dest-2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = i; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; k_src = i;
		} else {
		    printf("copy_into_south_boundary_3D(): ");
		    printf(" neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmin; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src);
		dest = bp->get_cell(i_dest,j_dest-2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == EAST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = bp_src->nnj - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    j_src = bp_src->nnj - k - 1; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    j_src = i; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = k; k_src = i;
		} else {
		    printf("copy_into_south_boundary_3D(): ");
		    printf(" neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imax; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src-1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == WEST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    j_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    j_src = k; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    j_src = bp_src->nnj - i - 1; k_src = bp_src->nnk - k - 1;
		} else if ( orientation == 3 ) {
		    j_src = bp_src->nnj - k - 1; k_src = i;
		} else {
		    printf("copy_into_south_boundary_3D(): ");
		    printf(" neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imin; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src+1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == TOP ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - k - 1; j_src = bp_src->nnj - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = k; j_src = i;
		} else {
		    printf("copy_into_south_boundary_3D(): ");
		    printf(" neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmax; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src-1);
		dest = bp->get_cell(i_dest,j_dest-2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    i_src = i; j_src = k;
		} else if ( orientation == 1 ) {
		    i_src = k; j_src = bp_src->nnj - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; j_src = bp_src->nnj - k - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - k - 1; j_src = i;
		} else {
		    printf("copy_into_south_boundary_3D(): ");
		    printf(" neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmin; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest-1,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src+1);
		dest = bp->get_cell(i_dest,j_dest-2,k_dest);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* k loop */
	}   /* i loop */
    } else {
	printf("copy_into_south_boundary_3D(): invlaid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }

    return SUCCESS;
}   /* end copy_into_south_boundary_3D() */


/** \brief Do a straight copy of boundary data into the top-boundary
 *         ghost cells.
 */
int copy_into_top_boundary_3D(Block *bp, Block *bp_src, int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t i_dest, i_src, j_dest, j_src, k_dest, k_src, i, j;
    FV_Cell *src, *dest;
    int neighbour_faceId, orientation;

    k_dest = bp->kmax;  /* index of the top-most plane of cells */
    neighbour_faceId = bp->bcp[TOP]->neighbour_face;
    orientation = bp->bcp[TOP]->neighbour_orientation;

    if ( neighbour_faceId == BOTTOM ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = i; j_src = j;
		} else if ( orientation == 1 ) {
		    i_src = j; j_src = bp_src->nnj - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; j_src = bp_src->nnj - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - j - 1; j_src = i;
		} else {
		    printf("copy_into_top_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmin; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src+1);
		dest = bp->get_cell(i_dest,j_dest,k_dest+2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == TOP ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; j_src = j;
		} else if ( orientation == 1 ) {
		    i_src = j; j_src = bp_src->nnj - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = i; j_src = bp_src->nnj - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = j; j_src = i;
		} else {
		    printf("copy_into_top_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmax; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src-1);
		dest = bp->get_cell(i_dest,j_dest,k_dest+2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = i; k_src = j;
		} else if ( orientation == 1 ) {
		    i_src = j; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - j - 1; k_src = i;
		} else {
		    printf("copy_into_top_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmax; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src-1,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - j - 1; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    i_src = i; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = j; k_src = i;
		} else {
		    printf("copy_into_top_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmin; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == EAST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    j_src = bp_src->nnj - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    j_src = bp->nnj - j - 1; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    j_src = i; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    j_src = j; k_src = i;
		} else {
		    printf("copy_into_top_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imax; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src-1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == WEST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    j_src = i; k_src = j;
		} else if ( orientation == 1 ) {
		    j_src = j; k_src = bp_src->nnk - i - 1;
		} else if ( orientation == 2 ) {
		    j_src = bp_src->nnj - i - 1; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    j_src = bp_src->nnj - j - 1; k_src = i;
		} else {
		    printf("copy_into_top_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imin; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src+1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest+2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else {
	printf("copy_into_top_boundary_3D(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }

    return SUCCESS;
}   /* end copy_into_top_boundary_3D() */


/** \brief Do a straight copy of boundary data into the bottom-boundary
 *         ghost cells.
 */
int copy_into_bottom_boundary_3D(Block *bp, Block *bp_src, int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t i_dest, i_src, j_dest, j_src, k_dest, k_src, i, j;
    FV_Cell *src, *dest;
    int neighbour_faceId, orientation;

    k_dest = bp->kmin;  /* index of the bottom-most plane of cells */
    neighbour_faceId = bp->bcp[BOTTOM]->neighbour_face;
    orientation = bp->bcp[BOTTOM]->neighbour_orientation;

    if ( neighbour_faceId == TOP ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = i; j_src = j;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - j - 1; j_src = i;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; j_src = bp_src->nnj - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = j; j_src = bp_src->nnj - i - 1;
		} else {
		    printf("copy_into_bottom_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmax; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src-1);
		dest = bp->get_cell(i_dest,j_dest,k_dest-2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; j_src = j;
		} else if ( orientation == 1 ) {
		    i_src = j; j_src = i;
		} else if ( orientation == 2 ) {
		    i_src = i; j_src = bp_src->nnj - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - j - 1; j_src = bp_src->nnj - i - 1;
		} else {
		    printf("copy_into_bottom_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		k_src = bp_src->kmin; 
		i_src += bp_src->imin; j_src += bp_src->jmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src,k_src+1);
		dest = bp->get_cell(i_dest,j_dest,k_dest-2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    i_src = j; k_src = i;
		} else if ( orientation == 2 ) {
		    i_src = i; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = bp_src->nni - j - 1; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_bottom_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmax; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src-1,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    i_src = bp_src->nni - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    i_src = bp_src->nni - j - 1; k_src = i;
		} else if ( orientation == 2 ) {
		    i_src = bp_src->nni - i - 1; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    i_src = j; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_bottom_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		j_src = bp_src->jmin; 
		i_src += bp_src->imin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src,j_src+1,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == EAST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    j_src = i; k_src = j;
		} else if ( orientation == 1 ) {
		    j_src = bp_src->nnj - j - 1; k_src = i;
		} else if ( orientation == 2 ) {
		    j_src = bp_src->nnj - i - 1; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    j_src = j; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_bottom_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imax; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src-1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else if ( neighbour_faceId == WEST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (i = 0; i < bp->nni; ++i) {
		i_dest = i + bp->imin;
		if ( orientation == 0 ) {
		    j_src = bp_src->nnj - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    j_src = j; k_src = i;
		} else if ( orientation == 2 ) {
		    j_src = i; k_src = bp_src->nnk - j - 1;
		} else if ( orientation == 3 ) {
		    j_src = bp_src->nnj - j - 1; k_src = bp_src->nnk - i - 1;
		} else {
		    printf("copy_into_bottom_boundary_3D(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		i_src = bp_src->imin; 
		j_src += bp_src->jmin; k_src += bp_src->kmin;
		src = bp_src->get_cell(i_src,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-1);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		src = bp_src->get_cell(i_src+1,j_src,k_src);
		dest = bp->get_cell(i_dest,j_dest,k_dest-2);
		dest->copy_values_from(*src, type_of_copy, gtl);
		dest->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    }   /* i loop */
	}   /* j loop */
    } else {
	printf("copy_into_bottom_boundary_3D(): unimplemented neighbour face %d.\n",
	       neighbour_faceId);
	exit(NOT_IMPLEMENTED_ERROR);
    }

    return SUCCESS;
}   /* end copy_into_bottom_boundary_3D() */


/*--------------------------------------------------------------------*/
/*                  Copy-via-buffer exchange functions.               */ 
/*--------------------------------------------------------------------*/


/** \brief Copy data into the send buffer from the appropriate
 *         boundary of the current block.
 *
 * \note Cells are always copied into the buffer in standard order.
 *       We will put the smarts for sorting out the orientation into
 *       the copy_from_receive_buffer_3D() function and its children.
 *
 * \param bd           : pointer to the source block
 * \param bndry        : value specifying the source block boundary
 * \param type_of_copy : sometimes we want to copy cell geometry,
 *                       other times the flow data
 */
int copy_into_send_buffer_3D(Block *bp, int bndry, int type_of_copy, 
			     double *send_buffer, size_t gtl)
{
    size_t i, j, k, i_src, j_src, k_src;
    size_t ib; /* position of cell in linear buffer. */
    size_t nv; /* number of double values transferred per cell */
    FV_Cell *cell;
    double *bufp;

    if ( bp->bcp[bndry]->neighbour_block < 0 ) {
	/* There is no neighbour, so we have no work to do. */
	return 0;
    }

    nv = number_of_values_in_cell_copy(type_of_copy);

    if ( bndry == NORTH ) {
	j_src = bp->jmax;
	for (i = 0; i < bp->nni; ++i) {
	    i_src = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_src = k + bp->kmin;
		/* Copy cell closest to boundary first. */
		cell = bp->get_cell(i_src,j_src,k_src);
		ib = (bp->nnk * i) + k;
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
		/* Copy inner cell second. */
		cell = bp->get_cell(i_src,j_src-1,k_src);
		ib += (bp->nnk * bp->nni);
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	    } /* k loop */
	} /* i loop */
    } else if ( bndry == SOUTH ) {
	j_src = bp->jmin;
	for (i = 0; i < bp->nni; ++i) {
	    i_src = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_src = k + bp->kmin;
		/* Copy cell closest to boundary first. */
		cell = bp->get_cell(i_src,j_src,k_src);
		ib = (bp->nnk * i) + k;
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
		/* Copy inner cell second. */
		cell = bp->get_cell(i_src,j_src+1,k_src);
		ib += (bp->nnk * bp->nni);
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	    } /* k loop */
	} /* i loop */
    } else if ( bndry == EAST ) {
	i_src = bp->imax;
	for (j = 0; j < bp->nnj; ++j) {
	    j_src = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_src = k + bp->kmin;
		/* Copy cell closest to boundary first. */
		cell = bp->get_cell(i_src,j_src,k_src);
		ib = (bp->nnk * j) + k;
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
		/* Copy inner cell second. */
		cell = bp->get_cell(i_src-1,j_src,k_src);
		ib += (bp->nnk * bp->nnj);
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	    } /* k loop */
	} /* j loop */
    } else if ( bndry == WEST ) {
	i_src = bp->imin;
	for (j = 0; j < bp->nnj; ++j) {
	    j_src = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_src = k + bp->kmin;
		/* Copy cell closest to boundary first. */
		cell = bp->get_cell(i_src,j_src,k_src);
		ib = (bp->nnk * j) + k;
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
		/* Copy inner cell second. */
		cell = bp->get_cell(i_src+1,j_src,k_src);
		ib += (bp->nnk * bp->nnj);
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	    } /* k loop */
	} /* j loop */
    } else if ( bndry == TOP ) {
	k_src = bp->kmax;
	for (i = 0; i < bp->nni; ++i) {
	    i_src = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_src = j + bp->jmin;
		/* Copy cell closest to boundary first. */
		cell = bp->get_cell(i_src,j_src,k_src);
		ib = (bp->nnj * i) + j;
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
		/* Copy inner cell second. */
		cell = bp->get_cell(i_src,j_src,k_src-1);
		ib += (bp->nnj * bp->nni);
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	    } /* j loop */
	} /* i loop */
    } else if ( bndry == BOTTOM ) {
	k_src = bp->kmin;
	for (i = 0; i < bp->nni; ++i) {
	    i_src = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_src = j + bp->jmin;
		/* Copy cell closest to boundary first. */
		cell = bp->get_cell(i_src,j_src,k_src);
		ib = (bp->nnj * i) + j;
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
		/* Copy inner cell second. */
		cell = bp->get_cell(i_src,j_src,k_src+1);
		ib += (bp->nnj * bp->nni);
		bufp = &(send_buffer[ib * nv]);
		cell->copy_values_to_buffer(bufp, type_of_copy, gtl);
	    } /* j loop */
	} /* i loop */
    } else {
        printf("copy_into_send_buffer(): invalid boundary %d\n", bndry);
        exit(VALUE_ERROR);
    }   /* end if */

    return SUCCESS;
}   /* end copy_into_send_buffer_3D() */


/** \brief Copy data from the receive buffer into the specified boundary.
 *         of the current block.
 *
 * This coordinating function exists so that we can have a set of smaller
 * functions contain the actual exchange code.  Having it all in one big 
 * function was a bit overwhelming.
 *
 * \param bp           : pointer to the target block
 * \param bndry        : value specifying the block boundary
 * \param type_of_copy :
 * \param receive_buffer:
 */
int copy_from_receive_buffer_3D(Block *bp, int bndry, int type_of_copy,
				double *receive_buffer, size_t gtl)
{
    if ( bp->bcp[bndry]->neighbour_block < 0 ) {
	/* There is no neighbour block so we have nothing to do. */
	return SUCCESS;
    }

    if ( bndry == NORTH ) {
	copy_from_receive_buffer_to_north(bp, type_of_copy, receive_buffer, gtl);
    } else if ( bndry == SOUTH ) {
	copy_from_receive_buffer_to_south(bp, type_of_copy, receive_buffer, gtl);
    } else if ( bndry == EAST ) {
	copy_from_receive_buffer_to_east(bp, type_of_copy, receive_buffer, gtl);
    } else if ( bndry == WEST ) {
	copy_from_receive_buffer_to_west(bp, type_of_copy, receive_buffer, gtl);
    } else if ( bndry == TOP ) {
	copy_from_receive_buffer_to_top(bp, type_of_copy, receive_buffer, gtl);
    } else if ( bndry == BOTTOM ) {
	copy_from_receive_buffer_to_bottom(bp, type_of_copy, receive_buffer, gtl);
    } else {
        cerr << "\ncopy_from_receive_buffer_3D(): invalid boundary\n" << endl;
        return FAILURE;
    }   /* end if */

    return SUCCESS;
}   /* end copy_from_receive_buffer_3D() */


/** \brief Copy data from the receive buffer into the NORTH boundary 
 *         of the current block.
 *
 * \param bp           : pointer to the target block
 * \param type_of_copy :
 * \param receive_buffer:
 */
int copy_from_receive_buffer_to_north(Block *bp, int type_of_copy,
				      double *receive_buffer, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int bndry, neighbour_faceId, orientation;
    size_t i, k, i_dest, j_dest, k_dest;
    size_t i_src, j_src, k_src, nni_src, nnj_src, nnk_src;
    size_t ib; /* position of cell in linear buffer. */
    size_t nv; /* number of double values transferred per cell */
    FV_Cell *cell;
    double *bufp;

    bndry = NORTH;
    neighbour_faceId = bp->bcp[bndry]->neighbour_face;
    orientation = bp->bcp[bndry]->neighbour_orientation;
    nv = number_of_values_in_cell_copy(type_of_copy);

    j_dest = bp->jmax;
    if ( neighbour_faceId == SOUTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = k; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = nni_src - i - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = nni_src - k - 1; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_north(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest+1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest+2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = nni_src - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = k; k_src = i;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = nni_src - k - 1; k_src = nnk_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_north(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest+1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest+2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == EAST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = nnj_src - k - 1; k_src = i;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = nnj_src - i - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = k; k_src = nnk_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_north(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest+1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest+2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == WEST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = nnj_src - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = k; k_src = i;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = nnj_src - k - 1; k_src = nnk_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_north(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest+1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest+2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == TOP ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = nni_src - k - 1; j_src = i;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = nni_src - i - 1; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = k; j_src = nnj_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_north(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest+1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest+2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = nni_src - i - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = k; j_src = i;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = nni_src - k - 1; j_src = nnj_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_north(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest+1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest+2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else {
	printf("copy_from_receive_buffer_to_north(): invalid neighbour face %d\n.",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }
    return SUCCESS;
}   /* end copy_from_receive_buffer_to_north() */


/** \brief Copy data from the receive buffer into the SOUTH boundary 
 *         of the current block.
 *
 * \param bp           : pointer to the target block
 * \param type_of_copy :
 * \param receive_buffer:
 */
int copy_from_receive_buffer_to_south(Block *bp, int type_of_copy,
				      double *receive_buffer, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int bndry, neighbour_faceId, orientation;
    size_t i, k, i_dest, j_dest, k_dest;
    size_t i_src, j_src, k_src, nni_src, nnj_src, nnk_src;
    size_t ib; /* position of cell in linear buffer. */
    size_t nv; /* number of double values transferred per cell */
    FV_Cell *cell;
    double *bufp;

    bndry = SOUTH;
    neighbour_faceId = bp->bcp[bndry]->neighbour_face;
    orientation = bp->bcp[bndry]->neighbour_orientation;
    nv = number_of_values_in_cell_copy(type_of_copy);

    j_dest = bp->jmin;
    if ( neighbour_faceId == NORTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = k; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = nni_src - i - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = nni_src - k - 1; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_south(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest-1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest-2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = nni_src - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = nni_src - k - 1; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nni;
		    i_src = k; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_south(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnk;
		    i_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest-1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest-2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == EAST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = nnj_src - i - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = nnj_src - k - 1; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = k; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_south(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest-1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest-2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == WEST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = k; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = nnj_src - i - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nni;
		    j_src = nnj_src - k - 1; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_south(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnk;
		    j_src = i; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest-1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest-2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == TOP ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = nni_src - i - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = nni_src - k - 1; j_src = nnj_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = k; j_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_south(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest-1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest-2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = k; j_src = nnj_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = nni_src - i - 1; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nni;
		    i_src = nni_src - k - 1; j_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_south(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnk;
		    i_src = i; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest-1,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest-2,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* i loop */
    } else {
	printf("copy_from_receive_buffer_to_south(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }

    return 0;
}   /* end copy_from_receive_buffer_to_south() */


/** \brief Copy data from the receive buffer into the EAST boundary 
 *         of the current block.
 *
 * \param bp           : pointer to the target block
 * \param type_of_copy :
 * \param receive_buffer:
 */
int copy_from_receive_buffer_to_east(Block *bp, int type_of_copy,
				     double *receive_buffer, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int bndry, neighbour_faceId, orientation;
    size_t j, k, i_dest, j_dest, k_dest;
    size_t i_src, j_src, k_src, nni_src, nnj_src, nnk_src;
    size_t ib; /* position of cell in linear buffer. */
    size_t nv; /* number of double values transferred per cell */
    FV_Cell *cell;
    double *bufp;

    bndry = EAST;
    neighbour_faceId = bp->bcp[bndry]->neighbour_face;
    orientation = bp->bcp[bndry]->neighbour_orientation;
    nv = number_of_values_in_cell_copy(type_of_copy);

    i_dest = bp->imax;
    if ( neighbour_faceId == WEST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = k; k_src = nnk_src - j - 1;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = nnj_src - j - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = nnj_src - k - 1; k_src = j;
		} else {
		    printf("copy_from_receive_buffer_to_east(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest+1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest+2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == EAST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = nnj_src - j - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = k; k_src = nnk_src - j - 1;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = k; k_src = j;
		} else {
		    printf("copy_from_receive_buffer_to_east(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest+1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest+2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = k; k_src = nnk_src - j - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = nni_src - j - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = nni_src - k - 1; k_src = j;
		} else {
		    printf("copy_from_receive_buffer_to_east(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest+1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest+2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = nni_src - j - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = nni_src - k - 1; k_src = nnk_src - j - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = k; k_src = j;
		} else {
		    printf("copy_from_receive_buffer_to_east(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest+1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest+2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == TOP ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = nni_src - j - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = nni_src - k - 1; j_src = nnj_src - j - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = k; j_src = j;
		} else {
		    printf("copy_from_receive_buffer_to_east(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest+1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest+2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = k; j_src = nnj_src - j - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = nni_src - j - 1; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = nni_src - k - 1; j_src = j;
		} else {
		    printf("copy_from_receive_buffer_to_east(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest+1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest+2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else {
	printf("copy_from_receive_buffer_to_east(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }
    return SUCCESS;
}   /* end copy_from_receive_buffer_to_east() */


/** \brief Copy data from the receive buffer into the WEST boundary 
 *         of the current block.
 *
 * \param bp           : pointer to the target block
 * \param type_of_copy :
 * \param receive_buffer:
 */
int copy_from_receive_buffer_to_west(Block *bp, int type_of_copy,
				     double *receive_buffer, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int bndry, neighbour_faceId, orientation;
    size_t j, k, i_dest, j_dest, k_dest;
    size_t i_src, j_src, k_src, nni_src, nnj_src, nnk_src;
    size_t ib; /* position of cell in linear buffer. */
    size_t nv; /* number of double values transferred per cell */
    FV_Cell *cell;
    double *bufp;

    bndry = WEST;
    neighbour_faceId = bp->bcp[bndry]->neighbour_face;
    orientation = bp->bcp[bndry]->neighbour_orientation;
    nv = number_of_values_in_cell_copy(type_of_copy);

    i_dest = bp->imin;
    if ( neighbour_faceId == EAST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = nnj_src - k - 1; k_src = j;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = nnj_src - j - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = k; k_src = nnk_src - j - 1;
		} else {
		    printf("copy_from_receive_buffer_to_west(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest-1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest-2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == WEST ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = nnj_src - j - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = k; k_src = j;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnk; nnk_src = bp->nnj;
		    j_src = nnj_src - k - 1; k_src = nnk_src - j - 1;
		} else {
		    printf("copy_from_receive_buffer_to_west(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nnj; nnk_src = bp->nnk;
		    j_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest-1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest-2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = nni_src - j - 1; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = k; k_src = j;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = nni_src - k - 1; k_src = nnk_src - j - 1;
		} else {
		    printf("copy_from_receive_buffer_to_west(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest-1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest-2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = nni_src - k - 1; k_src = j;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = nni_src - j - 1; k_src = nnk_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnk_src = bp->nnj;
		    i_src = k; k_src = nnk_src - j - 1;
		} else {
		    printf("copy_from_receive_buffer_to_west(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnk_src = bp->nnk;
		    i_src = j; k_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest-1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest-2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == TOP ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = nni_src - k - 1; j_src = j;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = nni_src - j - 1; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = k; j_src = j;
		} else {
		    printf("copy_from_receive_buffer_to_west(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest-1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest-2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (j = 0; j < bp->nnj; ++j) {
	    j_dest = j + bp->jmin;
	    for (k = 0; k < bp->nnk; ++k) {
		k_dest = k + bp->kmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = nni_src - j - 1; j_src = k;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = k; j_src = j;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = nnj_src - k - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnk; nnj_src = bp->nnj;
		    i_src = nni_src - k - 1; j_src = nnj_src - j - 1;
		} else {
		    printf("copy_from_receive_buffer_to_west(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nnj; nnj_src = bp->nnk;
		    i_src = j; j_src = k;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest-1,j_dest,k_dest);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest-2,j_dest,k_dest);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* k loop */
	} /* j loop */
    } else {
	printf("copy_from_receive_buffer_to_west(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }
    return SUCCESS;
}   /* end copy_from_receive_buffer_to_west() */


/** \brief Copy data from the receive buffer into the TOP boundary 
 *         of the current block.
 *
 * \param bp           : pointer to the target block
 * \param type_of_copy :
 * \param receive_buffer:
 */
int copy_from_receive_buffer_to_top(Block *bp, int type_of_copy,
				    double *receive_buffer, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int bndry, neighbour_faceId, orientation;
    size_t i, j, i_dest, j_dest, k_dest;
    size_t i_src, j_src, k_src, nni_src, nnj_src, nnk_src;
    size_t ib; /* position of cell in linear buffer. */
    size_t nv; /* number of double values transferred per cell */
    FV_Cell *cell;
    double *bufp;

    bndry = TOP;
    neighbour_faceId = bp->bcp[bndry]->neighbour_face;
    orientation = bp->bcp[bndry]->neighbour_orientation;
    nv = number_of_values_in_cell_copy(type_of_copy);

    k_dest = bp->kmax;
    if ( neighbour_faceId == BOTTOM ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = i; j_src = j;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = j; j_src = nnj_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = nni_src - i - 1; j_src = nnj_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = nni_src - j - 1; j_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_top(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest+1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest+2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == TOP ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = nni_src - i - 1; j_src = j;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = nni_src - j - 1; j_src = nnj_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = nni_src - i - 1; j_src = nnj_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = j; j_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_top(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest+1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest+2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = j; // PJ 2014-09-24, fix for Wilson.
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = j; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = nni_src - i - 1; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = nni_src - j - 1; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_top(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest+1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest+2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = nni_src - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = nni_src - j - 1; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = j; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_top(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest+1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest+2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == EAST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = nnj_src - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = nnj_src - j - 1; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = j; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_top(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest+1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest+2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == WEST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = j;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = j; k_src = nnk_src - i - 1;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = nnj_src - i - 1; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = nnj_src - j - 1; k_src = i;
		} else {
		    printf("copy_from_receive_buffer_to_top(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest+1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest+2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else {
	printf("copy_from_receive_buffer_to_top(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }
    return SUCCESS;
}   /* end copy_from_receive_buffer_to_top() */


/** \brief Copy data from the receive buffer into the BOTTOM boundary 
 *         of the current block.
 *
 * \param bp           : pointer to the target block
 * \param type_of_copy :
 * \param receive_buffer:
 */
int copy_from_receive_buffer_to_bottom(Block *bp, int type_of_copy,
				       double *receive_buffer, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    int bndry, neighbour_faceId, orientation;
    size_t i, j, i_dest, j_dest, k_dest;
    size_t i_src, j_src, k_src, nni_src, nnj_src, nnk_src;
    size_t ib; /* position of cell in linear buffer. */
    size_t nv; /* number of double values transferred per cell */
    FV_Cell *cell;
    double *bufp;

    bndry = BOTTOM;
    neighbour_faceId = bp->bcp[bndry]->neighbour_face;
    orientation = bp->bcp[bndry]->neighbour_orientation;
    nv = number_of_values_in_cell_copy(type_of_copy);

    k_dest = bp->kmin;
    if ( neighbour_faceId == TOP ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = i; j_src = j;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = nni_src - j - 1; j_src = i;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = nni_src - i - 1; j_src = nnj_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = j; j_src = nnj_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_bottom(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest-1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest-2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == BOTTOM ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = nni_src - i - 1; j_src = j;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = j; j_src = i;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = i; j_src = nnj_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnj_src = bp->nni;
		    i_src = nni_src - j - 1; j_src = nnj_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_bottom(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnj_src = bp->nnj;
		    i_src = i; j_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnj_src * i_src) + j_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest-1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnj_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest-2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == NORTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = nni_src - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = j; k_src = i;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = nni_src - j - 1; k_src = nnk_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_bottom(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest-1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest-2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == SOUTH ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = j;
		} else if ( orientation == 1 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = nni_src - j - 1; k_src = i;
		} else if ( orientation == 2 ) {
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = nni_src - i - 1; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nni_src = bp->nnj; nnk_src = bp->nni;
		    i_src = j; k_src = nnk_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_bottom(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nni_src = bp->nni; nnk_src = bp->nnj;
		    i_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * i_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest-1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nni_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest-2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == EAST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = j;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = nnj_src - j - 1; k_src = i;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = nnj_src - i - 1; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = j; k_src = nnk_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_bottom(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest-1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest-2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else if ( neighbour_faceId == WEST ) {
	for (i = 0; i < bp->nni; ++i) {
	    i_dest = i + bp->imin;
	    for (j = 0; j < bp->nnj; ++j) {
		j_dest = j + bp->jmin;
		if ( orientation == 0 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = nnj_src - i - 1; k_src = j;
		} else if ( orientation == 1 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = j; k_src = i;
		} else if ( orientation == 2 ) {
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = nnk_src - j - 1;
		} else if ( orientation == 3 ) {
		    nnj_src = bp->nnj; nnk_src = bp->nni;
		    j_src = nnj_src - j - 1; k_src = nnk_src - i - 1;
		} else {
		    printf("copy_from_receive_buffer_to_bottom(): ");
		    printf("neighbour_faceId %d, invalid orientation %d.\n",
			   neighbour_faceId, orientation);
		    nnj_src = bp->nni; nnk_src = bp->nnj;
		    j_src = i; k_src = j;
		    exit(VALUE_ERROR);
		}
		/* Fill the first line of ghost cells. */
		ib = (nnk_src * j_src) + k_src;
		cell = bp->get_cell(i_dest,j_dest,k_dest-1);
		bufp = &(receive_buffer[ib * nv]);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
		/* Fill the second line of ghost cells. */
		ib += (nnk_src * nnj_src);
		bufp = &(receive_buffer[ib * nv]);
		cell = bp->get_cell(i_dest,j_dest,k_dest-2);
		cell->copy_values_from_buffer(bufp, type_of_copy, gtl);
		cell->encode_conserved(gtl, 0, bp->omegaz, with_k_omega);
	    } /* j loop */
	} /* i loop */
    } else {
	printf("copy_from_receive_buffer_to_bottom(): invalid neighbour face %d.\n",
	       neighbour_faceId);
	exit(VALUE_ERROR);
    }
    return SUCCESS;
}   /* end copy_from_receive_buffer_to_bottom() */

/*--------------------------------------------------------------------*/
