// exch_mapped_cell_shmem.cxx
// Copy mapped-cell data in a shared-memory context.
// 
// PJ, 08-Mar-2014

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cell.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"

int copy_mapped_cell_data_via_shmem(int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    bool with_k_omega = (G.turbulence_model == TM_K_OMEGA);
    size_t src_blk, src_i, src_j, src_k, ghost_cell_index;
    size_t dest_i, dest_j, dest_k;
    FV_Cell * src;
    FV_Cell * dest;

    for ( Block *bdp : G.my_blocks ) {
	int number_faces = (G.dimensions == 3 ? 6: 4);
	for ( int iface = 0; iface < number_faces; ++iface ) {
	    if ( bdp->bcp[iface]->type_code == MAPPED_CELL ) {
		// Let's get the ghost-cell data.
		switch ( iface ) {
		case NORTH:
		case SOUTH:
		    // cout << "NORTH/SOUTH boundary" << endl;
		    ghost_cell_index = 0;
		    for ( size_t i = 0; i < bdp->nni; ++i ) {
			for ( size_t k = 0; k < bdp->nnk; ++k ) {
			    for ( size_t ghost_cell_count = 1; ghost_cell_count < 2; ++ghost_cell_count ) {
				src_blk = bdp->bcp[iface]->mapped_cells[ghost_cell_index][0];
				src_i = bdp->bcp[iface]->mapped_cells[ghost_cell_index][1];
				src_j = bdp->bcp[iface]->mapped_cells[ghost_cell_index][2];
				src_k = bdp->bcp[iface]->mapped_cells[ghost_cell_index][3];
				// cout << "ghost_cell_index=" << ghost_cell_index 
				//      << " src_blk=" << src_blk << " src_i=" << src_i
				//      << " src_j=" << src_j << " src_k=" << src_k << endl;
				src = get_block_data_ptr(src_blk)->get_cell(src_i, src_j, src_k);
				dest_i = i + bdp->imin;
				dest_j = iface == NORTH ? (bdp->jmax + ghost_cell_count) : (bdp->jmin - ghost_cell_count);
				dest_k = k + bdp->kmin;
				dest = bdp->get_cell(dest_i, dest_j, dest_k);
				dest->copy_values_from(*src, type_of_copy, gtl);
				dest->encode_conserved(gtl, 0, bdp->omegaz, with_k_omega);
			    }
			} // end for k
		    } // end for i
		    break;
		case EAST:
		case WEST:
		    for ( size_t j = 0; j < bdp->nnj; ++j ) {
			for ( size_t k = 0; k < bdp->nnk; ++k ) {
			    for ( size_t ghost_cell_count = 1; ghost_cell_count < 2; ++ghost_cell_count ) {
				src_blk = bdp->bcp[iface]->mapped_cells[ghost_cell_index][0];
				src_i = bdp->bcp[iface]->mapped_cells[ghost_cell_index][1];
				src_j = bdp->bcp[iface]->mapped_cells[ghost_cell_index][2];
				src_k = bdp->bcp[iface]->mapped_cells[ghost_cell_index][3];
				src = get_block_data_ptr(src_blk)->get_cell(src_i, src_j, src_k);
				dest_i = iface == EAST ? (bdp->imax + ghost_cell_count) : (bdp->imin - ghost_cell_count);
				dest_j = j + bdp->jmin;
				dest_k = k + bdp->kmin;
				dest = bdp->get_cell(dest_i, dest_j, dest_k);
				dest->copy_values_from(*src, type_of_copy, gtl);
				dest->encode_conserved(gtl, 0, bdp->omegaz, with_k_omega);
			    }
			} // end for k
		    } // end for i
		    break;
		case TOP:
		case BOTTOM:
		    for ( size_t i = 0; i < bdp->nni; ++i ) {
			for ( size_t j = 0; j < bdp->nnj; ++j ) {
			    for ( size_t ghost_cell_count = 1; ghost_cell_count < 2; ++ghost_cell_count ) {
				src_blk = bdp->bcp[iface]->mapped_cells[ghost_cell_index][0];
				src_i = bdp->bcp[iface]->mapped_cells[ghost_cell_index][1];
				src_j = bdp->bcp[iface]->mapped_cells[ghost_cell_index][2];
				src_k = bdp->bcp[iface]->mapped_cells[ghost_cell_index][3];
				src = get_block_data_ptr(src_blk)->get_cell(src_i, src_j, src_k);
				dest_i = i + bdp->imin;
				dest_j = j + bdp->jmin;
				dest_k = iface == TOP ? (bdp->kmax + ghost_cell_count) : (bdp->kmin - ghost_cell_count);
				dest = bdp->get_cell(dest_i, dest_j, dest_k);
				dest->copy_values_from(*src, type_of_copy, gtl);
				dest->encode_conserved(gtl, 0, bdp->omegaz, with_k_omega);
			    }
			} // end for j
		    } // end for i
		} // end switch iface
	    } // end if type_code == MAPPED_CELL
	} // end for face
    } // end for Block
    return SUCCESS;
} // end copy_mapped_cell_data_via_shmem()
