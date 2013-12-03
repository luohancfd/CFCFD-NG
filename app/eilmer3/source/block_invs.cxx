/** \file block_invs.cxx
 * \ingroup eilmer3
 * \brief Inviscid Flux functions
 *
 * \author PA Jacobs
 *
 * \version 03-Nov-96 : Updated multi-species code.
 * \version 26-Jun-97 : AUSMDV added
 * \version 12-Oct-97 : moved local_frame() and XY_frame() to cns_bc.c
 * \version 13-Oct-97 : put in temporary arrays for vectorisation (again)
 * \version 13-Feb-01 : Adaptive flux added.
 * \version 05-Aug-04 : Moved the generic flux calculation function to
 *                      ../../flux_calc/source/flux_calc.c
 *
 */

/*-----------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/geometry2/source/geom.hh"
#include "cell.hh"
#include "flux_calc.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "one_d_interp.hh"

/*-----------------------------------------------------------------*/

/* \brief  Given the cell-center values, compute the inviscid fluxes
 *         across the cell interfaces.
 *
 * The cells are treated one row (or column) at a time.
 * First, the left and right interface states are reconstructed
 * from the cell-centre data and then the fluxes across the
 * interfaces are calculated.
 */
int Block::inviscid_flux(size_t dimensions)
{
    global_data &G = *get_global_data_ptr();
    FV_Cell *cL1, *cL0, *cR0, *cR1;
    FV_Interface *IFace;
    Gas_model *gmodel = get_gas_model_ptr();
    // Maybe these two FlowState objects should be in the Block object
    // and initialised there so that we don't thrash the memory so much.
    FlowState Lft(gmodel);
    FlowState Rght(gmodel);
    
    // ifi interfaces are East-facing interfaces.
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t j = jmin; j <= jmax; ++j ) {
	    for ( size_t i = imin; i <= imax+1; ++i ) {
		IFace = get_ifi(i,j,k);
		cL1 = get_cell(i-2,j,k);
		cL0 = get_cell(i-1,j,k);
		cR0 = get_cell(i,j,k);
		cR1 = get_cell(i+1,j,k);
		if ( (bcp[WEST]->type_code == SHOCK_FITTING_IN) && (i == imin) ) {
		    // We're shock-fitting and we're on the shock boundary.
		    // Compute the flux at the boundary directly from the free-stream flow.
		    set_flux_vector_in_global_frame(*IFace, *(cL0->fs), this->omegaz);
		    // Second, save u, v, w, T for the viscous flux calculation.
		    IFace->fs->average_values_from(*(cL0->fs), *(cL0->fs), G.diffusion);
		} else {
		    // Compute the flux from data on either-side of the interface.
		    // First, interpolate LEFT and RIGHT interface states from cell-center properties.
		    if ( (i == imin+1) && (bcp[WEST]->ghost_cell_data_available == false) ) {
        		one_d_interp_right(*cL0, *cR0, *cR1, 
					   cL0->iLength, cR0->iLength, cR1->iLength,
					   Lft, Rght);
		    } else if ( (i == imax) && (bcp[EAST]->ghost_cell_data_available == false) ) {
        		one_d_interp_left(*cL1, *cL0, *cR0, 
					  cL1->iLength, cL0->iLength, cR0->iLength,
					  Lft, Rght);
		    } else { // General symmetric reconstruction.
			one_d_interp_both(*cL1, *cL0, *cR0, *cR1,
					  cL1->iLength, cL0->iLength, cR0->iLength, cR1->iLength,
					  Lft, Rght);
		    }
		    // Second, save u, v, w, T for the viscous flux calculation by making a local average.
		    // The values for u, v and T may be updated subsequently by the interface-flux function.
		    if ( (i == imin) && (bcp[WEST]->ghost_cell_data_available == false) ) {
			IFace->fs->average_values_from(Rght, Rght, G.diffusion);
		    } else if ( (i == imax+1) && (bcp[EAST]->ghost_cell_data_available == false) ) {
			IFace->fs->average_values_from(Lft, Lft, G.diffusion);
		    } else {
			IFace->fs->average_values_from(Lft, Rght, G.diffusion);
		    }
		    // Finally, the flux calculation itself.
		    if ( (i == imin && bcp[WEST]->sets_conv_flux()) ||
			 (i == imax+1 && bcp[EAST]->sets_conv_flux()) ) {
			// Retain the b.c. set flux at the boundary by doing nothing here.
		    } else {
			compute_interface_flux(Lft, Rght, *IFace, omegaz);
		    }
		} // end else Compute the flux from data on either-side of the interface.
	    } // i loop
	} // j loop
    } // for k

    // ifj interfaces are North-facing interfaces.
    for ( size_t k = kmin; k <= kmax; ++k ) {
	for ( size_t i = imin; i <= imax; ++i ) {
	    for ( size_t j = jmin; j <= jmax+1; ++j ) {
		IFace = get_ifj(i,j,k);
		cL1 = get_cell(i,j-2,k);
		cL0 = get_cell(i,j-1,k);
		cR0 = get_cell(i,j,k);
		cR1 = get_cell(i,j+1,k);
		// Interpolate LEFT and RIGHT interface states from the cell-center properties.
		if ( (j == jmin+1) && (bcp[SOUTH]->ghost_cell_data_available == false) ) {
		    one_d_interp_right(*cL0, *cR0, *cR1, 
				       cL0->jLength, cR0->jLength, cR1->jLength,
				       Lft, Rght);
		} else if ( (j == jmax) && (bcp[NORTH]->ghost_cell_data_available == false) ) {
		    one_d_interp_left(*cL1, *cL0, *cR0, 
				      cL1->jLength, cL0->jLength, cR0->jLength,
				      Lft, Rght);
		} else { // General symmetric reconstruction.
		    one_d_interp_both(*cL1, *cL0, *cR0, *cR1,
				      cL1->jLength, cL0->jLength, cR0->jLength, cR1->jLength,
				      Lft, Rght);
		}
		// Second, save u, v, w, T for the viscous flux calculation by making a local average.
		// The values for u, v and T may be updated subsequently by the interface-flux function.
		if ( (j == jmin) && (bcp[SOUTH]->ghost_cell_data_available == false) ) {
		    IFace->fs->average_values_from(Rght, Rght, G.diffusion);
		} else if ( (j == jmax+1) && (bcp[NORTH]->ghost_cell_data_available == false) ) {
		    IFace->fs->average_values_from(Lft, Lft, G.diffusion);
		} else {
		    IFace->fs->average_values_from(Lft, Rght, G.diffusion);
		}
		// Finally, the flux calculation.
	        if ( (j == jmin && bcp[SOUTH]->sets_conv_flux()) ||
	             (j == jmax+1 && bcp[NORTH]->sets_conv_flux()) ) {
	            // Retain the b.c. flux at the boundary by doing nothing here.
	        } else {
	            compute_interface_flux(Lft, Rght, *IFace, omegaz);
	        } // end if
	    } // j loop
	} // i loop
    } // for k
    
    if ( dimensions == 2 ) return SUCCESS;
    
    // ifk interfaces are TOP-facing interfaces.
    for ( size_t i = imin; i <= imax; ++i ) {
	for ( size_t j = jmin; j <= jmax; ++j ) {
	    for ( size_t k = kmin; k <= kmax+1; ++k ) {
		IFace = get_ifk(i,j,k);
		cL1 = get_cell(i,j,k-2);
		cL0 = get_cell(i,j,k-1);
		cR0 = get_cell(i,j,k);
		cR1 = get_cell(i,j,k+1);
		// Interpolate LEFT and RIGHT interface states from the cell-center properties.
		if ( (k == kmin+1) && (bcp[BOTTOM]->ghost_cell_data_available == false) ) {
		    one_d_interp_right(*cL0, *cR0, *cR1, 
				       cL0->kLength, cR0->kLength, cR1->kLength,
				       Lft, Rght);
		} else if ( (k == kmax) && (bcp[TOP]->ghost_cell_data_available == false) ) {
		    one_d_interp_left(*cL1, *cL0, *cR0, 
				      cL1->kLength, cL0->kLength, cR0->kLength,
				      Lft, Rght);
		} else { // General symmetric reconstruction.
		    one_d_interp_both(*cL1, *cL0, *cR0, *cR1,
				      cL1->kLength, cL0->kLength, cR0->kLength, cR1->kLength,
				      Lft, Rght);
		}
		// Second, save u, v, w, T for the viscous flux calculation by making a local average.
		// The values for u, v and T may be updated subsequently by the interface-flux function.
		if ( (k == kmin) && (bcp[BOTTOM]->ghost_cell_data_available == false) ) {
		    IFace->fs->average_values_from(Rght, Rght, G.diffusion);
		} else if ( (k == kmax+1) && (bcp[TOP]->ghost_cell_data_available == false) ) {
		    IFace->fs->average_values_from(Lft, Lft, G.diffusion);
		} else {
		    IFace->fs->average_values_from(Lft, Rght, G.diffusion);
		}
		// Finally, the flux calculation.
		if ( (k == kmin && bcp[BOTTOM]->sets_conv_flux()) ||
		     (k == kmax+1 && bcp[TOP]->sets_conv_flux()) ) {
		    // Retain the b.c. set flux at the boundary by doing nothing here.
		} else {
		    compute_interface_flux(Lft, Rght, *IFace, omegaz);
		} // end if
	    } // for k 
	} // j loop
    } // i loop
    
    return SUCCESS;
}   // end of inviscid_flux()
