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
    size_t i, j, k;
    FV_Cell *cL1, *cL0, *cR0, *cR1, *cR2;
    FV_Interface *IFace;
    Gas_model *gmodel = get_gas_model_ptr();
    // Maybe these two FlowState objects should be in the Block object
    // and initialised there so that we don't thrash the memory so much.
    FlowState Lft(gmodel);
    FlowState Rght(gmodel);
    
    // ifi interfaces are East-facing interfaces.
    for (k = kmin; k <= kmax; ++k) {
	for (j = jmin; j <= jmax; ++j) {
	    for (i = imin; i <= imax+1; ++i) {
		IFace = get_ifi(i,j,k);
		cL1 = get_cell(i-2,j,k);
		cL0 = get_cell(i-1,j,k);
		cR0 = get_cell(i,j,k);
		cR1 = get_cell(i+1,j,k);
		// Interpolate LEFT and RIGHT interface states from the cell-center properties.
		if ( ( i == imin ) && ( bcp[WEST]->type_code == SHOCK_FITTING_IN ) ) {
		    FV_Interface *IFaceL = get_ifi(i-1,j,k);
		    FV_Interface *IFaceR = get_ifi(i,j,k);
		    compute_boundary_flux(IFaceL, IFaceR, omegaz);
                    // We're shock-fitting and we're on the shock boundary.
		    // Retain the inflow defined flux at the boundary by doing nothing here.
		} else {
		    // Either we're not shock fitting or we're not on the shock boundary.
		    if ( get_adaptive_reconstruction_flag() == 1 ) {
			mach_weighted_one_d_interp(*cL1, *cL0, *cR0, *cR1,
				      cL1->iLength, cL0->iLength, cR0->iLength, cR1->iLength, Lft, Rght);
		    } else {
			one_d_interp(*cL1, *cL0, *cR0, *cR1,
				     cL1->iLength, cL0->iLength, cR0->iLength, cR1->iLength, Lft, Rght);
		    }
		    // Use special interpolation scheme for first interface after shock.
		    if ( ( i == imin+1 ) && ( bcp[WEST]->type_code == SHOCK_FITTING_IN ) ) {
			cR2 = get_cell(i+2,j,k);
        		one_d_interior_interp(*cL0, *cR0, *cR1, *cR2, 
					      cL0->iLength, cR0->iLength, cR1->iLength, cR2->iLength,
					      Lft, Rght);
		    }
		    // Save u, v, w, T for the viscous flux calculation by making a local average.
		    // The values for u, v and T may be updated subsequently by the interface-flux function.
		    IFace->fs->average_values_from(Lft, Rght, get_diffusion_flag()==1);
		    // Finally, the flux calculation.
		    if ( (i == imin && bcp[WEST]->use_udf_flux()) ||
			 (i == imax+1 && bcp[EAST]->use_udf_flux()) ) {
			// Retain the user-defined flux at the boundary
			// by doing nothing here.
		    } else {
			compute_interface_flux(Lft, Rght, *IFace, omegaz);
		    } // end if
		} // end if
		// DEBUGGING for ablating BC
#               if 0
		if ( i == imax+1 ) {
		    cout << "i = " << i << endl;
		    cout << "j = " << j << endl;
		    cout << "pos = " << IFace->pos.str() << endl;
		    cout << "total mass flux = " << IFace->F->mass << endl;
		    for ( size_t isp=0; isp<IFace->F->massf.size(); ++isp ) {
		    	cout << "species[" << isp << "] mass flux = " 
			     << IFace->F->mass * IFace->F->massf[isp] << endl;
		    }
		}
#               endif
	    } // i loop
	} // j loop
    } // for k

    // ifj interfaces are North-facing interfaces.
    for (k = kmin; k <= kmax; ++k) {
	for (i = imin; i <= imax; ++i) {
	    for (j = jmin; j <= jmax+1; ++j) {
		IFace = get_ifj(i,j,k);
		cL1 = get_cell(i,j-2,k);
		cL0 = get_cell(i,j-1,k);
		cR0 = get_cell(i,j,k);
		cR1 = get_cell(i,j+1,k);
		// Interpolate LEFT and RIGHT interface states from the cell-center properties.
		if ( get_adaptive_reconstruction_flag() == 1 ) {
		    // Use Mach number weighted 
	            mach_weighted_one_d_interp(*cL1, *cL0, *cR0, *cR1, 
	                          cL1->jLength, cL0->jLength, cR0->jLength, cR1->jLength, 
	                          Lft, Rght);
	        } else {
	            one_d_interp(*cL1, *cL0, *cR0, *cR1, 
				 cL1->jLength, cL0->jLength, cR0->jLength, cR1->jLength, 
				 Lft, Rght);
	        }
		// Save u, v, w, T for the viscous flux calculation by making a local average.
		// The values for u, v and T may be updated subsequently by the interface-flux function.
	        IFace->fs->average_values_from(Lft, Rght, get_diffusion_flag()==1);
		// Finally, the flux calculation.
	        if ( (j == jmin && bcp[SOUTH]->use_udf_flux()) ||
	             (j == jmax+1 && bcp[NORTH]->use_udf_flux()) ) {
	            // Retain the user-defined flux at the boundary
	            // by doing nothing here.
	        } else {
	            compute_interface_flux(Lft, Rght, *IFace, omegaz);
	        } // end if
	    } // j loop
	} // i loop
    } // for k
    
    if ( dimensions == 2 ) return SUCCESS;
    
    // ifk interfaces are TOP-facing interfaces.
    for (i = imin; i <= imax; ++i) {
	for (j = jmin; j <= jmax; ++j) {
	    for (k = kmin; k <= kmax+1; ++k) {
		IFace = get_ifk(i,j,k);
		cL1 = get_cell(i,j,k-2);
		cL0 = get_cell(i,j,k-1);
		cR0 = get_cell(i,j,k);
		cR1 = get_cell(i,j,k+1);
		// Interpolate LEFT and RIGHT interface states from the cell-center properties.
		if ( get_adaptive_reconstruction_flag() == 1 ) {
		    mach_weighted_one_d_interp(*cL1, *cL0, *cR0, *cR1, 
				  cL1->kLength, cL0->kLength, cR0->kLength, cR1->kLength, 
				  Lft, Rght);
		} else {
		    one_d_interp(*cL1, *cL0, *cR0, *cR1, 
				 cL1->kLength, cL0->kLength, cR0->kLength, cR1->kLength, 
				 Lft, Rght);
		}
		// Save u, v, w, T for the viscous flux calculation by making a local average.
		// The values for u, v and T may be updated subsequently by the interface-flux function.
		IFace->fs->average_values_from(Lft, Rght, get_diffusion_flag()==1);
		// Finally, the flux calculation.
		if ( (k == kmin && bcp[BOTTOM]->use_udf_flux()) ||
		     (k == kmax+1 && bcp[TOP]->use_udf_flux()) ) {
		    // Retain the user-defined flux at the boundary
		    // by doing nothing here.
		} else {
		    compute_interface_flux(Lft, Rght, *IFace, omegaz);
		} // end if
	    } // for k 
	} // j loop
    } // i loop
    
    return SUCCESS;
}   // end of inviscid_flux()
