// bc_menter_correction.cxx
/// \brief Apply Menter boundary correction to the cells near solid walls.
///
/// Menter's slightly-rough-surface boundary condition is described
/// in Wilcox's 2006 text, eqn 7.36.
/// For low-resolution grids, the k-omega model is reported to over-estimate
/// the magnitude of omega, well out into the boundary layer so,
/// to get reasonable values for omega close to the wall, we propagate
/// the 1/y**2 form of the omega data out a few cells from the wall.
///
/// PJ, October 2007

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"


double ideal_omega_at_wall(FV_Cell *cell)
{
    Gas_data *wall_gas = cell->cell_at_nearest_wall->fs->gas;
    double d0 = cell->half_cell_width_at_wall;
    return 400.0 * wall_gas->mu / wall_gas->rho / (d0 * d0);
}

double ideal_omega(FV_Cell *cell)
{
    double d0 = cell->half_cell_width_at_wall;
    double d = cell->distance_to_nearest_wall;
    return ideal_omega_at_wall(cell) * (d0 * d0) / ((d0 + d) * (d0 + d));
}

int apply_menter_boundary_correction(Block &bd, size_t ftl)
{
    size_t i, j, k;
    size_t layer_depth;
    size_t nominal_layer_depth=6; // Nominal number of cells over which we 
                               // will correct the omega value.
    FV_Cell *cell;
    FV_Interface *IFace;
    global_data &G = *get_global_data_ptr();

    // Step 1: Apply the ideal solution for omega to the layers of cells against viscous walls.
    //         Do not apply Menter boundary correction if we are in a laminar region.

    // NORTH boundary
    if ( bd.bcp[NORTH]->is_wall() && bd.bcp[NORTH]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the j-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnj/2, nominal_layer_depth);
	for ( i = bd.imin; i <= bd.imax; ++i ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    j =  bd.jmax - indx;
		    cell = bd.get_cell(i,j,k);
		    if ( cell->in_turbulent_zone == 0 ) continue;
		    cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
		    cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
	        } // j-loop
            } // k-loop
	} // i-loop
    }

    // SOUTH boundary
    if ( bd.bcp[SOUTH]->is_wall() && bd.bcp[SOUTH]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the j-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnj/2, nominal_layer_depth);
	for ( i = bd.imin; i <= bd.imax; ++i ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    j = bd.jmin + indx;
		    cell = bd.get_cell(i,j,k);
		    if ( cell->in_turbulent_zone == 0 ) continue;
		    cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
		    cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		}  // j-loop
	    } // k-loop
	} // i-loop
    }

    // EAST boundary
    if ( bd.bcp[EAST]->is_wall() && bd.bcp[EAST]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the i-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nni/2, nominal_layer_depth);
	for ( j = bd.jmin; j <= bd.jmax; ++j ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    i = bd.imax - indx;
		    cell = bd.get_cell(i,j,k);
		    if ( cell->in_turbulent_zone == 0 ) continue;
		    cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
		    cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
	        } // i-loop
            } // k-loop
	} // j-loop
    }

    // WEST boundary
    if ( bd.bcp[WEST]->is_wall() && bd.bcp[WEST]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the i-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nni/2, nominal_layer_depth);
	for (j = bd.jmin; j <= bd.jmax; ++j) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    i = bd.imin + indx;
		    cell = bd.get_cell(i,j,k);
		    if ( cell->in_turbulent_zone == 0 ) continue;
		    cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
		    cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
	        } // i-loop
            } // k-loop
	} // j-loop
    }

    if ( G.dimensions == 3 ) {
	// TOP boundary
	if ( bd.bcp[TOP]->is_wall() && bd.bcp[TOP]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the k-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnk/2, nominal_layer_depth);
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		    for ( size_t indx = 0; indx < layer_depth; ++indx ) {
			k = bd.kmax - indx;
			cell = bd.get_cell(i,j,k);
			if ( cell->in_turbulent_zone == 0 ) continue;
			cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		    } // k-loop
		} // j-loop
	    } // i-loop
	}
        
	// BOTTOM boundary
	if ( bd.bcp[BOTTOM]->is_wall() && bd.bcp[BOTTOM]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the k-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnk/2, nominal_layer_depth);
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		    for ( size_t indx = 0; indx < layer_depth; ++indx ) {
			k = bd.kmin + indx;
			cell = bd.get_cell(i,j,k);
			if ( cell->in_turbulent_zone == 0 ) continue;
			cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		    }  // k-loop
		} // j-loop
	    } // i-loop
	}
    } // end if G.dimensions == 3

    // Step 2: After doing all of the viscous walls,
    //         we need to go around and tidy up the faces on inviscid walls
    //         so that we have consistent omega values for cells that have
    //         been corrected at other boundaries.

    // NORTH boundary
    if ( bd.bcp[NORTH]->type_code == SLIP_WALL || 
	 bd.bcp[NORTH]->type_code == EXTRAPOLATE_OUT ) {
	for ( i = bd.imin; i <= bd.imax; ++i ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	        j = bd.jmax;
                cell = bd.get_cell(i,j,k);
	        IFace = cell->iface[NORTH];
		IFace->fs->omega = cell->fs->omega;
            } // k-loop
	} // i-loop
    }

    // SOUTH boundary
    if ( bd.bcp[SOUTH]->type_code == SLIP_WALL ||
	 bd.bcp[SOUTH]->type_code == EXTRAPOLATE_OUT ) {
	for ( i = bd.imin; i <= bd.imax; ++i ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
		j = bd.jmin;
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		IFace->fs->omega = cell->fs->omega;
	    } // k-loop
	} // i-loop
    }

    // EAST boundary
    if ( bd.bcp[EAST]->type_code == SLIP_WALL ||
	 bd.bcp[EAST]->type_code == EXTRAPOLATE_OUT ) {
	for ( j = bd.jmin; j <= bd.jmax; ++j ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	        i = bd.imax;
	        cell = bd.get_cell(i,j,k);
	        IFace = cell->iface[EAST];
		IFace->fs->omega = cell->fs->omega;
            } // k-loop
	} // j-loop
    }

    // WEST boundary
    if ( bd.bcp[WEST]->type_code == SLIP_WALL ||
	 bd.bcp[WEST]->type_code == EXTRAPOLATE_OUT ) {
	for (j = bd.jmin; j <= bd.jmax; ++j) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
	        i = bd.imin;
	        cell = bd.get_cell(i,j,k);
	        IFace = cell->iface[WEST];
		IFace->fs->omega = cell->fs->omega;
            } // k-loop
	} // j-loop
    }

    if ( G.dimensions == 3 ) {
	// TOP boundary
	if ( bd.bcp[TOP]->type_code == SLIP_WALL ||
	     bd.bcp[TOP]->type_code == EXTRAPOLATE_OUT ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		    k = bd.kmax;
		    cell = bd.get_cell(i,j,k);
		    IFace = cell->iface[TOP];
		    IFace->fs->omega = cell->fs->omega;
		} // j-loop
	    } // i-loop
	}
        
	// BOTTOM boundary
	if ( bd.bcp[BOTTOM]->type_code == SLIP_WALL ||
	     bd.bcp[BOTTOM]->type_code == EXTRAPOLATE_OUT ) {
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		    k = bd.kmin;
		    cell = bd.get_cell(i,j,k);
		    IFace = cell->iface[BOTTOM];
		    IFace->fs->omega = cell->fs->omega;
		} // j-loop
	    } // i-loop
	}
    } // end if G.dimensions == 3

    return SUCCESS;
} // end of apply_menter_boundary_correction()

int apply_wilson_omega_correction(Block &bd, size_t ftl)
{
    // bc_wilson_omega_correction.cxx
    /// \brief Apply Wilson's omega smoothing correction to the cells near solid walls.
    ///
    /// The near-wall oscillations seems to be from the amplification of numerical
    /// noise by the large magnitudes of omega values near wall. This correction essentially
    /// averages the omega values using averaged density values over the width of the block 
    /// (in the cross-stream direction) and then feeds the averaged omega values back to
    /// the cells. 
    ///
    /// I guess this is a good enough solution for our checkboard fluctuations at this moment
    /// despite it being an extremely brute-force and crude method. Also, it would be a rather
    /// justified correction since we are already fudging the omega values using Menter's
    /// method currently.
    ///
    /// Instead of averaging over a particular direction (we will obviously face issues when
    /// we have so many different cases), we shall try averaging the 8 cells surrounding a
    /// particular cell together with that cell (so average over 9 cells). This method of
    /// smoothing should be more general than our earlier fix.
    ///
    /// Out of convenience, we average the edge and corner cells the same way too (so it 
    /// uses values in the ghost cells for averaging).
    ///
    /// Only works for NORTH, SOUTH, TOP and EAST walls now.
    ///
    /// Wilson Chan, January 2011

    size_t i, j, k;
    size_t layer_depth;
    size_t nominal_layer_depth=10; // Nominal number of cells over which we 
                                // will correct the omega value.
    double omegaAvg, rhoAvg;
    FV_Cell *cell, *cN, *cE, *cS, *cW, *cNE, *cSE, *cSW, *cNW;
    global_data &G = *get_global_data_ptr();

    // NORTH boundary
    if ( bd.bcp[NORTH]->is_wall() && bd.bcp[NORTH]->type_code != SLIP_WALL ) {
        // Use the smaller value of either the number of cells in 
        // the j-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnj, nominal_layer_depth);
        for ( size_t indx = 0; indx < layer_depth; ++indx ) {
            j =  bd.jmax - indx;
            // Am too lazy to think of corner and edge cases, so
            // we will just average over the internal cells first
//            for ( k = bd.kmin+1; k <= bd.kmax-1; ++k ) {
//                for ( i = bd.imin+1; i <= bd.imax-1; ++i ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
                for ( i = bd.imin; i <= bd.imax; ++i ) {
                    // Step 1 : Average all density and omega values in 
                    //          surrounding cells in j-k plane.
                    cell = bd.get_cell(i,j,k);
                    cN = bd.get_cell(i,j,k+1);
                    cE = bd.get_cell(i+1,j,k);
                    cS = bd.get_cell(i,j,k-1);
                    cW = bd.get_cell(i-1,j,k);
                    cNE = bd.get_cell(i+1,j,k+1);
                    cSE = bd.get_cell(i+1,j,k-1);
                    cSW = bd.get_cell(i-1,j,k-1);
                    cNW = bd.get_cell(i-1,j,k+1);
                    if ( cell->in_turbulent_zone == 0 ) continue;
                    omegaAvg = 1.0 / 9.0 * (cell->fs->omega + cN->fs->omega + cE->fs->omega +
                               cS->fs->omega + cW->fs->omega + cNE->fs->omega +
                               cSE->fs->omega + cSW->fs->omega + cNW->fs->omega);
                    rhoAvg = 1.0 / 9.0 * (cell->fs->gas->rho + cN->fs->gas->rho +
                             cE->fs->gas->rho + cS->fs->gas->rho + cW->fs->gas->rho +
                             cNE->fs->gas->rho + cSE->fs->gas->rho +
                             cSW->fs->gas->rho + cNW->fs->gas->rho);
                    // Step 2 : Replace cell and conserved quantity values with
                    //          averaged (smoothed) values of omega.
                    cell->fs->omega = omegaAvg;
                    cell->U[ftl]->omega = rhoAvg * omegaAvg;
                }  // j-loop
            } // i-loop
        } // k-loop
    }

    // SOUTH boundary
    if ( bd.bcp[SOUTH]->is_wall() && bd.bcp[SOUTH]->type_code != SLIP_WALL ) {
        // Use the smaller value of either the number of cells in 
        // the j-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnj, nominal_layer_depth);
        for ( size_t indx = 0; indx < layer_depth; ++indx ) {
            j = bd.jmin + indx;
            // Am too lazy to think of corner and edge cases, so
            // we will just average over the internal cells first
//            for ( k = bd.kmin+1; k <= bd.kmax-1; ++k ) {
//                for ( i = bd.imin+1; i <= bd.imax-1; ++i ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
                for ( i = bd.imin; i <= bd.imax; ++i ) {
                    // Step 1 : Average all density and omega values in 
                    //          surrounding cells in j-k plane.
		    cell = bd.get_cell(i,j,k);
                    cN = bd.get_cell(i,j,k+1);
                    cE = bd.get_cell(i+1,j,k);
                    cS = bd.get_cell(i,j,k-1);
                    cW = bd.get_cell(i-1,j,k);
                    cNE = bd.get_cell(i+1,j,k+1);
                    cSE = bd.get_cell(i+1,j,k-1);
                    cSW = bd.get_cell(i-1,j,k-1);
                    cNW = bd.get_cell(i-1,j,k+1);
		    if ( cell->in_turbulent_zone == 0 ) continue;
                    omegaAvg = 1.0 / 9.0 * (cell->fs->omega + cN->fs->omega + cE->fs->omega +
                               cS->fs->omega + cW->fs->omega + cNE->fs->omega +
                               cSE->fs->omega + cSW->fs->omega + cNW->fs->omega);
                    rhoAvg = 1.0 / 9.0 * (cell->fs->gas->rho + cN->fs->gas->rho +
                             cE->fs->gas->rho + cS->fs->gas->rho + cW->fs->gas->rho +
                             cNE->fs->gas->rho + cSE->fs->gas->rho +
                             cSW->fs->gas->rho + cNW->fs->gas->rho);
                    // Step 2 : Replace cell and conserved quantity values with
                    //          averaged (smoothed) values of omega.
                    cell->fs->omega = omegaAvg;
                    cell->U[ftl]->omega = rhoAvg * omegaAvg;
		}  // j-loop
	    } // i-loop
	} // k-loop
    }

    // EAST boundary
    if ( bd.bcp[EAST]->is_wall() && bd.bcp[EAST]->type_code != SLIP_WALL ) {
        // Make a copy of the current flow state and conserved quantities
        // so that we have the actual values to average over and not average
        // over averaged values.
        // FIX-ME - HOW DO WE COPY ALL OMEGA AND RHO VALUES IN THE J-K PLANE??
        //   OBVIOUSLY IF WE DON'T COPY THEN WE WILL BE APPLYING A SMEARING FILTER!!!
        //tempCell->fs->copy_values_from(*(cell->fs)); 
        //tempCell->U[ftl]->copy_values_from(*(cell->U));
        // Use the smaller value of the number of cells in 
        // the i-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nni, nominal_layer_depth);
        for ( size_t indx = 0; indx < layer_depth; ++indx ) {
            i = bd.imax - indx;
            // Am too lazy to think of corner and edge cases, so
            // we will just average over the internal cells first
//            for ( k = bd.kmin+1; k <= bd.kmax-1; ++k ) {
//                for ( j = bd.jmin+1; j <= bd.jmax-1; ++j ) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
                for ( j = bd.jmin; j <= bd.jmax; ++j ) {
                    // Step 1 : Average all density and omega values in 
                    //          surrounding cells in j-k plane.
		    cell = bd.get_cell(i,j,k);
                    cN = bd.get_cell(i,j,k+1);
                    cE = bd.get_cell(i,j+1,k);
                    cS = bd.get_cell(i,j,k-1);
                    cW = bd.get_cell(i,j-1,k);
                    cNE = bd.get_cell(i,j+1,k+1);
                    cSE = bd.get_cell(i,j+1,k-1);
                    cSW = bd.get_cell(i,j-1,k-1);
                    cNW = bd.get_cell(i,j-1,k+1);
		    if ( cell->in_turbulent_zone == 0 ) continue;
                    omegaAvg = 1.0 / 9.0 * (cell->fs->omega + cN->fs->omega + cE->fs->omega + 
                               cS->fs->omega + cW->fs->omega + cNE->fs->omega + 
                               cSE->fs->omega + cSW->fs->omega + cNW->fs->omega);
                    rhoAvg = 1.0 / 9.0 * (cell->fs->gas->rho + cN->fs->gas->rho + 
                             cE->fs->gas->rho + cS->fs->gas->rho + cW->fs->gas->rho + 
                             cNE->fs->gas->rho + cSE->fs->gas->rho + 
                             cSW->fs->gas->rho + cNW->fs->gas->rho);
                    // Step 2 : Replace cell and conserved quantity values with
                    //          averaged (smoothed) values of omega.
                    cell->fs->omega = omegaAvg;
                    cell->U[ftl]->omega = rhoAvg * omegaAvg;
                } // j-loop
            } // k-loop
	} // i-loop
    }

/*  // WEST boundary
    if ( bd.bcp[WEST]->is_wall() && bd.bcp[WEST]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the i-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nni, nominal_layer_depth);
	for (j = bd.jmin; j <= bd.jmax; ++j) {
            for ( k = bd.kmin; k <= bd.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    i = bd.imin + indx;
		    cell = bd.get_cell(i,j,k);
		    if ( cell->in_turbulent_zone == 0 ) continue;
		    cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
		    cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
	        } // i-loop
            } // k-loop
	} // j-loop
    }
*/
    if ( G.dimensions == 3 ) {
	// TOP boundary
	if ( bd.bcp[TOP]->is_wall() && bd.bcp[TOP]->type_code != SLIP_WALL ) {
        // Use the smaller value of either the number of cells in 
        // the k-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnk, nominal_layer_depth);
        for ( size_t indx = 0; indx < layer_depth; ++indx ) {
            k = bd.kmax - indx;
            // Am too lazy to think of corner and edge cases, so
            // we will just average over the internal cells first
//            for ( j = bd.jmin+1; j <= bd.jmax-1; ++j ) {
//                for ( i = bd.imin+1; i <= bd.imax-1; ++i ) {
            for ( j = bd.jmin; j <= bd.jmax; ++j ) {
                for ( i = bd.imin; i <= bd.imax; ++i ) {
                    // Step 1 : Average all density and omega values in 
                    //          surrounding cells in j-k plane.
                    cell = bd.get_cell(i,j,k);
                    cN = bd.get_cell(i,j+1,k);
                    cE = bd.get_cell(i+1,j,k);
                    cS = bd.get_cell(i,j-1,k);
                    cW = bd.get_cell(i-1,j,k);
                    cNE = bd.get_cell(i+1,j+1,k);
                    cSE = bd.get_cell(i+1,j-1,k);
                    cSW = bd.get_cell(i-1,j-1,k);
                    cNW = bd.get_cell(i-1,j+1,k);
                    if ( cell->in_turbulent_zone == 0 ) continue;
                    omegaAvg = 1.0 / 9.0 * (cell->fs->omega + cN->fs->omega + cE->fs->omega +
                               cS->fs->omega + cW->fs->omega + cNE->fs->omega +
                               cSE->fs->omega + cSW->fs->omega + cNW->fs->omega);
                    rhoAvg = 1.0 / 9.0 * (cell->fs->gas->rho + cN->fs->gas->rho +
                             cE->fs->gas->rho + cS->fs->gas->rho + cW->fs->gas->rho +
                             cNE->fs->gas->rho + cSE->fs->gas->rho +
                             cSW->fs->gas->rho + cNW->fs->gas->rho);
                    // Step 2 : Replace cell and conserved quantity values with
                    //          averaged (smoothed) values of omega.
                    cell->fs->omega = omegaAvg;
                    cell->U[ftl]->omega = rhoAvg * omegaAvg;
                } // i-loop
            } // j-loop
        } // k-loop
    }
        
/*	// BOTTOM boundary
	if ( bd.bcp[BOTTOM]->is_wall() && bd.bcp[BOTTOM]->type_code != SLIP_WALL ) {
        // Use the smaller value of either half the number of cells in 
        // the k-direction or the specified nominal_layer_depth. 
        layer_depth = min(bd.nnk, nominal_layer_depth);
	    for ( i = bd.imin; i <= bd.imax; ++i ) {
		for ( j = bd.jmin; j <= bd.jmax; ++j ) {
		    for ( size_t indx = 0; indx < layer_depth; ++indx ) {
			k = bd.kmin + indx;
			cell = bd.get_cell(i,j,k);
			if ( cell->in_turbulent_zone == 0 ) continue;
			cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		    }  // k-loop
		} // j-loop
	    } // i-loop
	}
*/
    } // end if G.dimensions == 3

    return SUCCESS;
} // end of apply_wilson_omega_correction()
