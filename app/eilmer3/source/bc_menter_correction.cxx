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
		    if ( cell->in_turbulent_zone ) {
			cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		    }
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
		    if ( cell->in_turbulent_zone ) {
			cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		    }
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
		    if ( cell->in_turbulent_zone ) {
			cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		    }
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
		    if ( cell->in_turbulent_zone ) {
			cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
		    }
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
			if ( cell->in_turbulent_zone ) {
			    cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			    cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
			}
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
			if ( cell->in_turbulent_zone ) {
			    cell->fs->omega = MINIMUM(ideal_omega(cell), cell->fs->omega);
			    cell->U[ftl]->omega = cell->fs->gas->rho * cell->fs->omega;
			}
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
