// bc_menter_correction.d
// Apply Menter boundary correction to the cells near solid walls.
//
// Menter's slightly-rough-surface boundary condition is described
// in Wilcox's 2006 text, eqn 7.36.
// For low-resolution grids, the k-omega model is reported to over-estimate
// the magnitude of omega, well out into the boundary layer so,
// to get reasonable values for omega close to the wall, we propagate
// the 1/y**2 form of the omega data out a few cells from the wall.
//
// PJ, October 2007
// 2014-07-26 ported from C++ code.

import std.math;
import std.algorithm;
import fvcore;
import fvinterface;
import fvcell;
import bc;
import sblock;
import globalconfig;

double ideal_omega_at_wall(in FVCell cell)
{
    auto wall_gas = cell.cell_at_nearest_wall.fs.gas;
    double d0 = cell.half_cell_width_at_wall;
    return 400.0 * wall_gas.mu / wall_gas.rho / (d0 * d0);
}

double ideal_omega(in FVCell cell)
{
    double d0 = cell.half_cell_width_at_wall;
    double d = cell.distance_to_nearest_wall;
    return ideal_omega_at_wall(cell) * (d0 * d0) / ((d0 + d) * (d0 + d));
}

void apply_menter_boundary_correction(ref SBlock blk, size_t ftl)
{
    size_t i, j, k;
    size_t layer_depth;
    size_t nominal_layer_depth=6; // Nominal number of cells over which we 
                               // will correct the omega value.
    FVCell cell;
    FVInterface IFace;

    // Step 1: Apply the ideal solution for omega to the layers of cells against viscous walls.
    //         Do not apply Menter boundary correction if we are in a laminar region.

    // north boundary
    if ( blk.bc[Face.north].is_wall && blk.bc[Face.north].type_code != BCCode.slip_wall ) {
        // Use the smaller value of either half the number of cells in 
        // the j-direction or the specified nominal_layer_depth. 
        layer_depth = min(blk.njcell/2, nominal_layer_depth);
	for ( i = blk.imin; i <= blk.imax; ++i ) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    j =  blk.jmax - indx;
		    cell = blk.get_cell(i,j,k);
		    if ( cell.in_turbulent_zone ) {
			cell.fs.omega = fmin(ideal_omega(cell), cell.fs.omega);
			cell.U[ftl].omega = cell.fs.gas.rho * cell.fs.omega;
		    }
	        } // j-loop
            } // k-loop
	} // i-loop
    }

    // south boundary
    if ( blk.bc[Face.south].is_wall && blk.bc[Face.south].type_code != BCCode.slip_wall ) {
        // Use the smaller value of either half the number of cells in 
        // the j-direction or the specified nominal_layer_depth. 
        layer_depth = min(blk.njcell/2, nominal_layer_depth);
	for ( i = blk.imin; i <= blk.imax; ++i ) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    j = blk.jmin + indx;
		    cell = blk.get_cell(i,j,k);
		    if ( cell.in_turbulent_zone ) {
			cell.fs.omega = fmin(ideal_omega(cell), cell.fs.omega);
			cell.U[ftl].omega = cell.fs.gas.rho * cell.fs.omega;
		    }
		}  // j-loop
	    } // k-loop
	} // i-loop
    }

    // east boundary
    if ( blk.bc[Face.east].is_wall && blk.bc[Face.east].type_code != BCCode.slip_wall ) {
        // Use the smaller value of either half the number of cells in 
        // the i-direction or the specified nominal_layer_depth. 
        layer_depth = min(blk.nicell/2, nominal_layer_depth);
	for ( j = blk.jmin; j <= blk.jmax; ++j ) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    i = blk.imax - indx;
		    cell = blk.get_cell(i,j,k);
		    if ( cell.in_turbulent_zone ) {
			cell.fs.omega = fmin(ideal_omega(cell), cell.fs.omega);
			cell.U[ftl].omega = cell.fs.gas.rho * cell.fs.omega;
		    }
	        } // i-loop
            } // k-loop
	} // j-loop
    }

    // west boundary
    if ( blk.bc[Face.west].is_wall && blk.bc[Face.west].type_code != BCCode.slip_wall ) {
        // Use the smaller value of either half the number of cells in 
        // the i-direction or the specified nominal_layer_depth. 
        layer_depth = min(blk.nicell/2, nominal_layer_depth);
	for (j = blk.jmin; j <= blk.jmax; ++j) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
		for ( size_t indx = 0; indx < layer_depth; ++indx ) {
		    i = blk.imin + indx;
		    cell = blk.get_cell(i,j,k);
		    if ( cell.in_turbulent_zone ) {
			cell.fs.omega = fmin(ideal_omega(cell), cell.fs.omega);
			cell.U[ftl].omega = cell.fs.gas.rho * cell.fs.omega;
		    }
	        } // i-loop
            } // k-loop
	} // j-loop
    }

    if ( GlobalConfig.dimensions == 3 ) {
	// top boundary
	if ( blk.bc[Face.top].is_wall && blk.bc[Face.top].type_code != BCCode.slip_wall ) {
        // Use the smaller value of either half the number of cells in 
        // the k-direction or the specified nominal_layer_depth. 
        layer_depth = min(blk.nkcell/2, nominal_layer_depth);
	    for ( i = blk.imin; i <= blk.imax; ++i ) {
		for ( j = blk.jmin; j <= blk.jmax; ++j ) {
		    for ( size_t indx = 0; indx < layer_depth; ++indx ) {
			k = blk.kmax - indx;
			cell = blk.get_cell(i,j,k);
			if ( cell.in_turbulent_zone ) {
			    cell.fs.omega = fmin(ideal_omega(cell), cell.fs.omega);
			    cell.U[ftl].omega = cell.fs.gas.rho * cell.fs.omega;
			}
		    } // k-loop
		} // j-loop
	    } // i-loop
	}
        
	// bottom boundary
	if ( blk.bc[Face.bottom].is_wall && blk.bc[Face.bottom].type_code != BCCode.slip_wall ) {
        // Use the smaller value of either half the number of cells in 
        // the k-direction or the specified nominal_layer_depth. 
        layer_depth = min(blk.nkcell/2, nominal_layer_depth);
	    for ( i = blk.imin; i <= blk.imax; ++i ) {
		for ( j = blk.jmin; j <= blk.jmax; ++j ) {
		    for ( size_t indx = 0; indx < layer_depth; ++indx ) {
			k = blk.kmin + indx;
			cell = blk.get_cell(i,j,k);
			if ( cell.in_turbulent_zone ) {
			    cell.fs.omega = fmin(ideal_omega(cell), cell.fs.omega);
			    cell.U[ftl].omega = cell.fs.gas.rho * cell.fs.omega;
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

    // north boundary
    if ( blk.bc[Face.north].type_code == BCCode.slip_wall || 
	 blk.bc[Face.north].type_code == BCCode.extrapolate_out ) {
	for ( i = blk.imin; i <= blk.imax; ++i ) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
	        j = blk.jmax;
                cell = blk.get_cell(i,j,k);
	        IFace = cell.iface[Face.north];
		IFace.fs.omega = cell.fs.omega;
            } // k-loop
	} // i-loop
    }

    // south boundary
    if ( blk.bc[Face.south].type_code == BCCode.slip_wall ||
	 blk.bc[Face.south].type_code == BCCode.extrapolate_out ) {
	for ( i = blk.imin; i <= blk.imax; ++i ) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
		j = blk.jmin;
		cell = blk.get_cell(i,j,k);
		IFace = cell.iface[Face.south];
		IFace.fs.omega = cell.fs.omega;
	    } // k-loop
	} // i-loop
    }

    // east boundary
    if ( blk.bc[Face.east].type_code == BCCode.slip_wall ||
	 blk.bc[Face.east].type_code == BCCode.extrapolate_out ) {
	for ( j = blk.jmin; j <= blk.jmax; ++j ) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
	        i = blk.imax;
	        cell = blk.get_cell(i,j,k);
	        IFace = cell.iface[Face.east];
		IFace.fs.omega = cell.fs.omega;
            } // k-loop
	} // j-loop
    }

    // west boundary
    if ( blk.bc[Face.west].type_code == BCCode.slip_wall ||
	 blk.bc[Face.west].type_code == BCCode.extrapolate_out ) {
	for (j = blk.jmin; j <= blk.jmax; ++j) {
            for ( k = blk.kmin; k <= blk.kmax; ++k ) {
	        i = blk.imin;
	        cell = blk.get_cell(i,j,k);
	        IFace = cell.iface[Face.west];
		IFace.fs.omega = cell.fs.omega;
            } // k-loop
	} // j-loop
    }

    if ( GlobalConfig.dimensions == 3 ) {
	// top boundary
	if ( blk.bc[Face.top].type_code == BCCode.slip_wall ||
	     blk.bc[Face.top].type_code == BCCode.extrapolate_out ) {
	    for ( i = blk.imin; i <= blk.imax; ++i ) {
		for ( j = blk.jmin; j <= blk.jmax; ++j ) {
		    k = blk.kmax;
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    IFace.fs.omega = cell.fs.omega;
		} // j-loop
	    } // i-loop
	}
        
	// bottom boundary
	if ( blk.bc[Face.bottom].type_code == BCCode.slip_wall ||
	     blk.bc[Face.bottom].type_code == BCCode.extrapolate_out ) {
	    for ( i = blk.imin; i <= blk.imax; ++i ) {
		for ( j = blk.jmin; j <= blk.jmax; ++j ) {
		    k = blk.kmin;
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    IFace.fs.omega = cell.fs.omega;
		} // j-loop
	    } // i-loop
	}
    } // end if G.dimensions == 3
} // end of apply_menter_boundary_correction()
