// bc_adiabatic_wall.d
//
// Solid-wall with no-slip velocity and zero-temperature-gradient.
// Peter J. 2014-07-26

import fvcore;
import flowstate;
import fvinterface;
import fvcell;
import bc;
import block;
import sblock;
import bc_menter_correction;

class AdiabaticWallBC: BoundaryCondition {

    this(ref SBlock blk, int which_boundary, double emissivity=0.0) 
    {
	type_code = BCCode.adiabatic_wall;
	is_wall = true;
	this.emissivity = emissivity;
	this.which_boundary = which_boundary;
	blk.bc[which_boundary] = this;
    }

    // Let the base class implementation do the work.
    // apply_convective -- mirror velocity

    override void apply_viscous(double t)
    // Notes:
    // Menter's slightly-rough-surface boundary condition as described
    // in Wilcox 2006 text, eqn 7.36.
    // We assume that the y2 in eqn 7.16 is the same as
    // the height of our finite-volume cell.
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;

	final switch ( which_boundary ) {
	case north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[north];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[east];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[south];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[west];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[top];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for i
	    break;
	case bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[bottom];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply_viscous()
} // end class SlipWallBC
