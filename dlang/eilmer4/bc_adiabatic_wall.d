// bc_adiabatic_wall.d
//
// Solid-wall with no-slip velocity and zero-temperature-gradient.
// Peter J. 2014-07-26

import std.conv;

import fvcore;
import flowstate;
import fvinterface;
import fvcell;
import bc;
import block;
import sblock;
import bc_menter_correction;
import globaldata;

class AdiabaticWallBC: BoundaryCondition {

    this(int id, int boundary, double emissivity=0.0) 
    {
	blk_id = id;
	which_boundary = boundary;
	type_code = BCCode.adiabatic_wall;
	is_wall = true;
	ghost_cell_data_available = true;
	sets_conv_flux_directly = false;
	sets_visc_flux_directly = false;
	this.emissivity = emissivity;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "AdiabaticWallBC(";
	repr ~= "emissivity=" ~ to!string(emissivity);
	repr ~= ")";
	return to!string(repr);
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
	auto blk = allBlocks[blk_id];

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
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
