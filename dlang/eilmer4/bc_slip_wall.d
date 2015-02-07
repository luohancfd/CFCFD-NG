// bc_slip_wall.d
//
// Solid-wall which allows the fluid to slip along it with no shear stress.
// Peter J. 2014-07-26

import std.conv;

import bc;
import block;
import sblock;

class SlipWallBC: BoundaryCondition {

    this(ref SBlock blk_, int which_boundary_, double emissivity=0.0) 
    {
	blk = blk_;
	which_boundary = which_boundary_;
	type_code = BCCode.slip_wall;
	is_wall = true;
	this.emissivity = emissivity;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "SlipWallBC(";
	repr ~= "emissivity=" ~ to!string(emissivity);
	repr ~= ")";
	return to!string(repr);
    }

    // Let the base class implementations do the work.
    // apply_convective
    // apply_viscous

} // end class SlipWallBC
