// bc/slip_wall.d
//
// Solid-wall which allows the fluid to slip along it with no shear stress.
// Peter J. 2014-07-26

import std.conv;

import bc;
import block;
import sblock;

class SlipWallBC: BoundaryCondition {

    this(int id, int boundary, double emissivity=0.0) 
    {
	blk_id = id;
	which_boundary = boundary;
	type_code = BCCode.slip_wall;
	is_wall = true;
	ghost_cell_data_available = true;
	sets_conv_flux_directly = false;
	sets_visc_flux_directly = false;
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
