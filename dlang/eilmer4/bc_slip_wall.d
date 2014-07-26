// bc_slip_wall.d
//
// Solid-wall which allows the fluid to slip along it with no shear stress.
// Peter J. 2014-07-26

import bc;
import block;
import sblock;

class SlipWallBC: BoundaryCondition {

    this(ref SBlock blk, int which_boundary, double emissivity=0.0) 
    {
	type_code = BCCode.slip_wall;
	is_wall = true;
	this.emissivity = emissivity;
	this.which_boundary = which_boundary;
	blk.bc[which_boundary] = this;
    }

    // Let the base class implementations do the work.
    // apply_convective
    // apply_viscous

} // end class SlipWallBC
