// bc_slip_wall.d
//
// Solid-wall which allows the fluid to slip along it with no shear stress.
// Peter J. 2014-07-26

import bc;

class SlipWallBC: BoundaryCondition {

    this() 
    {
	type_code = BCCode.slip_wall;
    }

    // Let the base class implementations do the work.
    // apply_convective
    // apply_viscous

} // end class SlipWallBC
