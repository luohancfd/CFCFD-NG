// bc_slip_wall.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_slip_wall.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

SlipWallBC::SlipWallBC( Block *bdp, int which_boundary )
    : BoundaryCondition(bdp, which_boundary, SLIP_WALL, "SlipWallBC",
			0, true, false, -1, -1) 
{}

SlipWallBC::SlipWallBC( const SlipWallBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face) 
{}

SlipWallBC::SlipWallBC()
    : BoundaryCondition(0, 0, SLIP_WALL, "SlipWallBC",
			0, true, false, -1, -1) 
{}

SlipWallBC & SlipWallBC::operator=(const SlipWallBC &bc)
{
    BoundaryCondition::operator=(bc);
    return *this;
}

SlipWallBC::~SlipWallBC() {}

