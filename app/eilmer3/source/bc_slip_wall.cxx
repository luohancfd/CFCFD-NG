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

SlipWallBC::SlipWallBC(Block *bdp, int which_boundary, double emissivity)
    : BoundaryCondition(bdp, which_boundary, SLIP_WALL, emissivity)
{
    is_wall_flag = true;
}

SlipWallBC::SlipWallBC(const SlipWallBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.emissivity)
{
    is_wall_flag = bc.is_wall_flag;
}

SlipWallBC::SlipWallBC()
    : BoundaryCondition(0, 0, SLIP_WALL, 1.0)
{
    is_wall_flag = true;
}

SlipWallBC & SlipWallBC::operator=(const SlipWallBC &bc)
{
    BoundaryCondition::operator=(bc);
    return *this;
}

SlipWallBC::~SlipWallBC() {}

