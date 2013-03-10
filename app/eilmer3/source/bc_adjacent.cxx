// bc_adjacent.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_adjacent.hh"

//------------------------------------------------------------------------

AdjacentBC::AdjacentBC( Block *bdp, int which_boundary, 
			int other_block, int other_face,
			int neighbour_orientation)
    : BoundaryCondition(bdp, which_boundary, ADJACENT, "AdjacentBC", 0, false, false, 
			other_block, other_face, neighbour_orientation) 
{}

AdjacentBC::AdjacentBC( const AdjacentBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation) 
{}

AdjacentBC::AdjacentBC()
    : BoundaryCondition(0, 0, ADJACENT, "AdjacentBC", 0, false, false, -1, -1, 0)
{}

AdjacentBC & AdjacentBC::operator=(const AdjacentBC &bc)
{
    BoundaryCondition::operator=(bc);
    return *this;
}

AdjacentBC::~AdjacentBC() {}

int AdjacentBC::apply_inviscid( double t )
{
    // Do nothing here; the real work 
    // is delegated to the exchange functions.
    return SUCCESS;
}
