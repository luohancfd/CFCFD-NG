// bc_adjacent.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_adjacent.hh"

//------------------------------------------------------------------------

AdjacentBC::AdjacentBC(Block *bdp, int which_boundary, 
		       int other_block, int other_face,
		       int _neighbour_orientation)
    : BoundaryCondition(bdp, which_boundary, ADJACENT) 
{
    neighbour_block = other_block;
    neighbour_face = other_face;
    neighbour_orientation = _neighbour_orientation;
}

AdjacentBC::AdjacentBC(const AdjacentBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code) 
{
    neighbour_block = bc.neighbour_block;
    neighbour_face = bc.neighbour_face;
    neighbour_orientation = bc.neighbour_orientation;
}

AdjacentBC::AdjacentBC()
    : BoundaryCondition(0, 0, ADJACENT)
{}

AdjacentBC & AdjacentBC::operator=(const AdjacentBC &bc)
{
    BoundaryCondition::operator=(bc);
    return *this;
}

AdjacentBC::~AdjacentBC() {}

void AdjacentBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "neighbour_block= " << neighbour_block << endl;
    cout << lead_in << "neighbour_face= " << neighbour_face 
	 << " (" << get_face_name(neighbour_face) << ")" << endl;
    cout << lead_in << "neighbour_orientation= " << neighbour_orientation << endl;
    return;
}

int AdjacentBC::apply_convective(double t)
{
    // Do nothing here; the real work 
    // is delegated to the exchange functions.
    
    // RJG note (25-Jul-2013): This method must remain in place
    // despite the fact that is appears to be "doing nothing".
    // What it actually does is PREVENT the default application
    // of a slip wall at the boundary interface.
    // If this method is removed, the ghost cells
    // are filled with reflections of the interior cells
    // and the reconstruction/flux computation gives the effect of a slip wall.
    // This is because boundary data is exchanged BEFORE the
    // boundary conditions are applied.
    return SUCCESS;
}
