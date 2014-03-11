// bc_adjacent.cxx
#include <math.h> 
#include <vector>
#include <algorithm>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "kernel.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_adjacent.hh"

//------------------------------------------------------------------------

AdjacentBC::AdjacentBC(Block *bdp, int which_boundary, 
		       int other_block, int other_face, int _neighbour_orientation,
		       bool _reorient_vector_quantities, vector<double>& _Rmatrix)
    : BoundaryCondition(bdp, which_boundary, ADJACENT) 
{
    neighbour_block = other_block;
    neighbour_face = other_face;
    neighbour_orientation = _neighbour_orientation;
    reorient_vector_quantities = _reorient_vector_quantities;
    Rmatrix = _Rmatrix;
}

AdjacentBC::AdjacentBC(const AdjacentBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code) 
{
    neighbour_block = bc.neighbour_block;
    neighbour_face = bc.neighbour_face;
    neighbour_orientation = bc.neighbour_orientation;
    reorient_vector_quantities = bc.reorient_vector_quantities;
    Rmatrix = bc.Rmatrix;
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
    cout << lead_in << "reorient_vector_quantities= " << reorient_vector_quantities << endl;
    cout << lead_in << "Rmatrix= " << Rmatrix << endl;
    return;
}

int AdjacentBC::apply_convective(double t)
{
    // At this point, we assume that the exchange functions have
    // copied the other-block data into the appropriate ghost cells.
    
    // RJG note (25-Jul-2013): This method must remain in place
    // despite the fact that is appears to be "doing nothing".
    // What it actually does is PREVENT the default application
    // of a slip wall at the boundary interface.
    // If this method is removed, the ghost cells
    // are filled with reflections of the interior cells
    // and the reconstruction/flux computation gives the effect of a slip wall.
    // This is because boundary data is exchanged BEFORE the
    // boundary conditions are applied.

    if ( !reorient_vector_quantities ) return SUCCESS;

    // 28-Nov-2013 Put the reorientation code in here.
    size_t i, j, k;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		reorient_vector_quantities_in_cell(bd.get_cell(i,j+1,k), Rmatrix); // ghost cell 1.
		reorient_vector_quantities_in_cell(bd.get_cell(i,j+2,k), Rmatrix); // ghost cell 2.
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		reorient_vector_quantities_in_cell(bd.get_cell(i+1,j,k), Rmatrix); // ghost cell 1.
		reorient_vector_quantities_in_cell(bd.get_cell(i+2,j,k), Rmatrix); // ghost cell 2.
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		reorient_vector_quantities_in_cell(bd.get_cell(i,j-1,k), Rmatrix); // ghost cell 1.
		reorient_vector_quantities_in_cell(bd.get_cell(i,j-2,k), Rmatrix); // ghost cell 2.
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		reorient_vector_quantities_in_cell(bd.get_cell(i-1,j,k), Rmatrix); // ghost cell 1.
		reorient_vector_quantities_in_cell(bd.get_cell(i-2,j,k), Rmatrix); // ghost cell 2.
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		reorient_vector_quantities_in_cell(bd.get_cell(i,j,k+1), Rmatrix); // ghost cell 1.
		reorient_vector_quantities_in_cell(bd.get_cell(i,j,k+2), Rmatrix); // ghost cell 2.
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		reorient_vector_quantities_in_cell(bd.get_cell(i,j,k-1), Rmatrix); // ghost cell 1.
		reorient_vector_quantities_in_cell(bd.get_cell(i,j,k-2), Rmatrix); // ghost cell 2.
	    } // end j loop
	} // for i
    } // end switch...

    return SUCCESS;
}

// Helper functions
void apply_matrix_transform(const std::vector<double>& Rmatrix, 
			    const std::vector<double>& oldv,
			    std::vector<double> &newv)
{
    // Write out the matrix multiplication, long-hand.
    newv[0] = Rmatrix[0]*oldv[0] + Rmatrix[1]*oldv[1] + Rmatrix[2]*oldv[2];
    newv[1] = Rmatrix[3]*oldv[0] + Rmatrix[4]*oldv[1] + Rmatrix[5]*oldv[2];
    newv[2] = Rmatrix[6]*oldv[0] + Rmatrix[7]*oldv[1] + Rmatrix[8]*oldv[2];
}

void reorient_vector_quantities_in_cell(FV_Cell *c, const std::vector<double>& Rmatrix) 
{
    global_data &G = *get_global_data_ptr();
    std::vector<double>oldv(3);
    std::vector<double>newv(3);
    oldv[0] = c->fs->vel.x; oldv[1] = c->fs->vel.y; oldv[2] = c->fs->vel.z;
    apply_matrix_transform(Rmatrix, oldv, newv);
    c->fs->vel.x = newv[0]; c->fs->vel.y = newv[1]; c->fs->vel.z = newv[2];
    if ( G.MHD ) {
	oldv[0] = c->fs->B.x; oldv[1] = c->fs->B.y; oldv[2] = c->fs->B.z;
	apply_matrix_transform(Rmatrix, oldv, newv);
	c->fs->B.x = newv[0]; c->fs->B.y = newv[1]; c->fs->B.z = newv[2];
    }
}
