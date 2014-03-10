// bc_adjacent.cxx
#include <math.h> 
#include <vector>
#include <algorithm>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "../../../lib/geometry2/source/geom.hh"
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
    double pi = 3.141592653589793238463;
    std::vector<double>matrix(9);
    std::vector<double>nB(3);
    std::vector<double>t1B(3);
    std::vector<double>t2B(3);
    std::vector<double>nA(3);
    std::vector<double>t1A(3);
    std::vector<double>t2A(3);
    double a11,a12,a13,a21,a22,a23,a31,a32,a33;
    double A11,A12,A13,A21,A22,A23,A31,A32,A33;
    double b11,b21,b31;
    double DET;

    nB[0] = -sin(Rmatrix[0]*pi/180.0); nB[1] = cos(Rmatrix[0]*pi/180.0); nB[2] = 0.0;
    t1B[0] = cos(Rmatrix[0]*pi/180.0); t1B[1] = sin(Rmatrix[0]*pi/180.0); t1B[2] = 0.0;
    t2B[0] = 0.0; t2B[1] = 0.0; t2B[2] = 1.0;
    nA[0] = -sin(Rmatrix[1]*pi/180.0); nA[1] = cos(Rmatrix[1]*pi/180.0); nA[2] = 0.0;
    t1A[0] = cos(Rmatrix[1]*pi/180.0); t1A[1] = sin(Rmatrix[1]*pi/180.0); t1A[2] = 0.0;
    t2A[0] = 0.0; t2A[1] = 0.0; t2A[2] = 1.0;

    a11 = nB[0]; a12 = nB[1]; a13 = nB[2];
    a21 = t1B[0]; a22 = t1B[1]; a23 = t1B[2];
    a31 = t2B[0]; a32 = t2B[1]; a33 = t2B[2];

    DET = a11*(a33*a22-a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12-a22*a13);
    A11 = (a33*a22-a32*a23) / DET;
    A12 = -(a33*a12-a32*a13) / DET;
    A13 = (a23*a12-a22*a13) /DET;
    A21 = -(a33*a21-a31*a23) / DET;
    A22 = (a33*a11-a31*a13) / DET;
    A23 = -(a23*a11-a21*a13) / DET;
    A31 = (a32*a21-a31*a22) /DET;
    A32 = -(a32*a11-a31*a12) / DET;
    A33 = (a22*a11-a21*a12) / DET;

    b11 = nA[0]; b21 = t1A[0]; b31 = t2A[0];
    matrix[0] = A11*b11 + A12*b21 + A13*b31;
    matrix[1] = A21*b11 + A22*b21 + A23*b31;
    matrix[2] = A31*b11 + A32*b21 + A33*b31;

    b11 = nA[1]; b21 = t1A[1]; b31 = t2A[1];
    matrix[3] = A11*b11 + A12*b21 + A13*b31;
    matrix[4] = A21*b11 + A22*b21 + A23*b31;
    matrix[5] = A31*b11 + A32*b21 + A33*b31;

    b11 = nA[2]; b21 = t1A[2]; b31 = t2A[2];
    matrix[6] = A11*b11 + A12*b21 + A13*b31;
    matrix[7] = A21*b11 + A22*b21 + A23*b31;
    matrix[8] = A31*b11 + A32*b21 + A33*b31;

    newv[0] = matrix[0]*oldv[0] + matrix[1]*oldv[1] + matrix[2]*oldv[2];
    newv[1] = matrix[3]*oldv[0] + matrix[4]*oldv[1] + matrix[5]*oldv[2];
    newv[2] = matrix[6]*oldv[0] + matrix[7]*oldv[1] + matrix[8]*oldv[2];

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
