// bc_mapped_cell.cxx
// PJ -- 06-march-2014 finally got serious about implementation

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
#include "bc_mapped_cell.hh"

//------------------------------------------------------------------------

MappedCellBC::MappedCellBC(Block *bdp, int which_boundary, 
			   std::string filename,
			   bool _reorient_vector_quantities, vector<double>& _Rmatrix)
    : BoundaryCondition(bdp, which_boundary, MAPPED_CELL) 
{
    global_data &G = *get_global_data_ptr();
    reorient_vector_quantities = _reorient_vector_quantities;
    Rmatrix = _Rmatrix;
    // Read mapping of ghost-cells to active-cells from a file.
    std::ifstream fstrm(filename);
    std::string text;
    int src_blk, src_i, src_j, src_k;
    std::vector<int> mapping;
    getline(fstrm, text); // Discard comment line.
    switch ( which_boundary ) {
    case NORTH:
    case SOUTH:
	// cout << "NORTH/SOUTH boundary" << endl;
	for ( size_t i = 0; i < bdp->nni; ++i ) {
	    for ( size_t k = 0; k < bdp->nnk; ++k ) {
		for ( size_t ghost_cell_count = 1; ghost_cell_count < 2; ++ghost_cell_count ) {
		    fstrm >> src_blk >> src_i >> src_j >> src_k;
		    // Translate the cell index to account for ghost cells around
		    // the periphery of the cell storage arrays.
		    src_i += G.nghost;
		    src_j += G.nghost;
		    if ( bdp->nnk > 1 ) src_k += G.nghost;
		    // cout << "ghost-cell=" << ghost_cell_count 
		    // 	 << " src_blk=" << src_blk << " src_i=" << src_i
		    // 	 << " src_j=" << src_j << " src_k=" << src_k << endl;
		    mapping.clear();
		    mapping.push_back(src_blk);
		    mapping.push_back(src_i);
		    mapping.push_back(src_j);
		    mapping.push_back(src_k);
		    mapped_cells.push_back(mapping);
		}
	    } // end for k
	} // end for i
	break;
    case EAST:
    case WEST:
	for ( size_t j = 0; j < bdp->nnj; ++j ) {
	    for ( size_t k = 0; k < bdp->nnk; ++k ) {
		for ( size_t ghost_cell_count = 1; ghost_cell_count < 2; ++ghost_cell_count ) {
		    fstrm >> src_blk >> src_i >> src_j >> src_k;
		    // Translate the cell index to account for ghost cells around
		    // the periphery of the cell storage arrays.
		    src_i += G.nghost;
		    src_j += G.nghost;
		    if ( bdp->nnk > 1 ) src_k += G.nghost;
		    mapping.clear();
		    mapping.push_back(src_blk);
		    mapping.push_back(src_i);
		    mapping.push_back(src_j);
		    mapping.push_back(src_k);
		    mapped_cells.push_back(mapping);
		}
	    } // end for k
	} // end for i
	break;
    case TOP:
    case BOTTOM:
	for ( size_t i = 0; i < bdp->nni; ++i ) {
	    for ( size_t j = 0; j < bdp->nnj; ++j ) {
		for ( size_t ghost_cell_count = 1; ghost_cell_count < 2; ++ghost_cell_count ) {
		    fstrm >> src_blk >> src_i >> src_j >> src_k;
		    // Translate the cell index to account for ghost cells around
		    // the periphery of the cell storage arrays.
		    src_i += G.nghost;
		    src_j += G.nghost;
		    if ( bdp->nnk > 1 ) src_k += G.nghost;
		    mapping.clear();
		    mapping.push_back(src_blk);
		    mapping.push_back(src_i);
		    mapping.push_back(src_j);
		    mapping.push_back(src_k);
		    mapped_cells.push_back(mapping);
		}
	    } // end for j
	} // end for i
    }
    fstrm.close();
}

MappedCellBC::MappedCellBC(const MappedCellBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code) 
{
    mapped_cells = bc.mapped_cells;
    incoming_mapped_cells = bc.incoming_mapped_cells;
    outgoing_mapped_cells = bc.outgoing_mapped_cells;
    reorient_vector_quantities = bc.reorient_vector_quantities;
    Rmatrix = bc.Rmatrix;
}

MappedCellBC::MappedCellBC()
    : BoundaryCondition(0, 0, MAPPED_CELL)
{}

MappedCellBC & MappedCellBC::operator=(const MappedCellBC &bc)
{
    BoundaryCondition::operator=(bc);
    return *this;
}

MappedCellBC::~MappedCellBC() {}

void MappedCellBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "mapped_cells.size()= " << mapped_cells.size() << endl;
    cout << lead_in << "reorient_vector_quantities= " << reorient_vector_quantities << endl;
    cout << lead_in << "Rmatrix= " << Rmatrix << endl;
    return;
}

int MappedCellBC::apply_convective(double t)
{
    // [TODO] -- copy data from source cells, probably at same level as exchange is done.

    if ( !reorient_vector_quantities ) return SUCCESS;

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
