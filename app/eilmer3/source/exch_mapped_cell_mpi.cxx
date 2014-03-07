// exch_mapped_cell_mpi.cxx
// Copy mapped-cell data in a distributed-memory context.
// 
// PJ, 08-Mar-2014

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cell.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"

int copy_mapped_cell_data_via_mpi(int type_of_copy, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    // [TODO] -- PJ 2014-03-08 -- mapped-cell BCs are a work in progress.
    // Need to do something sensible to copy the flow data.
    bool found_mapped_cell_bc = false;
    for ( Block *bdp : G.my_blocks ) {
	int number_faces = (G.dimensions == 3 ? 6: 4);
	for ( int iface = 0; iface < number_faces; ++iface ) {
	    if ( bdp->bcp[iface]->type_code == MAPPED_CELL ) found_mapped_cell_bc = true;
	}
    }
    if ( found_mapped_cell_bc ) {
	throw std::runtime_error("Have not yet implemented mapped-cell exchange for MPI.");
    }
    return SUCCESS;
}
