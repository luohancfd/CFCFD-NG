// bc_supersonic_in.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_supersonic_in.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

SupersonicInBC::SupersonicInBC( Block *bdp, int which_boundary, 
				int inflow_condition_id )
    : BoundaryCondition(bdp, which_boundary, SUP_IN, "SupersonicIn", 
			0, false, false, -1, -1, 0), 
      inflow_condition_id(inflow_condition_id) {}

SupersonicInBC::SupersonicInBC( const SupersonicInBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation),
      inflow_condition_id(bc.inflow_condition_id) {}

SupersonicInBC::~SupersonicInBC() {}

int SupersonicInBC::apply_inviscid( double t )
{
    // Set up ghost cells with inflow state. 
    int i, j, k;
    FV_Cell *dest_cell;
    FV_Interface *dest_face;
    global_data &gd = *get_global_data_ptr();
    CFlowCondition *gsp = gd.gas_state[inflow_condition_id];
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bd.get_ifj(i,j+1,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bd.get_ifi(i+1,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bd.get_ifj(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bd.get_ifi(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bd.get_ifk(i,j,k+1);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bd.get_ifk(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    } // end switch

    return SUCCESS;
}

