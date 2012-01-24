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

SupersonicInBC::SupersonicInBC( Block &bdp, int which_boundary, 
				int inflow_condition_id )
    : BoundaryCondition(bdp, which_boundary, SUP_IN, "SupersonicIn", 
			false, false, -1, -1, 0), 
      inflow_condition_id(inflow_condition_id) {}

SupersonicInBC::SupersonicInBC( const SupersonicInBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.is_wall_flag, bc.use_udf_flux_flag,
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
    global_data &gdp = *get_global_data_ptr();
    CFlowCondition *gsp = gdp.gas_state[inflow_condition_id];

    switch ( which_boundary ) {
    case NORTH:
	j = bdp.jmax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		dest_cell = bdp.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifj(i,j+1,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bdp.imax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifi(i+1,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bdp.jmin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		dest_cell = bdp.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifj(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*gsp);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bdp.imin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifi(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bdp.kmax;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifk(i,j,k+1);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*gsp);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bdp.kmin;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		dest_cell = bdp.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*gsp);
		// Although this is principally an inviscid BC,
		// we need the face values for derivatives.
		dest_face = bdp.get_ifk(i,j,k);
		dest_face->fs->copy_values_from(*gsp);
		dest_cell = bdp.get_cell(i,j,k-2);
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

