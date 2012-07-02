// bc_fixed_p_out.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_fixed_p_out.hh"
#include "kernel.hh"

//------------------------------------------------------------------------

FixedPOutBC::FixedPOutBC( Block &bdp, int which_boundary, double Pout, int x_order )
    : BoundaryCondition(bdp, which_boundary, FIXED_P_OUT, "FixedPOutBC",
			x_order, false, false, -1, -1, 0), 
      Pout(Pout) 
{}

FixedPOutBC::FixedPOutBC( const FixedPOutBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation), 
      Pout(Pout) 
{}

FixedPOutBC::~FixedPOutBC() {}

int FixedPOutBC::apply_inviscid( double t )
{
    // Fill ghost cells with data from just inside the boundary
    // using zero-order extrapolation (i.e. just copy the data)
    // and then impose a specified pressure.
    //
    // We assume that this boundary is an outflow boundary.
    int i, j, k;
    FV_Cell *src_cell, *dest_cell;
    Gas_model *gmodel = get_gas_model_ptr();

    switch ( which_boundary ) {
    case NORTH:
	j = bdp.jmax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		src_cell = bdp.get_cell(i,j,k);
		dest_cell = bdp.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bdp.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bdp.imax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		dest_cell = bdp.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bdp.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bdp.jmin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		src_cell = bdp.get_cell(i,j,k);
		dest_cell = bdp.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bdp.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bdp.imin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		dest_cell = bdp.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bdp.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bdp.kmax;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		dest_cell = bdp.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bdp.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bdp.kmin;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		dest_cell = bdp.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bdp.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
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
