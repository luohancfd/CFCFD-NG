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

FixedPOutBC::FixedPOutBC( Block *bdp, int which_boundary, double Pout, int x_order )
    : BoundaryCondition(bdp, which_boundary, FIXED_P_OUT, "FixedPOutBC",
			x_order, false, false, -1, -1, 0), 
      Pout(Pout) 
{}

FixedPOutBC::FixedPOutBC( const FixedPOutBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation), 
      Pout(bc.Pout) 
{}

FixedPOutBC::FixedPOutBC()
    : BoundaryCondition(0, 0, FIXED_P_OUT, "FixedPOutBC",
			0, false, false, -1, -1, 0), 
      Pout(0.0) 
{}

FixedPOutBC & FixedPOutBC::operator=(const FixedPOutBC &bc)
{
    BoundaryCondition::operator=(bc);
    Pout = bc.Pout; // Ok for self-assignment.
    return *this;
}

FixedPOutBC::~FixedPOutBC() {}

int FixedPOutBC::apply_inviscid( double t )
{
    // Fill ghost cells with data from just inside the boundary
    // using zero-order extrapolation (i.e. just copy the data)
    // and then impose a specified pressure.
    //
    // We assume that this boundary is an outflow boundary.
    size_t i, j, k;
    FV_Cell *src_cell, *dest_cell;
    Gas_model *gmodel = get_gas_model_ptr();
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		dest_cell->fs->gas->p = Pout;
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j,k-2);
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
