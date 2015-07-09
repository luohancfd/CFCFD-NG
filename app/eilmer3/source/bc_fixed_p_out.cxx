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

FixedPOutBC::FixedPOutBC(Block *bdp, int which_boundary, double Pout_, 
			 double Tout_, bool use_Tout_, int x_order_)
    : BoundaryCondition(bdp, which_boundary, FIXED_P_OUT), 
      Pout(Pout_), Tout(Tout_), use_Tout(use_Tout_), x_order(x_order_) 
{}

FixedPOutBC::FixedPOutBC(const FixedPOutBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      Pout(bc.Pout), Tout(bc.Tout), use_Tout(bc.use_Tout), x_order(bc.x_order) 
{}

FixedPOutBC::FixedPOutBC()
    : BoundaryCondition(0, 0, FIXED_P_OUT),
      Pout(0.0), Tout(0.0), use_Tout(false), x_order(0) 
{}

FixedPOutBC & FixedPOutBC::operator=(const FixedPOutBC &bc)
{
    BoundaryCondition::operator=(bc);
    Pout = bc.Pout; // Ok for self-assignment.
    Tout = bc.Tout;
    use_Tout = bc.use_Tout;
    x_order = bc.x_order;
    return *this;
}

FixedPOutBC::~FixedPOutBC() {}

void FixedPOutBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "Pout= " << Pout << endl;
    cout << lead_in << "Tout= " << Tout << endl;
    cout << lead_in << "use_Tout=" << use_Tout << endl;
    cout << lead_in << "x_order= " << x_order << endl;
    return;
}

int FixedPOutBC::apply_convective(double t)
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
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();    

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= nmodes; ++imode )
			src_cell->iface[NORTH]->fs->gas->T[imode] = Tout;			
		}		
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;		
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= nmodes; ++imode )
			src_cell->iface[EAST]->fs->gas->T[imode] = Tout;			
		}		
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;		
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= nmodes; ++imode )
			src_cell->iface[SOUTH]->fs->gas->T[imode] = Tout;			
		}		
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;			
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= nmodes; ++imode )
			src_cell->iface[WEST]->fs->gas->T[imode] = Tout;			
		}		
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;		
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= nmodes; ++imode )
			src_cell->iface[TOP]->fs->gas->T[imode] = Tout;			
		}		
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;			
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= nmodes; ++imode )
			src_cell->iface[BOTTOM]->fs->gas->T[imode] = Tout;			
		}		
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;			
		}
		gmodel->eval_thermo_state_pT(*(dest_cell->fs->gas));  // make density consistent
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		dest_cell->fs->gas->p = Pout;
		if ( use_Tout ) {
		    for ( size_t imode=0; imode <= dest_cell->fs->gas->T.size(); ++imode )
			dest_cell->fs->gas->T[imode] = Tout;
		}
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
