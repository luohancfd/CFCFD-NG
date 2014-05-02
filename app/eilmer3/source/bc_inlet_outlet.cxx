// bc_inletoutlet.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_inlet_outlet.hh"
#include "kernel.hh"
#include "math.h"

//------------------------------------------------------------------------

InletOutletBC::InletOutletBC(Block *bdp, int which_boundary, double Pout_, 
                         double I_turb_, double u_turb_lam_, 
                         double Tout_, bool use_Tout_, int x_order_)
    : BoundaryCondition(bdp, which_boundary, INLET_OUTLET), 
      Pout(Pout_), I_turb(I_turb_), u_turb_lam(u_turb_lam_), 
      Tout(Tout_), use_Tout(use_Tout_), x_order(x_order_) 
{}

InletOutletBC::InletOutletBC(const InletOutletBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      Pout(bc.Pout), I_turb(bc.I_turb), u_turb_lam(bc.u_turb_lam), 
      Tout(bc.Tout), use_Tout(bc.use_Tout), x_order(bc.x_order) 
{}

InletOutletBC::InletOutletBC()
    : BoundaryCondition(0, 0, INLET_OUTLET),
      Pout(0.0), I_turb(0.0), u_turb_lam(0.0), 
      Tout(0.0), use_Tout(false), x_order(0) 
{}

InletOutletBC & InletOutletBC::operator=(const InletOutletBC &bc)
{
    BoundaryCondition::operator=(bc);
    Pout = bc.Pout; // Ok for self-assignment.
    I_turb = bc.I_turb;
    u_turb_lam = bc.u_turb_lam;
    Tout = bc.Tout;
    use_Tout = bc.use_Tout;
    x_order = bc.x_order;
    return *this;
}

InletOutletBC::~InletOutletBC() {}

void InletOutletBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "Pout= " << Pout << endl;
    cout << lead_in << "I_turb= " << I_turb << endl;
    cout << lead_in << "u_turb_lam= " << u_turb_lam << endl;
    cout << lead_in << "Tout= " << Tout << endl;
    cout << lead_in << "use_Tout=" << use_Tout << endl;
    cout << lead_in << "x_order= " << x_order << endl;
    return;
}

int InletOutletBC::apply_convective(double t)
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
} // end InletOutletBC::apply_convective(double t)


int InletOutletBC::apply_viscous(double t)
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    Block & bd = *bdp;
    double U;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                ConservedQuantities &F = *(IFace->F);
                U = sqrt( pow(fs.vel.x,2)+pow(fs.vel.y,2)+pow(fs.vel.z,2) );
                if (F.mass < 0) {
                     fs.tke = 1.5*pow((I_turb*U),2) ;
                     fs.omega = fs.gas->rho*fs.tke/fs.gas->mu/u_turb_lam ;           
                } // end if
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                ConservedQuantities &F = *(IFace->F);
                U = sqrt( pow(fs.vel.x,2)+pow(fs.vel.y,2)+pow(fs.vel.z,2) );
                if (F.mass < 0) {
                     fs.tke = 1.5*pow((I_turb*U),2) ;
                     fs.omega = fs.gas->rho*fs.tke/fs.gas->mu/u_turb_lam ;           
                } // end if
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                ConservedQuantities &F = *(IFace->F);
                U = sqrt( pow(fs.vel.x,2)+pow(fs.vel.y,2)+pow(fs.vel.z,2) );
                if (F.mass > 0) {
                     fs.tke = 1.5*pow((I_turb*U),2) ;
                     fs.omega = fs.gas->rho*fs.tke/fs.gas->mu/u_turb_lam ;           
                } // end if
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                ConservedQuantities &F = *(IFace->F);
                U = sqrt( pow(fs.vel.x,2)+pow(fs.vel.y,2)+pow(fs.vel.z,2) );
                if (F.mass > 0) {
                     fs.tke = 1.5*pow((I_turb*U),2) ;
                     fs.omega = fs.gas->rho*fs.tke/fs.gas->mu/u_turb_lam ;           
                } // end if
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                ConservedQuantities &F = *(IFace->F);
                U = sqrt( pow(fs.vel.x,2)+pow(fs.vel.y,2)+pow(fs.vel.z,2) );
                if (F.mass < 0) {
                     fs.tke = 1.5*pow((I_turb*U),2) ;
                     fs.omega = fs.gas->rho*fs.tke/fs.gas->mu/u_turb_lam ;           
                } // end if
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                ConservedQuantities &F = *(IFace->F);
                U = sqrt( pow(fs.vel.x,2)+pow(fs.vel.y,2)+pow(fs.vel.z,2) );
                if (F.mass > 0) {
                     fs.tke = 1.5*pow((I_turb*U),2) ;
                     fs.omega = fs.gas->rho*fs.tke/fs.gas->mu/u_turb_lam ;           
                } // end if
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_viscous not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
} // end InletOutletBC::apply_viscous()


