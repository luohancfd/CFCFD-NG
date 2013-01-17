// bc_shock_fitting_in.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_shock_fitting_in.hh"
#include "kernel.hh"
#include "one_d_interp_scalar.hh"

//------------------------------------------------------------------------

ShockFittingInBC::ShockFittingInBC( Block &bdp, int which_boundary, 
				    int inflow_condition_id )
    : BoundaryCondition(bdp, which_boundary, SHOCK_FITTING_IN, "ShockFittingIn", 
			0, false, false, -1, -1, 0), 
      inflow_condition_id(inflow_condition_id) {}

ShockFittingInBC::ShockFittingInBC( const ShockFittingInBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation),
      inflow_condition_id(bc.inflow_condition_id) {}

ShockFittingInBC::~ShockFittingInBC() {}

int ShockFittingInBC::apply_inviscid( double t )
{
    // Set up ghost cells with inflow state. 
    int i, j, k;
    FV_Cell *cL1;
    FV_Cell *cL0;
    FV_Cell *cR0;
    FV_Cell *cR1;
    FV_Cell *cR2;
    FV_Interface *IFaceL;
    FV_Interface *IFaceR;
    FV_Cell *dest_cell;
    FV_Interface *dest_face;
    FV_Vertex *vtx;
    FV_Vertex *wvtx;
    Vector3 trv;
    global_data &gdp = *get_global_data_ptr();
    CFlowCondition *gsp = gdp.gas_state[inflow_condition_id];

    switch ( which_boundary ) {
    case NORTH:
	printf( "Error: ShockFittingInBC not implemented for boundary %d\n, please use West boundary.", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
	break;
    case EAST:
	printf( "Error: ShockFittingInBC not implemented for boundary %d\n, please use West boundary.", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
	break;
    case SOUTH:
	printf( "Error: ShockFittingInBC not implemented for boundary %d\n, please use West boundary.", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
	break;
    case WEST:
	i = bdp.imin;
	//dest_face = bdp.get_ifi(i,bdp.jmin-1,bdp.kmin);
	//dest_face->fs->copy_values_from(*gsp);
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		cL1 = bdp.get_cell(i-2,j,k);
		cL1->copy_values_from(*gsp);
		cL0 = bdp.get_cell(i-1,j,k);
		cL0->copy_values_from(*gsp);
		cR0 = bdp.get_cell(i,j,k);
		cR1 = bdp.get_cell(i+1,j,k);
		cR2 = bdp.get_cell(i+2,j,k);
		IFaceL = bdp.get_ifi(i-1,j,k);
		IFaceL->fs->copy_values_from(*gsp);
		IFaceR = bdp.get_ifi(i,j,k);
		IFaceR->fs->copy_values_from(*gsp);
		calculate_shock_speed(cL1, cL0, cR0, cR1, cR2, 
				      cL1->iLength, cL0->iLength, cR0->iLength, cR1->iLength, cR2->iLength, 
				      IFaceL, IFaceR);
		compute_boundary_flux(IFaceL, IFaceR, bdp.omegaz);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	printf( "Error: ShockFittingInBC not implemented for boundary %d\n, please use West boundary.", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
	break;
    case BOTTOM:
	printf( "Error: ShockFittingInBC not implemented for boundary %d\n, please use West boundary.", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	exit(NOT_IMPLEMENTED_ERROR);
    } // end switch

    return SUCCESS;
}

/// \brief Calculate shock speed at interface. 
///
int ShockFittingInBC::calculate_shock_speed(FV_Cell *cL1, FV_Cell *cL0, FV_Cell *cR0, FV_Cell *cR1, FV_Cell *cR2, 
					    double lenL1, double lenL0, double lenR0, double lenR1, double lenR2, 
					    const FV_Interface *IFaceL, FV_Interface *IFaceR)
{
#   define EPSILON 1.0e-12
#   define ALPHA 0.5
    
    global_data &gdp = *get_global_data_ptr();
    double time_weight = gdp.shock_fitting_speed_factor;
    if ( get_shock_fitting_decay_flag() ) {
        time_weight = time_weight - time_weight*(exp(-gdp.sim_time/gdp.max_time) - 1.0) / 
                             (exp(-1.0) - 1.0);
    }
    // Downwind pressure recontruction as per Ian Johnston's thesis.
    //double delRminus = 2.0 * (cR1->fs->gas->p - cR0->fs->gas->p) / (lenR1 + lenR0);
    //double delRplus = 2.0 * (cR2->fs->gas->p - cR1->fs->gas->p) / (lenR2 + lenR1);
    //double delLminus = 2.0 * (cR0->fs->gas->p - cL0->fs->gas->p) / (lenR0 + lenL0);
    //double sR = (delRminus * delRplus + fabs(delRminus * delRplus)) /
    //    (delRminus * delRminus + delRplus * delRplus + EPSILON);
    //double delR = 0.5 * sR * ( delRplus + delRminus );
    Gas_data *gL = IFaceL->fs->gas;
    Gas_data *gR = IFaceR->fs->gas;
    
    if ( gdp.sim_time >= gdp.t_shock ) {
        //if ( delLminus > 1.1 * delR ) {
        if ( fabs(cR0->fs->gas->rho - cL0->fs->gas->rho) / 
             max(cL0->fs->gas->rho, cR0->fs->gas->rho) > 0.2 ) {
	    // Stronger pressure gradient on left suggests shock, use downwind reconstruction
	    // to find shock speed.
	    onesided_interp(*cR0, *cR1, *cR2,
			    lenR0, lenR1, lenR2,
			    *IFaceR->fs);
            double ws1 = (gL->rho * dot(IFaceL->fs->vel, IFaceR->n) - 
                          gR->rho * dot(IFaceR->fs->vel, IFaceR->n)) / ( gL->rho - gR->rho);
            double ws2 = dot(IFaceL->fs->vel, IFaceR->n) - (1.0 / gL->rho) *
		sqrt( fabs( (gL->p - gR->p) / (1.0/gL->rho - 1.0/gR->rho) ) );
	    double ws = (ALPHA*ws1 + (1-ALPHA)*ws2);
	    // Set interface velocity to the lower of the calculated shock speed and
	    // the lower of the two flow velocities, for stability
	    double flowv = min(vabs(IFaceL->fs->vel), vabs(IFaceR->fs->vel));
            IFaceR->vel.x = time_weight * copysign(min( fabs(ws), flowv ), ws);
            IFaceR->vel.y = 0.0;
            IFaceR->vel.z = 0.0;
            IFaceR->vel.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);
        } else {
	    // Probably no shock so move boundary at flow speed.
	    //mach_weighted_one_d_interp(*cL1, *cL0, *cR0, *cR1,
	    //			       lenL1, lenL0, lenR0, lenR1,
	    //			       *IFaceL->fs, *IFaceR->fs);
            IFaceR->vel.x = time_weight * min(vabs(IFaceL->fs->vel), vabs(IFaceR->fs->vel));
            IFaceR->vel.y = 0.0;
            IFaceR->vel.z = 0.0;
            IFaceR->vel.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);
        }
    } else {
	IFaceR->vel.x = 0.0;
        IFaceR->vel.y = 0.0;
        IFaceR->vel.z = 0.0;
        IFaceR->vel.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);
    }
    return SUCCESS;
} // end ShockFittingInBC::calculate_post_shock_properties()

/// \brief Calculate inviscid fluxes based only on freestream conditions. 
///

int ShockFittingInBC::compute_boundary_flux(FV_Interface *IFaceL, FV_Interface *IFaceR, double omegaz)
{
    ConservedQuantities &F = *(IFaceR->F);
    FlowState *IFace_flow_state = IFaceL->fs;

    IFaceR->vel.transform_to_local(IFaceR->n, IFaceR->t1, IFaceR->t2);
    IFace_flow_state->vel.transform_to_local(IFaceR->n, IFaceR->t1, IFaceR->t2);
    // also transform the magnetic field
    if (get_mhd_flag() == 1) {
	IFace_flow_state->B.transform_to_local(IFaceR->n, IFaceR->t1, IFaceR->t2);
    }
    set_interface_flux(*IFaceR, IFace_flow_state);
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	double x = IFaceR->pos.x;
	double y = IFaceR->pos.y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy,
	// so we need to take -(u**2)/2 off the total energy Shock.
	// Note that rotating frame velocity u = omegaz * r.
	F.total_energy -= F.mass * 0.5*omegaz*omegaz*rsq;
    }
    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);
    // also transform the interface velocities
    IFaceR->vel.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);	  
    // also transform the magnetic field
    if (get_mhd_flag() == 1) {
	F.B.transform_to_global(IFaceR->n, IFaceR->t1, IFaceR->t2);
    }
    
    return SUCCESS;
} // end ShockFittingInBC::inflow_fluxes()
