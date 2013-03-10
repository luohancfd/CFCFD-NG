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

ShockFittingInBC::ShockFittingInBC( Block *bdp, int which_boundary, 
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
    FV_Cell *cL1, *cL0, *cR0, *cR1, *cR2;
    FV_Interface *IFaceL, *IFaceR;
    global_data &gdp = *get_global_data_ptr();
    CFlowCondition *gsp = gdp.gas_state[inflow_condition_id];
    Block & bd = *bdp;

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
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cL1 = bd.get_cell(i-2,j,k);
		cL1->copy_values_from(*gsp);
		cL0 = bd.get_cell(i-1,j,k);
		cL0->copy_values_from(*gsp);
		cR0 = bd.get_cell(i,j,k);
		cR1 = bd.get_cell(i+1,j,k);
		cR2 = bd.get_cell(i+2,j,k);
		IFaceL = bd.get_ifi(i-1,j,k);
		IFaceL->fs->copy_values_from(*gsp);
		IFaceR = bd.get_ifi(i,j,k);
		IFaceR->fs->copy_values_from(*gsp);
		calculate_shock_speed(*cL0, *cR0, *cR1, *cR2, 
				      cL0->iLength, cR0->iLength, cR1->iLength, cR2->iLength, 
				      *IFaceL, *IFaceR);
		// Necessary to stop edge effects when appying the spatial filter.
		cL1->copy_values_from(*cR1, COPY_FLOW_STATE);
		cL0->copy_values_from(*cR0, COPY_FLOW_STATE);
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

int ShockFittingInBC::apply_viscous( double t )
{
    int i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    Block & bd = *bdp;
    i = bd.imin;
    for (k = bd.kmin; k <= bd.kmax; ++k) {
	for (j = bd.jmin; j <= bd.jmax; ++j) {
	    cell = bd.get_cell(i,j,k);
	    IFace = cell->iface[WEST];
	    IFace->fs->copy_values_from(*(cell->fs));
	} // end j loop
    } // for k
    return SUCCESS;
} // end ShockFittingInBC::apply_viscous()


/// \brief Calculate shock speed at interface. 
///
int ShockFittingInBC::calculate_shock_speed(const FV_Cell &cL0, const FV_Cell &cR0,
					    const FV_Cell &cR1, const FV_Cell &cR2, 
					    double lenL0, double lenR0, double lenR1, double lenR2, 
					    const FV_Interface &IFaceL, FV_Interface &IFaceR)
{
#   define EPSILON 1.0e-12
#   define ALPHA 0.5
    global_data &gdp = *get_global_data_ptr();
    double time_weight = gdp.shock_fitting_speed_factor;
    if ( get_shock_fitting_decay_flag() ) {
        time_weight = time_weight - time_weight*(exp(-gdp.sim_time/gdp.max_time) - 1.0) / 
                             (exp(-1.0) - 1.0);
    }
    Gas_data &gL = *(IFaceL.fs->gas);
    Gas_data &gR = *(IFaceR.fs->gas);
    
    if ( gdp.sim_time >= gdp.t_shock ) {
	// Detect shock from density jump.
        if ( fabs(cR0.fs->gas->rho - cL0.fs->gas->rho) / 
             max(cL0.fs->gas->rho, cR0.fs->gas->rho) > 0.2 ) {
	    onesided_interp(cR0, cR1, cR2,
	    		    lenR0, lenR1, lenR2,
	    		    *(IFaceR.fs));
            double ws1 = (gL.rho * dot(IFaceL.fs->vel, IFaceR.n) - 
                          gR.rho * dot(IFaceR.fs->vel, IFaceR.n)) / (gL.rho - gR.rho);
            double ws2 = dot(IFaceL.fs->vel, IFaceR.n) - copysign(1.0, gR.p - gL.p) / gL.rho *
		         sqrt( fabs( (gR.p - gL.p) / (1.0/gL.rho - 1.0/gR.rho) ) );
	    double ws = (ALPHA*ws1 + (1-ALPHA)*ws2);
	    // Set interface normal velocity to the lower of the calculated shock speed and
	    // the lower of the two flow velocities either side of the interface, for stability.
	    double flowv = min(vabs(IFaceL.fs->vel), vabs(IFaceR.fs->vel));
            IFaceR.vel.x = time_weight * copysign(min( fabs(ws), flowv ), ws);
            IFaceR.vel.y = 0.0;
            IFaceR.vel.z = 0.0;
            IFaceR.vel.transform_to_global(IFaceR.n, IFaceR.t1, IFaceR.t2);
        } else {
	    // Probably no shock so move boundary at the lower of the two flow velocities
	    // either side of the interface.
            IFaceR.vel.x = time_weight * min(vabs(IFaceL.fs->vel), vabs(IFaceR.fs->vel));
            IFaceR.vel.y = 0.0;
            IFaceR.vel.z = 0.0;
            IFaceR.vel.transform_to_global(IFaceR.n, IFaceR.t1, IFaceR.t2);
        }
    } else {
	IFaceR.vel.x = 0.0;
        IFaceR.vel.y = 0.0;
        IFaceR.vel.z = 0.0;
        IFaceR.vel.transform_to_global(IFaceR.n, IFaceR.t1, IFaceR.t2);
    }
    return SUCCESS;
} // end ShockFittingInBC::calculate_post_shock_properties()
