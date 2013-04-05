// bc_shock_fitting_in.cxx

#include <stdexcept>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "cell.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_shock_fitting_in.hh"
#include "kernel.hh"
#include "one_d_interp_scalar.hh"
#include "one_d_interp.hh"

//------------------------------------------------------------------------

ShockFittingInBC::
ShockFittingInBC( Block *bdp, int which_boundary, int inflow_condition_id )
    : BoundaryCondition(bdp, which_boundary, SHOCK_FITTING_IN, "ShockFittingIn", 
			0, false, false, -1, -1, 0), 
      inflow_condition_id(inflow_condition_id) 
{}

ShockFittingInBC::
ShockFittingInBC( const ShockFittingInBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation),
      inflow_condition_id(bc.inflow_condition_id) 
{}

ShockFittingInBC::
ShockFittingInBC()
    : BoundaryCondition(0, 0, SHOCK_FITTING_IN, "ShockFittingIn", 
			0, false, false, -1, -1, 0), 
      inflow_condition_id(0) 
{}

ShockFittingInBC & 
ShockFittingInBC::operator=(const ShockFittingInBC &bc)
{
    BoundaryCondition::operator=(bc);
    inflow_condition_id = bc.inflow_condition_id; // Benign for self-assignment.
    return *this;
}

ShockFittingInBC::~ShockFittingInBC() {}

int ShockFittingInBC::apply_inviscid( double t )
// Copies from FlowCondition to ghost cells.
{
    // Set up ghost cells with inflow state. 
    size_t i, j, k;
    FV_Cell *cL1, *cL0, *cR0, *cR1, *cR2;
    FV_Interface *IFaceR;
    global_data &gd = *get_global_data_ptr();
    CFlowCondition *gsp = gd.gas_state[inflow_condition_id];
    Block & bd = *bdp;

    switch ( which_boundary ) {
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
		IFaceR = bd.get_ifi(i,j,k);
		IFaceR->fs->copy_values_from(*gsp);
		calculate_shock_speed(*cL0, *cR0, *cR1, *cR2, 
				      cL0->iLength, cR0->iLength, cR1->iLength, cR2->iLength, 
				      *IFaceR);
		// Necessary to provide a symmetry boundary condition when appying the spatial filter.
		cL1->copy_values_from(*cR1, COPY_FLOW_STATE, 0);
		cL0->copy_values_from(*cR0, COPY_FLOW_STATE, 0);
	    } // end j loop
	} // for k
 	break;
    default:
	cout << "ShockFittingInBC not implemented for " 
	     << get_face_name(which_boundary) << " boundary." << endl;
        cout << "    Please use West boundary." << endl;
	throw runtime_error("ShockFittingInBC not implemented for this boundary.");
    } // end switch

    return SUCCESS;
}

int ShockFittingInBC::apply_viscous( double t )
// Copies interior-cell flow-properties to interface.
{
    size_t i, j, k;
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
// FIX-ME moving-grid: Andrew, the function name seems misleading to me.
// It's really calculating the desired grid velocity, (which should 
// normally be the shock speed) and that will eventually be zero, yes?
// Finally, should it retain the tangential velocity components at the interface?
// My guess is no.
///
int ShockFittingInBC::calculate_shock_speed(const FV_Cell &cL0, const FV_Cell &cR0,
					    const FV_Cell &cR1, const FV_Cell &cR2, 
					    double lenL0, double lenR0, double lenR1, double lenR2, 
					    FV_Interface &IFaceR)
{
#   define EPSILON 1.0e-12
#   define ALPHA 0.5
    global_data &gd = *get_global_data_ptr();
    CFlowCondition *gsp = gd.gas_state[inflow_condition_id];
    Gas_model *gmodel = get_gas_model_ptr();
    FlowState fsL(gmodel);
    fsL.copy_values_from(*gsp);
    Gas_data &gL = *(fsL.gas);
    Gas_data &gR = *(IFaceR.fs->gas);
    double time_weight = gd.shock_fitting_speed_factor;
    if ( get_shock_fitting_decay_flag() ) {
        time_weight = time_weight - time_weight*(exp(-gd.sim_time/gd.max_time) - 1.0) / 
                             (exp(-1.0) - 1.0);
    }
    
    if ( gd.sim_time >= gd.t_shock ) {
	// Detect shock from density jump.
	double rel_density_jump = fabs(cR0.fs->gas->rho - cL0.fs->gas->rho) / 
	    max(cL0.fs->gas->rho, cR0.fs->gas->rho);
        if ( rel_density_jump > 0.2 ) {
	    onesided_interp(cR0, cR1, cR2, lenR0, lenR1, lenR2, *(IFaceR.fs));
            double ws1 = (gL.rho * dot(fsL.vel, IFaceR.n) - 
                          gR.rho * dot(IFaceR.fs->vel, IFaceR.n)) / (gL.rho - gR.rho);
            double ws2 = dot(fsL.vel, IFaceR.n) - copysign(1.0, gR.p - gL.p) / gL.rho *
		         sqrt( fabs( (gR.p - gL.p) / (1.0/gL.rho - 1.0/gR.rho) ) );
	    double ws = (ALPHA*ws1 + (1-ALPHA)*ws2);
	    // Set interface normal velocity to the lower of the calculated shock speed and
	    // the lower of the two flow velocities either side of the interface, for stability.
	    double flowv = min(vabs(fsL.vel), vabs(IFaceR.fs->vel));
            IFaceR.vel.x = time_weight * copysign(min( fabs(ws), flowv ), ws);
            IFaceR.vel.y = 0.0;
            IFaceR.vel.z = 0.0;
            IFaceR.vel.transform_to_global(IFaceR.n, IFaceR.t1, IFaceR.t2);
        } else {
	    // Probably no shock so move boundary at the lower of the two flow velocities
	    // either side of the interface.
            IFaceR.vel.x = time_weight * min(vabs(fsL.vel), vabs(IFaceR.fs->vel));
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
} // end ShockFittingInBC::calculate_shock_speed()
