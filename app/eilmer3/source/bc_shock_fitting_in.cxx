// bc_shock_fitting_in.cxx

#include <stdexcept>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "cell.hh"
#include "block.hh"
#include "bc.hh"
#include "bc_shock_fitting_in.hh"
#include "kernel.hh"
#include "one_d_interp.hh"

//------------------------------------------------------------------------

ShockFittingInBC::
ShockFittingInBC(Block *bdp, int which_boundary, int inflow_condition_id)
    : BoundaryCondition(bdp, which_boundary, SHOCK_FITTING_IN), 
      inflow_condition_id(inflow_condition_id) 
{
    ghost_cell_data_available = false;
}

ShockFittingInBC::
ShockFittingInBC(const ShockFittingInBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      inflow_condition_id(bc.inflow_condition_id) 
{
    ghost_cell_data_available = bc.ghost_cell_data_available;
}

ShockFittingInBC::
ShockFittingInBC()
    : BoundaryCondition(0, 0, SHOCK_FITTING_IN), 
      inflow_condition_id(0) 
{
    ghost_cell_data_available = false;
}

ShockFittingInBC & 
ShockFittingInBC::operator=(const ShockFittingInBC &bc)
{
    BoundaryCondition::operator=(bc);
    inflow_condition_id = bc.inflow_condition_id; // Benign for self-assignment.
    ghost_cell_data_available = bc.ghost_cell_data_available;
    return *this;
}

ShockFittingInBC::~ShockFittingInBC() {}

void ShockFittingInBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "inflow_condition_id= " << inflow_condition_id << endl;
    return;
}

int ShockFittingInBC::apply_convective(double t)
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
	throw std::runtime_error("ShockFittingInBC not implemented for this boundary.");
    } // end switch

    return SUCCESS;
}

int ShockFittingInBC::apply_viscous(double t)
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

//-----------------------------------------------------------------------------
// FIX-ME
// Andrew, I've moved the interpolation functions that only your shock-fitting uses.
// Does it really need such high-order reconstruction? 
// Can it do something simpler to get the shock speed? 

/// \brief One-sided one-dimensional reconstruction of a scalar quantity.
///
/// Implemented by Andrew Pastrello (I believe, PJ).  See Ian Johnston's thesis.
///
inline int onesided_interp_scalar( double qR0, double qR1, double qR2, 
				   double lenR0, double lenR1, double lenR2, 
				   double &qR)
{
    const double epsilon = 1.0e-12;

    double aR0, delRminus, delRplus, sR;
    
    // Set up differences and limiter values.
    aR0 = 0.5 * lenR0 / (lenR0 + 2.0*lenR1 + lenR2);
    delRminus = 2.0 * (qR1 - qR0) / (lenR1 + lenR0);
    delRplus = 2.0 * (qR2 - qR1) / (lenR2 + lenR1);
    if ( get_apply_limiter_flag() ) {
	// val Albada limiter as per Ian Johnston's thesis.
	sR = (delRminus * delRplus + fabs(delRminus * delRplus)) /
	    (delRminus * delRminus + delRplus * delRplus + epsilon);
    } else {
	sR = 1.0;
    }
    // The high-order reconstruction, possibly limited.
    qR = qR0 - 0.5 * aR0 * sR * ( (1 - sR) * delRplus * lenR1 + (1 + sR) * 
				  delRminus * (lenR1 + lenR2 + lenR0) );

    return SUCCESS;
} // end of onesided_interp_scalar()

/// \brief Reconstruct flow properties at an interface from FV_Cell properties from one side only.
///
/// This scheme uses the three cells to the right of the interface to do
/// a one-sided extrapolation of the flow properties. This is done to determine the shock
/// speed when shock-fitting.
int onesided_interp(const FV_Cell &cL0, const FV_Cell &cR0, const FV_Cell &cR1,
		    double cL0Length, double cR0Length, double cR1Length,
		    FlowState &Rght )
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    Gas_data &gL0 = *(cL0.fs->gas);
    Gas_data &gR0 = *(cR0.fs->gas);
    Gas_data &gR1 = *(cR1.fs->gas);
    int nsp = gmodel->get_number_of_species();
    
    // Low-order reconstruction just copies data from adjacent FV_Cell.
    Rght.copy_values_from(*(cR0.fs));
    if ( G.Xorder > 1 ) {
	// High-order reconstruction for some properties.
	onesided_interp_scalar(cL0.fs->vel.x, cR0.fs->vel.x, cR1.fs->vel.x,
			       cL0Length, cR0Length, cR1Length, Rght.vel.x);
	onesided_interp_scalar(cL0.fs->vel.y, cR0.fs->vel.y, cR1.fs->vel.y,
			       cL0Length, cR0Length, cR1Length, Rght.vel.y);
	onesided_interp_scalar(cL0.fs->vel.z, cR0.fs->vel.z, cR1.fs->vel.z,
			       cL0Length, cR0Length, cR1Length, Rght.vel.z);
	if ( G.MHD ) {
	    onesided_interp_scalar(cL0.fs->B.x, cR0.fs->B.x, cR1.fs->B.x,
				   cL0Length, cR0Length, cR1Length, Rght.B.x);
	    onesided_interp_scalar(cL0.fs->B.y, cR0.fs->B.y, cR1.fs->B.y,
				   cL0Length, cR0Length, cR1Length, Rght.B.y);
	    onesided_interp_scalar(cL0.fs->B.z, cR0.fs->B.z, cR1.fs->B.z,
				   cL0Length, cR0Length, cR1Length, Rght.B.z);
	}
	if ( G.turbulence_model == TM_K_OMEGA ) {
	    onesided_interp_scalar(cL0.fs->tke, cR0.fs->tke, cR1.fs->tke,
				   cL0Length, cR0Length, cR1Length, Rght.tke);
	    onesided_interp_scalar(cL0.fs->omega, cR0.fs->omega, cR1.fs->omega,
				   cL0Length, cR0Length, cR1Length, Rght.omega);
	}
        for ( int isp = 0; isp < nsp; ++isp ) {
	    onesided_interp_scalar(cL0.fs->gas->massf[isp],
				   cR0.fs->gas->massf[isp], cR1.fs->gas->massf[isp],
				   cL0Length, cR0Length, cR1Length, Rght.gas->massf[isp]);
        }
	
	// Make the thermodynamic properties consistent.
	// Pressure, Local Speed of Sound and Temperature.
        // The value of 1 indicates that the old temperature
        // should be used as an initial guess for the iterative
        // EOS functions.
	// If the EOS call fouls up, just copy the cell data, low-order.
	if ( nsp > 1 ) {
	    if ( scale_mass_fractions( Rght.gas->massf ) != 0 ) {
		for ( size_t isp=0; isp<Rght.gas->massf.size(); ++isp )
		    cR0.fs->gas->massf[isp] = Rght.gas->massf[isp];
	    }
	}
	
	// Interpolate on two of the thermodynamic quantities, and fill
	// in the rest based on an EOS call.
	
	onesided_interp_scalar(gL0.rho, gR0.rho, gR1.rho,
	                       cL0Length, cR0Length, cR1Length, Rght.gas->rho);
	for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	    onesided_interp_scalar(gL0.e[i], gR0.e[i], gR1.e[i],
				   cL0Length, cR0Length, cR1Length, Rght.gas->e[i]);
	}
	int status = gmodel->eval_thermo_state_rhoe(*(Rght.gas));
	
	if ( status != SUCCESS ) {
	    // Rght state failed.
	    Rght.copy_values_from(*(cL0.fs));
	}			      
    } // end of high-order reconstruction
    return SUCCESS;
} // end of onesided_interp()

//-----------------------------------------------------------------------------


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
    const double ALPHA = 0.5;
    global_data &gd = *get_global_data_ptr();
    CFlowCondition *gsp = gd.gas_state[inflow_condition_id];
    Gas_model *gmodel = get_gas_model_ptr();
    FlowState fsL(gmodel);
    fsL.copy_values_from(*gsp);
    Gas_data &gL = *(fsL.gas);
    Gas_data &gR = *(IFaceR.fs->gas);
    double time_weight = gd.shock_fitting_speed_factor;
    if ( gd.shock_fitting_decay ) {
        time_weight = time_weight - time_weight*(exp(-gd.sim_time/gd.max_time) - 1.0) / 
                             (exp(-1.0) - 1.0);
    }
    
    if ( gd.sim_time >= gd.t_moving ) {
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
