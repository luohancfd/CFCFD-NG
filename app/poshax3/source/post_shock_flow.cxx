/** \file post_shock_flow.cxx
 *
 *  \brief Defintions for the post-shock flow classes.
 *
 *  \author Rowan J. Gollan
 *  \version 07-Jan-07 : ported from Python version
 *
 **/

#include <iostream>
#include <iomanip>

#include <string>
#include <valarray>
#include <sstream>

#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/util/source/useful.h"
#include "post_shock_flow.hh"

using namespace std;

double p2_p1(double M1, double g)
{
    return 1.0 + 2.0 * g / (g + 1.0) * (M1*M1 - 1.0);
}

double r2_r1(double M1, double g)
{
    double numer = (g + 1.0) * M1*M1;
    double denom = 2.0 + (g - 1.0) * M1*M1;
    return numer / denom;
}

double T2_T1(double M1, double g)
{
    return  p2_p1(M1, g) / r2_r1(M1, g);
}

double u2_u1(double M1, double g)
{
    return 1.0 / r2_r1(M1, g);
}

Post_shock_flow::
Post_shock_flow() {}

Post_shock_flow::
Post_shock_flow( Flow_state &ic, Gas_model * gm, Reaction_update * ru, 
    		 Energy_exchange_update * eeu,
    		 PoshaxRadiationTransportModel * rtm )
: gmodel_( gm ), rupdate_( ru ), eeupdate_( eeu ), rtmodel_( rtm )
{
    // gas-model dimensions
    nsp_ = gmodel_->get_number_of_species();
    ntm_ = gmodel_->get_number_of_modes();
    
    // Set icflow
    icflow.set_flow_state(*ic.Q,ic.u);
    
    // make a new gas-data structure with the initial conditions for the zero solver
    Gas_data Q(*ic.Q);
    
    //    Make some guesses about the flow behind shock
    //    based on normal shock relations (for ideal gas)
    //    then use Newton iterations for the the post-shock state.
    int status;
    double T_inf = icflow.Q->T[0];
    double rho_inf = icflow.Q->rho;
    double u_inf = icflow.u;
    double M_inf = u_inf / icflow.Q->a;
    double gamma = gmodel_->gamma(*icflow.Q,status);
    double T2 = T_inf * T2_T1(M_inf, gamma);
    double rho2 = rho_inf * r2_r1(M_inf, gamma);
    double u2 = u_inf * u2_u1(M_inf, gamma);
    
    valarray<double> yguess(3), yout(3);
    yguess[0] = rho2; yguess[1] = T2; yguess[2] = u2;
    
    FrozenConservationSystem con_sys(gm,Q,u_inf);
    zero_solver_.set_constants(3, 1.0e-6, 100, true);
    zero_solver_.solve(con_sys, yguess, yout);
    
    Q.rho = yout[0];
    Q.T[0] = yout[1];
    gmodel_->eval_thermo_state_rhoT( Q );
    double u_s = yout[2];

    psflow.set_flow_state(Q,u_s);

    // Initialise the ODE solver pieces
    dcdt_.resize( nsp_, 0.0 );
    dedt_.resize( ntm_, 0.0 );
    int ndim = nsp_ + 1 + ntm_;
    yin_.resize( ndim, 0.0 );
    yout_.resize( ndim, 0.0 );
    yguess_.resize( ndim, 0.0 );
    ode_solver_.set_constants( "poshax3 noneq ODE system", ndim,
    	                       "rkf", 20, 1.15, 1.0e-2, 0.333 );
    ode_sys_ptr_ = dynamic_cast<OdeSystem*>(this);
    
    // Initialise the conservation system pieces
    zero_solver_.set_constants(ndim, 1.0e-6, 100, true);
    con_sys_.initialise( gmodel_, *psflow.Q, psflow.u );
}

Post_shock_flow::
~Post_shock_flow() {}

double 
Post_shock_flow::
increment_in_space(double x, double delta_x)
{
    // 1. Integrate the ODE system in space
    double new_x = ode_solve(x, delta_x);
    
    // 2. Solve for the flow state
    zero_solve( yout_ );
    
    return new_x;
}

double 
Post_shock_flow::
ode_solve(double x, double delta_x)
{
    double dx = delta_x;
    double dx_suggest = dx;
    
    // 1.  Encode the flux quantities from the current flow-state
    con_sys_.encode_conserved( yin_, *psflow.Q, psflow.u );
    
    // 2. Submit to ODE solver
    int flag = ode_solver_.solve_over_interval( *(ode_sys_ptr_), 0.0, dx, 
    	                                        &dx_suggest, yin_, yout_ );
    
    if ( ! flag ) {
    	cout << "Post_shock_flow::ode_solve()" << endl
    	     << "The ODE solver failed with the following flow-state:" << endl
    	     << psflow.str( bool(rtmodel_) ) << endl
    	     << "Bailing out!" << endl;
    	exit( FAILURE );
    }
    
    return dx_suggest;
}

void
Post_shock_flow::
zero_solve( const valarray<double> &A )
{
    // 0. Set constants in the conservation system
    con_sys_.set_constants( A );
    
    // 1. Create a guess
    for ( int isp=0; isp<nsp_; ++isp )
    	yguess_[isp] = psflow.Q->rho * psflow.Q->massf[isp];
    for ( int itm=0; itm<ntm_; ++itm )
    	yguess_[nsp_+itm] = psflow.Q->T[itm];
    yguess_[nsp_+ntm_] = psflow.u;
    
    // 2. Solve
    zero_solver_.solve(con_sys_, yguess_, yout_);
    
    // 3. Map
    psflow.Q->rho = 0.0;
    for ( int isp=0; isp<nsp_; ++isp ) {
    	if ( yout_[isp] < 0.0 ) yout_[isp] = 0.0;
    	psflow.Q->rho += yout_[isp];
    }
    for ( int isp=0; isp<nsp_; ++isp )
    	psflow.Q->massf[isp] = yout_[isp] / psflow.Q->rho;
    for ( int itm=0; itm<ntm_; ++itm )
    	psflow.Q->T[itm] = yout_[nsp_+itm];
    psflow.u = yout_[nsp_+ntm_];

    // 4. Rescale the mass-fractions to sum to one
    scale_mass_fractions( psflow.Q->massf );
    
    // 5. Equation of state evaluation
    gmodel_->eval_thermo_state_rhoT(*psflow.Q);
}

int
Post_shock_flow::
eval( const valarray<double> &y, valarray<double> &ydot )
{
    // NOTE: may need to omit this step if its too slow
    // 1. Evaluate the flow state from the y vector
    zero_solve( y );
    
    // 2. Fill out the ydot vector
    int iy=0;
    // 2a. Species mass production
    if ( rupdate_ ) rupdate_->rate_of_change( *psflow.Q, dcdt_ );
    for ( int isp=0; isp<nsp_; ++isp ) {
	ydot[iy] = dcdt_[isp] * gmodel_->molecular_weight(isp);
	++iy;
    }
    // 2b. Total momentum
    ydot[iy] = 0.0;
    ++iy;
    // 2c. Total energy
    psflow.Q_rad = 0.0;
    if ( rtmodel_ ) psflow.Q_rad = rtmodel_->eval_Q_rad( *psflow.Q );
    ydot[iy] = psflow.Q_rad;
    ++iy;
    // 2d. Modal energies
    if ( eeupdate_ ) eeupdate_->rate_of_change( *psflow.Q, dedt_ );
    if ( rupdate_ )
    	rupdate_->eval_chemistry_energy_coupling_source_terms( *psflow.Q, dedt_ );
    for ( int itm=1; itm<ntm_; ++itm ) {
    	ydot[iy] = dedt_[itm] * psflow.Q->rho;
    	if ( itm==(ntm_-1) ) ydot[iy] += psflow.Q_rad;
    	++iy;
    }
    
    // print_valarray( ydot );
    
    return 0;
}

