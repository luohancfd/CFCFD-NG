/** \file reaction_scheme.cxx
 *  \brief Definitions for the ReactionScheme class (and derivatives)
 *  \ingroup libgas2
 *
 *  \author Rowan J Gollan
 *  \version 23-Feb-2006
 **/

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <valarray>

#include "../models/gas.hh"

#include "../../util/source/useful.h"
#include "../../util/source/config_parser.hh"
// #include "../../util/source/dbc_assert.hh"

#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "../../nm/source/ode_system.hh"
#include "../../nm/source/ode_solver.hh"

#include "reaction.hh"
#include "reacting_species.hh"
#include "reaction_pieces.hh"
#include "reaction_scheme.hh"
#include "species_pieces.hh"

using namespace std;

#define RETRY_WARNING 0

ReactionScheme::ReactionScheme( const string name, const vector<Reaction*> reactions,
				int nsp, double err_tol, double T_trigger )
    : name_( name ), err_tol_( err_tol ), T_trigger_( T_trigger )
{
    g_ = get_gas_model_pointer();
    if( get_number_of_vibrational_species() >= 1 ) {
	pgm_ = new ThermallyPerfectGasMix("pgm", "mTg.pgm");
    }
    else {
	pgm_ = dynamic_cast<ThermallyPerfectGasMix*>(g_);
    }

    for( size_t i = 0; i < reactions.size(); ++i ) {
	reactions_.push_back( reactions[i]->clone() );
    }
    
    update_reacting_species();

}

ReactionScheme::ReactionScheme( const string name,
				int nsp, double err_tol, double T_trigger )
    : name_( name ), err_tol_( err_tol ), T_trigger_( T_trigger )
{
    g_ = get_gas_model_pointer();
    if( get_number_of_vibrational_species() >= 1 ) {
	pgm_ = new ThermallyPerfectGasMix("pgm", "mTg.pgm");
    }
    else {
	pgm_ = dynamic_cast<ThermallyPerfectGasMix*>(g_);
    }
    reactions_.clear();
    
}

ReactionScheme::ReactionScheme( const string name, const string input_file )
    : name_( name ) 
{
    ConfigParser cfg( input_file );
    int nr;
    if( ! cfg.parse_int( "reaction_scheme", "number_of_reactions", nr, 0 ) ) {
	cout << "ReactionScheme::ReactionScheme() --- \n";
	cout << "Error reading number_of_reactions in reaction_scheme section of: "
	     << input_file << endl;
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    return;
}

ReactionScheme::ReactionScheme( const ReactionScheme &r )
    : name_( r.name_ ), err_tol_( r.err_tol_ ), T_trigger_( r.T_trigger_ )
{
    g_ = get_gas_model_pointer();
    if( get_number_of_vibrational_species() >= 1 ) {
	pgm_ = new ThermallyPerfectGasMix("pgm", "mTg.pgm");
    }
    else {
	pgm_ = dynamic_cast<ThermallyPerfectGasMix*>(g_);
    }
    for( size_t i = 0; i < r.reactions_.size(); ++i ) {
	reactions_.push_back( r.reactions_[i]->clone() );
    }
    for( size_t i = 0; i < r.reacting_species_.size(); ++i ) {
	reacting_species_.push_back( r.reacting_species_[i]->clone() );
    }

    update_reacting_species();
	
}

ReactionScheme::~ReactionScheme()
{
    for( size_t i = 0; i < reactions_.size(); ++i ) {
	delete reactions_[i];
    }
    for( size_t i = 0; i < reacting_species_.size(); ++i ) {
	delete reacting_species_[i];
    }
    if( get_number_of_vibrational_species() >= 1 ) {
	delete pgm_;
    }
}

string ReactionScheme::str() const
{
    ostringstream ost;

    for( size_t i = 0; i < reactions_.size(); ++i ) {
	ost << "[reac-" << i << "]\n"
	    << reactions_[i]->str()
	    << endl;
    }
    return ost.str();
}

void ReactionScheme::add_reaction( Reaction *r )
{
    reactions_.push_back( r->clone() );
    update_reacting_species();
}

void ReactionScheme::update_reacting_species()
{
    // Assemble the reacting species
    reacting_species_.clear();

    vector<ReactionPieces*> r_pieces;

    for( int isp = 0; isp < g_->nsp; ++isp ) {
	for( size_t ir = 0; ir < reactions_.size(); ++ir ) {
	    int nu = reactions_[ir]->get_nu( isp );
	    if ( nu == 0 )
		continue;
	    if ( nu >= 1 ) {
		r_pieces.push_back( new ReactionPiecesA( ir, nu ) );
	    }
	    else {
		r_pieces.push_back( new ReactionPiecesB( ir, nu ) );
	    }
	}

	// Now we can create the ReactingSpecies
	reacting_species_.push_back( new ReactingSpecies( r_pieces ) );
	// Clean up the allocated pieces
	for( size_t i = 0; i < r_pieces.size(); ++i ) {
	    delete r_pieces[i];
	}
	r_pieces.clear();
    }
}

double ReactionScheme::get_k_f_for_reaction(int ir, Gas_data &Q)
{
    reactions_[ir]->fill_in_forward_coeff(Q);
    return reactions_[ir]->k_f();
}

ReactionSchemeODE::ReactionSchemeODE( const string name, vector<Reaction*> reactions,
				      int nsp, double err_tol, double T_trigger,
				      OdeSolver *ode_solver)
    : ReactionScheme( name, reactions, nsp, err_tol, T_trigger ),
      OdeSystem( nsp, true ), Q_( 0 )
{

    ode_solver_ = ode_solver->clone();

    w_f_.resize( reactions_.size() );
    w_b_.resize( reactions_.size() );
    
    yin_.resize( g_->nsp );
    yout_.resize( g_->nsp );
    q_.resize( g_->nsp );
    L_.resize( g_->nsp );

    mf_.resize(ndim_);
    c_.resize(ndim_);
    called_at_least_once = false;

}

ReactionSchemeODE::ReactionSchemeODE( const string name,
				      int nsp, double err_tol, double T_trigger,
				      OdeSolver *ode_solver)
    : ReactionScheme( name, nsp, err_tol, T_trigger ),
      OdeSystem( nsp, true ), Q_( 0 )
{

    ode_solver_ = ode_solver->clone();

    w_f_.resize( reactions_.size() );
    w_b_.resize( reactions_.size() );

    yin_.resize( g_->nsp );
    yout_.resize( g_->nsp );
    q_.resize( g_->nsp );
    L_.resize( g_->nsp );

    mf_.resize(ndim_);
    c_.resize(ndim_);
    called_at_least_once = false;

}

ReactionSchemeODE::ReactionSchemeODE( const ReactionSchemeODE &r )
    : ReactionScheme( r.name_, r.reactions_, r.ndim_, r.err_tol_, r.T_trigger_ ),
      OdeSystem( r.ndim_, true ), Q_( 0 )
{
    ode_solver_ = r.ode_solver_->clone();

    w_f_.resize( r.reactions_.size() );
    w_b_.resize( r.reactions_.size() );
    
    yin_.resize( g_->nsp );
    yout_.resize( g_->nsp );
    q_.resize( g_->nsp );
    L_.resize( g_->nsp );

    mf_.resize(r.ndim_);
    c_.resize(r.ndim_);
    called_at_least_once = false;
}

ReactionSchemeODE::~ReactionSchemeODE()
{
    delete ode_solver_;
}

ReactionSchemeODE*
ReactionSchemeODE::
clone()
{
    return new ReactionSchemeODE(*this);
}


string
ReactionSchemeODE::
str() const
{
    return string("");
	
}

void ReactionSchemeODE::update_reacting_species()
{
    ReactionScheme::update_reacting_species();
    
    w_f_.resize( reactions_.size() );
    w_b_.resize( reactions_.size() );
}

int
ReactionSchemeODE::
update_gas_state( Gas_data &Q, double dt_flow, bool include_evib_exchange )
{
    // 0. We can keep moving if the temperature
    //    is too cold for chemical reactions (set by user)
    
    if( Q.T <= T_trigger_ ) {
	Q.dt_chem = -1.0;
	return SUCCESS;
    }

    // 1. Setup for solution
    Q_ = &Q;  // Point to the current Gas_data struture:
              // this allows access in the OdeSystem eval functions

    for( int isp = 0; isp < g_->nsp; ++isp ) {
	yin_[isp] = Q.c[isp];
    }
    
    // Reaction rate coefficients have not been evaluated
    called_at_least_once = false;

    double h = Q.dt_chem;
    bool flag = false;

    // 2. Solve the system

    if( h > 0.0 ) {  // then we have a guess for the timestep
	flag = ode_solver_->solve_over_interval( *this, 0.0, dt_flow, &h,
						 yin_, yout_ );
	if( ! flag ) {
	    // then we retry with the timestep selected by our function
	    h = ReactionSchemeODE::stepsize_select( yin_ );
	    if( h > (0.5 * Q.dt_chem) ) {
		// If we're going to reduce the timestep, it's probably
		// best to do so drastically.  Anything less than half
		// what we started with is not really worth it.
		// Let's just try 10% of first timestep size
		h = 0.1 * Q.dt_chem ;
	    }
	    flag = ode_solver_->solve_over_interval( *this, 0.0, dt_flow, &h,
						     yin_, yout_ );
	}

	if( ! flag ) {
	    // cout << "ReactionSchemeODE::update_gas_state()\n"
// 		 << "The ode solver has failed twice:\n"
// 		 << "   1. with previously used dt_chem, and\n"
// 		 << "   2. with dt_chem based on stepsize selection algorithm.\n"
// 		 << "The following gas state was not solved over dt_flow= " << dt_flow << endl
// 		 << "with input dt_chem= " << Q.dt_chem << endl;
// 	    Q.print_values();
	    return( NUMERICAL_ERROR );
	}
    }
    else { // it's probably our first step (or after T_trigger invocation)	
	 h = ReactionSchemeODE::stepsize_select( yin_ );
	 flag = ode_solver_->solve_over_interval( *this, 0.0, dt_flow, &h,
						 yin_, yout_ );

	 if( ! flag ) {
	    // cout << "ReactionSchemeODE::update_gas_state()\n"
// 		 << "The ode solver has failed with dt_chem based on stepsize selection algorithm.\n"
// 		 << "The following gas state was not solved over dt_flow= " << dt_flow << endl
// 	    	 << "with input dt_chem= " << h << endl;
// 	    Q.print_values();
	    return( NUMERICAL_ERROR );
	}

    }

    // 3. If we've made it this far than we're doing well.
    //    Let's update the gas state and leave.
    

    for( int isp = 0; isp < g_->nsp; ++isp ) {
	Q.c[isp] = yout_[isp];
    }

    fill_in_mass_fractions( Q.rho, Q.c, g_->mol_weight, Q.f );

    // Normalise the mass fractions
    double f_sum = 0.0;
    for( int isp = 0; isp < g_->nsp; ++isp) {
	Q.f[isp] = Q.f[isp] >= 0.0 ? Q.f[isp] : 0.0;
	f_sum += Q.f[isp];
    }
    for( int isp = 0; isp < g_->nsp; ++isp) Q.f[isp] /= f_sum;

    Q.dt_chem = h;
    return SUCCESS;
}

int
ReactionSchemeODE::
eval_split( const valarray<double> &y, 
	    valarray<double> &q, valarray<double> &L )
{
    if( ! called_at_least_once ) { // We need to compute the rate coefficients
	for( size_t ir = 0; ir < reactions_.size(); ++ir ) {
	    if( reactions_[ir]->compute_kf_followed_by_kb ) {
		reactions_[ir]->fill_in_forward_coeff( (*Q_) );
		reactions_[ir]->fill_in_backward_coeff( (*Q_) );
	    }
	    else {
		reactions_[ir]->fill_in_backward_coeff( (*Q_) );
		reactions_[ir]->fill_in_backward_coeff( (*Q_) );
	    }
	}
	called_at_least_once = true;
    }

    // Compute all the forward and backward production rates
    // for individual reactions.

    for( size_t ir = 0; ir < reactions_.size(); ++ir) {
	w_f_[ir] = reactions_[ir]->forward_rate( y );
	w_b_[ir] = reactions_[ir]->backward_rate( y );
	//	cout << "ir= " << ir << " w_f= " << _w_f[ir] << " w_b= " << _w_b[ir] << endl;
    }


    // Now collect the production rates for the reactions
    // and apply them to the species.
    for( size_t isp = 0; isp < reacting_species_.size(); ++isp ) {
	q[isp] = reacting_species_[isp]->production( w_f_, w_b_ );
	L[isp] = reacting_species_[isp]->loss( w_f_, w_b_ );
    }

    return 0;
}

int
ReactionSchemeODE::
eval( const valarray<double> &y, valarray<double> &ydot )
{
    
    ReactionSchemeODE::eval_split( y, q_, L_ );
    
    for( size_t isp = 0; isp < y.size(); ++isp ) {
	ydot[isp] = q_[isp] - L_[isp];
    }
    return 0;
}

const double eps1 = 0.001;
const double chem_step_upper_limit = 1.0e-3;
const double chem_step_lower_limit = 1.0e-20;
const double zero_tol = 1.0e-30;

double
ReactionSchemeODE::
stepsize_select( const valarray<double> &y )
{
    valarray<double> ydot = y;

    ReactionSchemeODE::eval( y, ydot );

    double min_dt = chem_step_upper_limit; // to get us started
    double old_dt = 0.0;

    for( int isp = 0; isp < ndim_; ++isp ) {
	if( (y[isp] > 0.0) && (fabs(ydot[isp]) > zero_tol) ) {
	    old_dt = fabs( y[isp] / ydot[isp] );
	    if( old_dt < min_dt ) {
		min_dt = old_dt;
	    }
	}
    }

    double dt_chem = eps1 * min_dt;

    // Impose upper and lower chem_step limits
    if( dt_chem > chem_step_upper_limit )
	dt_chem = chem_step_upper_limit;
    else if ( dt_chem < chem_step_lower_limit )
	dt_chem = chem_step_lower_limit;

    return dt_chem;

}

const double min_conc = 1.0e-30;

bool
ReactionSchemeODE::
passes_system_test( valarray<double> &y )
{
    double f_tot = 0.0;
    for( int isp = 0; isp < g_->nsp; ++isp ) {
	// cout << "isp: " << isp << " y= " << y[isp];
	y[isp] = y[isp] < min_conc ? 0.0 : y[isp];
	c_[isp] = y[isp];
	// cout << " c_= " << c_[isp] << endl;
    }
	 
    f_tot = fill_in_mass_fractions( Q_->rho, c_, g_->mol_weight, mf_ );
   //  cout << "f_tot= " << f_tot << endl;

    double lower_lim = 1.0 - ( err_tol_ / 100.0 );
    double upper_lim = 1.0 + ( err_tol_ / 100.0 );

    if( (f_tot < lower_lim) || (f_tot > upper_lim) )
	return false;
    else {
	// We can normalise before moving on.
	if( f_tot != 1.0) {
	    for( int isp = 0; isp < g_->nsp; ++isp) mf_[isp] /= f_tot;
	    fill_in_concentrations( Q_->rho, mf_, g_->mol_weight, c_ );
	    for( int isp = 0; isp < g_->nsp; ++isp) y[isp] = c_[isp];
	}

	return true;
    }

}

void
ReactionSchemeODE::evolution_rates(int isp, valarray<double> &W, Gas_data &Q)
{
    if( W.size() < size_t(reacting_species_[isp]->nreac()) ) {
	cerr << "ReactionSchemeODE::evolution_rates()\n"
	     << "Not enough space to store results of species evolution for each reaction.\n"
	     << "Increase size of valarray<double> W.\n"
	     << "Bailing out!\n";
	exit(MISMATCHED_DIMENSIONS);
    }


    for( size_t ir = 0; ir < reactions_.size(); ++ir ) {
	if( reactions_[ir]->compute_kf_followed_by_kb ) {
	    reactions_[ir]->fill_in_forward_coeff( Q );
	    reactions_[ir]->fill_in_backward_coeff( Q );
	}
	else {
	    reactions_[ir]->fill_in_backward_coeff( Q );
	    reactions_[ir]->fill_in_backward_coeff( Q );
	}
    }


    // Compute all the forward and backward production rates
    // for individual reactions.

    valarray<double> y( ndim_ );
    for( int i = 0; i < ndim_; ++i ) {
	y[i] = Q.c[i];
    }
    for( size_t ir = 0; ir < reactions_.size(); ++ir) {
	w_f_[ir] = reactions_[ir]->forward_rate( y );
	w_b_[ir] = reactions_[ir]->backward_rate( y );
    }

    for( int ir = 0; ir < reacting_species_[isp]->nreac(); ++ir ) {
	W[ir] = reacting_species_[isp]->rate_for_reac( ir, w_f_, w_b_ );
    }

}


// -----------------------------------------------------------
// Functions that provide access to the "managed" finite-rate
// chemistry module.

static ReactionScheme *rscheme = 0;

/// \brief Selects and sets up the managed ReactionScheme
int set_reaction_scheme( const string name,
			 const string input_file)
{
    ConfigParser cfg( input_file );
    int nr = 0;

    if( ! cfg.parse_int( "reaction_scheme", "number_of_reactions",  nr, 0 ) ) {
	cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
	     << "Error reading number_of_reactions in [reaction_scheme] section of: "
	     << input_file << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( nr <= 0 ) {
	cerr << "set_reaction_scheme() --- \n"
	     << "Error for number of reactions, nr= " << nr << endl
	     << "A positive integer (> 0) expected.\n"
	     << "Check number_of_reactions in section [reaction_scheme] of file: " 
	     << input_file << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    vector<Reaction*> reactions;
    Reaction *r_pointer;
    for( int ir = 0; ir < nr; ++ir ) {
	r_pointer = initialize_reaction( cfg, ir );
	if( r_pointer != 0 ) {
	    reactions.push_back( r_pointer->clone() );
	    delete r_pointer;
	}
	else {
	    cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
		 << "Error initializing reaction: " << ir << endl
		 << "Check input file: " << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
    }
    
    // Gather the extra information.
    double err_tol = 0.0;
    if( ! cfg.parse_double( "reaction_scheme", "err_tol", err_tol, 0.1 ) ) {
	cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
	     << "Error reading err_tol in [reaction_scheme] section of: "
	     << input_file << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    
    double T_trigger = 0.0;
    if( ! cfg.parse_double( "reaction_scheme", "T_trigger", T_trigger, 0.0 ) ) {
	cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
	     << "Error reading T_trigger in [reaction_scheme] section of: "
	     << input_file << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    GasModel *g = get_gas_model_pointer();

    string method;
    if( ! cfg.parse_string( "reaction_scheme", "method", method, "none" ) ) {
	cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
	     << "Error reading method in [reaction_scheme] section of: "
	     << input_file << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( method == "ode" ) {
	// Read in values
	string step_routine;
	if( ! cfg.parse_string( "ode_method", "step_routine", step_routine, "qss" ) ) {
	    cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
		 << "Error reading step_routine in [ode_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	int max_step_attempts = 0;
	if( ! cfg.parse_int( "ode_method", "max_step_attempts", max_step_attempts, 4 ) ) {
	    cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
		 << "Error reading max_step_attempts in [ode_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	double max_increase_factor = 0.0;
	if( ! cfg.parse_double( "ode_method", "max_increase_factor", max_increase_factor, 1.0 ) ) {
	    cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
		 << "Error reading max_increase_factor in [ode_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	double max_decrease_factor = 0.0;
	if( ! cfg.parse_double( "ode_method", "max_decrease_factor", max_decrease_factor, 0.01 ) ) {
	    cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
		 << "Error reading max_decrease_factor in [ode_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	double decrease_factor = 0.0;
	if( ! cfg.parse_double( "ode_method", "decrease_factor", decrease_factor, 0.01 ) ) {
	    cerr << "set_reaction_scheme() --- " << __FILE__ << (__LINE__ - 1) << endl
		 << "Error reading decrease_factor in [ode_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	OdeSolver ode_solver( "reacting system", g->nsp, step_routine,
			      max_step_attempts, max_increase_factor,
			      max_decrease_factor, decrease_factor );

	rscheme = new ReactionSchemeODE( name, reactions, g->nsp, err_tol, T_trigger,
					 &ode_solver );


    }
    else if ( method == "ode_mc" ) {
	string step_routine;
	if( ! cfg.parse_string( "ode_mc_method", "step_routine", step_routine, "rkf" ) ) {
	    cerr << "set_reaction_scheme() --- \n"
		 << "Error reading step_routine in [ode_mc_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	int max_step_attempts = 0;
	if( ! cfg.parse_int( "ode_mc_method", "max_step_attempts", max_step_attempts, 4 ) ) {
	    cerr << "set_reaction_scheme() --- \n"
		 << "Error reading max_step_attempts in [ode_mc_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	double max_increase_factor = 0.0;
	if( ! cfg.parse_double( "ode_mc_method", "max_increase_factor", max_increase_factor, 1.0 ) ) {
	    cerr << "set_reaction_scheme() --- \n"
		 << "Error reading max_increase_factor in [ode_mc_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	double max_decrease_factor = 0.0;
	if( ! cfg.parse_double( "ode_mc_method", "max_decrease_factor", max_decrease_factor, 0.01 ) ) {
	    cerr << "set_reaction_scheme() --- \n"
		 << "Error reading max_decrease_factor in [ode_mc_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	double decrease_factor = 0.0;
	if( ! cfg.parse_double( "ode_mc_method", "decrease_factor", decrease_factor, 0.01 ) ) {
	    cerr << "set_reaction_scheme() --- \n"
		 << "Error reading decrease_factor in [ode_mc_method] section of: "
		 << input_file << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	OdeSolver ode_solver( "reacting system (mass conserving)", 2*reactions.size(), step_routine,
			      max_step_attempts, max_increase_factor,
			      max_decrease_factor, decrease_factor );

	rscheme = new ReactionSchemeODE_MC( name, reactions, g->nsp, err_tol, T_trigger,
					    &ode_solver );
	// Now that the pointer is set, it's possible to initialise the
	// evib coupling component.
	ReactionSchemeODE_MC *rs = dynamic_cast<ReactionSchemeODE_MC*>(rscheme);
	rs->initialise_evib_coupling(cfg);
    }
    else {
	cerr << "set_reaction_scheme() --- " << endl
	     << "Error with selection of method in [reaction_scheme] section of: "
	     << input_file << endl
	     << "The method: " << method << " is undefined or not implemented.\n"
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    for( size_t ir = 0; ir < reactions.size(); ++ir ) {
	delete reactions[ir];
    }

    return SUCCESS;
}

void clear_reaction_scheme_pointer()
{
    delete rscheme;
    rscheme = 0;
}


ReactionScheme *get_reaction_scheme_pointer()
{
    return rscheme;
}

ThermallyPerfectGasMix* get_tpgm_pointer()
{
    return rscheme->get_tpgm_pointer();
}

int get_nu_for_reaction(int ir, int index)
{
    return rscheme->get_nu_for_reaction(ir, index);
}

double get_k_f_for_reaction(int ir, Gas_data &Q)
{
    return rscheme->get_k_f_for_reaction(ir, Q);
}

static Gas_data *Q_save;

int estimate_appropriate_subcycle(double dt_flow, double dt_chem, double *dt, int *no_substeps)
{
    const double allowable_t_factor = 2.0;
    const int default_substeps = 50;

    if( dt_chem < 0.0 ) {
	*dt = dt_flow/default_substeps;
	*no_substeps = default_substeps;
    }
    else {
	double target_t = allowable_t_factor * dt_chem;
	*no_substeps = int(dt_flow/target_t) + 1;
	*dt = dt_flow/(*no_substeps);
    }
    return SUCCESS;
}


int perform_chemical_increment(Gas_data &Q, double t_interval)
{
    int flag = SUCCESS;
    
    if( Q_save == NULL ) {
	Q_save = new Gas_data(get_gas_model_pointer());
    }

    if( rscheme == 0 ) {
	cerr << "Reaction scheme has not yet been initialized.\n";
	cerr << "Bailing out!\n";
	exit( NOT_INITIALIZED_ERROR );
    }

    // 1. Attempt the increment with evib exchange allowed.
    //    Subcycle if necessary.

    Q_save->copy_values_from(Q);
    if( rscheme->update_gas_state(Q, t_interval) == SUCCESS && EOS_rhoe(Q, true, false, false) == SUCCESS ) {
        // All is well.
        return SUCCESS;
    }
    else {
#       if RETRY_WARNING == 1
        cout << "Not succesful on attempt of perform_chemical_increment which\n";
	cout << "includes the transfer of vibrational energy.\n";
	cout << "Repeat with more frequent thermo updates.\n";
#       endif       
	Q.copy_values_from(*Q_save);

	double dt;
	int no_substeps;
	estimate_appropriate_subcycle(t_interval, Q.dt_chem, &dt, &no_substeps);
// 	cout << "Retrying chemical increment with...\n";
// 	cout << "dt= " << dt << " no_substeps= " << no_substeps << endl;
// 	cout << "Given t_interval was: " << t_interval << " and Q.dt_chem= " << Q.dt_chem << endl;

	for( int i = 0; i < no_substeps; ++i ) {

	    if( rscheme->update_gas_state(Q, dt) != SUCCESS || EOS_rhoe(Q, true, false, false) != SUCCESS ) {
#               if RETRY_WARNING == 1		
		cout << "Failed in perform chemical increment when including\n"
		     << "evib exchange.  Retry without transfer of vibrational energies.\n";
#               endif	    

		flag = FAILURE;

		//		cout << "t_interval= " << t_interval << endl;
		//		cout << "Gas state at beginning of increment: " << endl;
		//		print_gas_data(Q_save);
		    
		break;
	    }
	}
	if ( flag == SUCCESS )
	    return SUCCESS;
	// otherwise proceed to iterations WITHOUT evib exchange.
    }
    
    // 2. If we are still having problems, we might need
    //    to repeat without including the energy exchange.
    //    If only if we are dealing with multi-temperature flows.

    if( get_number_of_vibrational_species() == 0 ) {
	// If we got here in a single-temp sim, then we're in trouble.
	return FAILURE;
    }
    else {
	Q.copy_values_from(*Q_save);
	if( rscheme->update_gas_state(Q, t_interval, false) == SUCCESS && EOS_rhoe(Q, true, false, false) == SUCCESS ) {
	    // Things are OK if evib exhnage is left out.
	    return SUCCESS;
	}
	else {
#           if RETRY_WARNING == 1
	    cout << "Not succesful on attempt of perform_chemical_increment which\n";
	    cout << "ignores the transfer of vibrational energy.\n";
	    cout << "Repeat with more frequent thermo updates.\n";
#           endif       
	    Q.copy_values_from(*Q_save);

	    double dt;
	    int no_substeps;
	    estimate_appropriate_subcycle(t_interval, Q.dt_chem, &dt, &no_substeps);
	
	    for( int i = 0; i < no_substeps; ++i ) {
		if( rscheme->update_gas_state(Q, dt, false) != SUCCESS ||
		    EOS_rhoe(Q, true, false, false) != SUCCESS ) { 
		    cout << "Failed in perform chemical increment when ignoring\n"
			 << "evib exchange.\n";
		    cout << "The initial condition was: \n";
		    Q_save->print_values();
		    cout << "Bailing out!\n";
		    exit(NUMERICAL_ERROR);
		}
	    }
	}
    }
    // If we get this far...
    return SUCCESS;
}

ReactionSchemeODE* initialize_ReactionSchemeODE(const string name,
						const string input_file)
{
    ReactionSchemeODE *rs;
    
    if ( set_reaction_scheme( name, input_file) == SUCCESS ) {
	rs = dynamic_cast< ReactionSchemeODE* >(get_reaction_scheme_pointer());
	return rs;
    }
    // Otherwise return empty pointer
    return 0;
}


ReactionSchemeODE_MC::ReactionSchemeODE_MC(const string name, const vector<Reaction*> reactions,
					   int nsp, double err_tol, double T_trigger,
					   OdeSolver *ode_solver )
    : ReactionScheme( name, reactions, nsp, err_tol, T_trigger ),
      OdeSystem( 2*reactions.size(), true ), Q_( 0 )
{
    ode_solver_ = ode_solver->clone();

    yin_.resize( 2*reactions_.size() );
    yout_.resize( 2*reactions_.size() );
    w_.resize( reactions_.size() );
    c_.resize(g_->nsp);
    cinit_.resize(g_->nsp);
    
    vector<int> reac_index;
    vector<int> nu;

    for( int isp = 0; isp < g_->nsp; ++isp ) {
	for( size_t ir = 0; ir < reactions_.size(); ++ir ) {
	    int nu_test = reactions_[ir]->get_nu( isp );
	    if ( nu_test == 0 )
		continue;
	    else {
		reac_index.push_back( ir );
		nu.push_back( nu_test );
	    }
	}

	// Now we can create the ReactingSpecies
	spec_.push_back( new SpeciesPieces( isp, reac_index, nu ) );
	// Clean up the allocated pieces
	reac_index.clear();
	nu.clear();
    }

    called_at_least_once = false;

}

ReactionSchemeODE_MC::ReactionSchemeODE_MC( const ReactionSchemeODE_MC &r )
    : ReactionScheme( r.name_, r.reactions_, r.ndim_/2, r.err_tol_, r.T_trigger_ ),
      OdeSystem( r.ndim_, true ), Q_( 0 )
{
    ode_solver_ = r.ode_solver_->clone();

    yin_.resize(r.ndim_, 0.0);
    yout_.resize(r.ndim_, 0.0);
    w_.resize(reactions_.size(), 0.0);
    c_.resize(g_->nsp, 0.0);
    cinit_.resize(g_->nsp, 0.0);

    for( int isp = 0; isp < g_->nsp; ++isp )
	spec_.push_back( r.spec_[isp]->clone() );

    called_at_least_once = false;
}
    

ReactionSchemeODE_MC::~ReactionSchemeODE_MC()
{
    delete ode_solver_;
    for( int isp = 0; isp < g_->nsp; ++isp )
	delete spec_[isp];
    for( int ivib = 0; ivib < get_number_of_vibrational_species(); ++ivib)
	delete evc_[ivib];
}

ReactionSchemeODE_MC*
ReactionSchemeODE_MC::clone()
{

    return new ReactionSchemeODE_MC(*this);
}

string
ReactionSchemeODE_MC::str() const
{
    return string("");
}


int
ReactionSchemeODE_MC::
initialise_evib_coupling(ConfigParser &cfg)
{

    evc_.resize(get_number_of_vibrational_species());
    vector<int> r_list;
    vector<int> pvib;
    vector<int> not_found;
    not_found.push_back(-1);

    for( int ivib = 0; ivib < get_number_of_vibrational_species(); ++ivib) {
	r_list.clear();
	for( int ir = 0; ir < int(reactions_.size()); ++ir) {
	    ostringstream section;
	    section << "reac-" << ir;
	    cfg.parse_vector_of_ints(section.str(), "ivib_list", pvib, not_found);

	    for( size_t i = 0; i < pvib.size(); ++i ) {
		if( pvib[i] == ivib )
		    r_list.push_back(ir);
	    }
	}
	evc_[ivib] = new Evib_coupling(cfg, r_list, ivib);
    }
	
    return 0;
}

int
ReactionSchemeODE_MC::
update_gas_state( Gas_data &Q, double dt_flow, bool include_evib_exchange )
{
    // 0. We can keep moving if the temperature
    //    is too cold for chemical reactions (set by user)

    if( Q.T <= T_trigger_ ) {
	Q.dt_chem = -1.0;
	return SUCCESS;
    }

     // 1. Setup for solution
    Q_ = &Q;  // Point to the current gas_data struture:
              // this allows access in the OdeSystem eval functions

    for( int i = 0; i < ndim_; ++i ) yin_[i] = 0.0;

    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	spec_[isp]->set_init_conc( Q.c[isp] );
	cinit_[isp] = Q.c[isp];
    }

    // Reaction rate coefficients have NOT been evaluated.
    called_at_least_once = false;
    
    double h = Q.dt_chem;
    bool flag = false;

    // 2. Solve the system

    if( h > 0.0 ) {  // then we have a guess for the timestep
	flag = ode_solver_->solve_over_interval( *this, 0.0, dt_flow, &h,
						 yin_, yout_ );
	if( ! flag ) {
	    // then we retry with the timestep selected by our function
	    h = ReactionSchemeODE_MC::stepsize_select( yin_ );
	    if( h > (0.5 * Q.dt_chem) ) {
		// If we're going to reduce the timestep, it's probably
		// best to do so drastically.  Anything less than half
		// what we started with is not really worth it.
		// Let's just try 10% of first timestep size
		h = 0.1 * Q.dt_chem ;
	    }
	    flag = ode_solver_->solve_over_interval( *this, 0.0, dt_flow, &h,
						     yin_, yout_ );
	}

	if( ! flag ) {
//	    cerr << "ReactionSchemeODE_MC::update_gas_state()\n"
//		 << "The ode solver has failed twice:\n"
//		 << "   1. with previously used dt_chem, and\n"
//		 << "   2. with dt_chem based on stepsize selection algorithm.\n"
//		 << "The following gas state was not solved over dt_flow= " << dt_flow << endl
//		 << "with input dt_chem= " << Q.dt_chem << endl;
//	    print_gas_data( Q );
	    return( NUMERICAL_ERROR );
	}
    }
    else { // it's probably our first step (or after T_trigger invocation)
	 h = ReactionSchemeODE_MC::stepsize_select( yin_ );
	 flag = ode_solver_->solve_over_interval( *this, 0.0, dt_flow, &h,
						  yin_, yout_ );

	 if( ! flag ) {
	    cerr << "ReactionSchemeODE_MC::update_gas_state()\n"
		 << "The ode solver has failed with dt_chem based on stepsize selection algorithm.\n"
		 << "The following gas state was not solved over dt_flow= " << dt_flow << endl
		 << "with input dt_chem= " << h << endl;
	    print_gas_data( Q );
	    exit( NUMERICAL_ERROR );
	}

    }

    // 3. Assemble the solution

    for( int ir = 0; ir < ndim_; ir += 2 ) {
	w_[ir/2] = yout_[ir] - yout_[ir+1];
    }

    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	Q.c[isp] = spec_[isp]->eval_conc( w_ );
    }

    fill_in_mass_fractions( Q.rho, Q.c, g_->mol_weight, Q.f );

    // 4. Update vibrational energy change (if appropriate)
    if ( include_evib_exchange ) {
	for( int ivib = 0; ivib < get_number_of_vibrational_species(); ++ivib ) {
	    Q.e_vib[ivib] = evc_[ivib]->compute_new_e_vib(Q, yout_, cinit_, Q.c);
	    // ASSERT( !isnan(Q.e_vib[ivib]) && !isinf(Q.e_vib[ivib]) );
	}
	if( separate_electron_energy() && !get_tTg_flag() ) {
	    Q.T_e = Q.T_vib[0];
	    Q.e_e = compute_e_elec(Q);
	}
	else if ( get_tTg_flag() ) {
	    // Sometimes the electronic energy is lost (?), needs to be re-evaluated
	    Q.e_e = compute_e_elec(Q);
	}
    }
    Q.dt_chem = h;

    return SUCCESS;

}

int
ReactionSchemeODE_MC::
eval( const valarray<double> &y, valarray<double> &ydot )
{
    if( ! called_at_least_once ) { // We need to compute the rate coefficients
	for( size_t ir = 0; ir < reactions_.size(); ++ir ) {
	    if( reactions_[ir]->compute_kf_followed_by_kb ) {
		reactions_[ir]->fill_in_forward_coeff( (*Q_) );
		reactions_[ir]->fill_in_backward_coeff( (*Q_) );
	    }
	    else {
		reactions_[ir]->fill_in_backward_coeff( (*Q_) );
		reactions_[ir]->fill_in_backward_coeff( (*Q_) );
	    }
	}
	called_at_least_once = true;
    }

    // Calculate w

    for( int ir = 0; ir < ndim_; ir += 2 ) {
	w_[ir/2] = y[ir] - y[ir+1];
    }

    // Now calculate concentration vector
    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	c_[isp] = spec_[isp]->eval_conc( w_ );
    }

    // Finally assemble the rate.
    for( int ir = 0; ir < ndim_; ir += 2 ) {
	ydot[ir] = reactions_[ir/2]->forward_rate( c_ );
	ydot[ir+1] = reactions_[ir/2]->backward_rate( c_ );
    }

    return 0;
}


double
ReactionSchemeODE_MC::
stepsize_select( const valarray<double> &y )
{
    valarray<double> ydot = y;
    double rate = 0.0;

    ReactionSchemeODE_MC::eval( y, ydot );

    double min_dt = chem_step_upper_limit; // to get us started
    double old_dt = 0.0;

    for( size_t isp = 0; isp < spec_.size(); ++isp ) {
	rate = 0.0;
	for( size_t i = 0; i < spec_[isp]->nu_.size(); ++i ) {
	    rate += spec_[isp]->nu_[i] * ( ydot[2*i] - ydot[2*i + 1] );
	}

	if( (spec_[isp]->get_init_conc() > 0.0) && (fabs(rate) > zero_tol) ) {
	    old_dt = fabs( spec_[isp]->get_init_conc() / rate );
	    if( old_dt < min_dt ) {
		min_dt = old_dt;
	    }
	}
    }

    double dt_chem = eps1 * min_dt;

    // Impose upper and lower chem_step limits
    if( dt_chem > chem_step_upper_limit )
	dt_chem = chem_step_upper_limit;
    else if ( dt_chem < chem_step_lower_limit )
	dt_chem = chem_step_lower_limit;

    return dt_chem;

}


    
bool
ReactionSchemeODE_MC::
passes_system_test( valarray<double> &y )
{
    // We should never have a problem 
    return true;
}


int
ReactionSchemeODE_MC::
 update_evib(gas_data &Q, valarray<double> &y,
 	    valarray<double> &c_old, valarray<double> &c_new)
{
    for(size_t ivib = 0; ivib < evc_.size(); ++ivib) {
	Q.e_vib[ivib] = evc_[ivib]->compute_new_e_vib(Q, y, c_old, c_new);
    }
    if( separate_electron_energy() ) {
	Q.T_e = Q.T_vib[0];
	Q.e_e = compute_e_elec(Q);
    }
    return 0;
}
    




