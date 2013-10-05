// Author: Daniel F. Potter
// Date: 18-Nov-2009

#include <iostream>
#include <sstream>
#include <numeric>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "energy-exchange-ODE-update.hh"
#include "ode_setup.hh"

#define PERSISTENT_EE_UPDATE 0
#define ALLOW_EE_SUBCYCLING 1

using namespace std;

Energy_exchange_ODE_update::
Energy_exchange_ODE_update(lua_State *L, Gas_model &g)
 : g_( &g )
{
    lua_getglobal(L, "scheme_t");
    lua_getfield(L, -1, "temperature_limits");
    T_lower_limit_ = get_positive_number(L, -1, "lower");
    T_upper_limit_ = get_positive_number(L, -1, "upper");
    lua_pop(L, 1);

    double error_tol = get_positive_number(L, -1, "error_tolerance");

    lua_pop(L, 1); // pop scheme_t
  
    if ( T_upper_limit_ <= T_lower_limit_ ) {
	ostringstream ost;
	ost << "Energy_exchange_ODE_update::Energy_exchange_ODE_update():\n";
	ost << "Error in input: T_upper_limit must be greater than T_lower_limit.\n";
	ost << "T_upper_limit= " << T_upper_limit_ << " T_lower_limit= " << T_lower_limit_ << endl;
	input_error(ost);
    }
    
    // We don't include the total energy in the energy exchange system
    int nmodes_m1 = g_->get_number_of_modes() - 1;
    if ( nmodes_m1==0 ) {
	ostringstream ost;
	ost << "Energy_exchange_ODE_update::Energy_exchange_ODE_update():\n";
	ost << "Error in input: gas-model considers only one thermal mode.\n";
	input_error(ost);
    }

    lua_getglobal(L, "ode_t");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Energy_exchange_ODE_update::Energy_exchange_ODE_update():\n";
	ost << "Error in input: ode_solver table is missing.\n";
	input_error(ost);
    }
    ode_solver_ = create_ode_solver(L, nmodes_m1, "energy exchange ODE update system");
    lua_pop(L, 1);

    ees_ = new Energy_exchange_system(L, g, error_tol);

    yin_.resize(nmodes_m1, 0.0);
    yout_.resize(nmodes_m1, 0.0);
    ydot_.resize(nmodes_m1, 0.0);
    Q_save_ = new Gas_data(g_);
    molef_.resize( g_->get_number_of_species() );
    
    // Create equilibrium mechanisms
    lua_getglobal(L, "equilibriation_mechanisms");
    if ( lua_istable(L, -1) ) {
	for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i);
	    eq_mechs_.push_back( new Thermal_equilibrium_mechanism(L) );
	    lua_pop(L, 1);
	}
    }
    lua_pop(L, 1);
}

Energy_exchange_ODE_update::
~Energy_exchange_ODE_update()
{
    delete ode_solver_;
    delete ees_;
    delete Q_save_;
    yin_.resize(0);
    yout_.resize(0);
    ydot_.resize(0);
    molef_.resize(0);
}

int
Energy_exchange_ODE_update::
s_update_state(Gas_data &Q, double t_interval, double &dt_suggest, Gas_model *gm)
{
#   if ALLOW_EE_SUBCYCLING
    int flag = SUCCESS;
#   endif
    
    if ( Q.T[0] <= T_lower_limit_  || Q.T[0] >= T_upper_limit_ ) {
	dt_suggest = -1.0;
	return SUCCESS;
    }

    // Keep a copy in case something goes wrong
    // and we need to retry
    Q_save_->copy_values_from(Q);
#   if ALLOW_EE_SUBCYCLING
    double dt_suggest_save = dt_suggest;
#   endif
    
    // 0. Compute mole-fractions
    convert_massf2molef( Q.massf, g_->M(), molef_ );

    // 1. Attempt to solve normally.
    if ( perform_increment(Q, t_interval, dt_suggest) == SUCCESS ) {
	// all is well
	return SUCCESS;
    }
#   if ALLOW_EE_SUBCYCLING == 1
    else { // Repeat the attempt by subcycling.
	Q.copy_values_from(*Q_save_);
	double dt_sub;
	int no_substeps;
	estimate_appropriate_subcycle(t_interval, dt_suggest_save,
				      dt_sub, no_substeps);
	cout << "Energy_exchange_ODE_update::s_update_state()" << endl
	     << "repeating the attempt by subcycling: no_substeps = " << no_substeps << endl;
	
	for( int i = 0; i < no_substeps; ++i ) {
	    // Update the gas-state assuming constant density and energy
	    if ( i > 0 && gm!=0 ) {
	    	gm->eval_thermo_state_rhoe( Q );
	    }
	    if ( perform_increment(Q, dt_sub, dt_suggest) != SUCCESS ) {
		flag = FAILURE;
		break;
	    }
	}
	if ( flag == SUCCESS )
	    return SUCCESS;
	else {
	    cout << "Failed to successfully update gas state due to thermal energy exchange.\n";
#           if PERSISTENT_EE_UPDATE == 1
            cout << "Freezing the gas state and attempting to continue.\n";
            Q_save_->copy_values_from( Q );
            return SUCCESS;
#           else
	    cout << "The initial condition was: \n";
	    Q_save_->print_values(false);
	    cout << "t_interval = " << dt_sub << ", dt_suggest = " << dt_suggest << endl;
	    cout << "Bailing Out!\n";
	    exit(NUMERICAL_ERROR);
#           endif
	}
    }
#   else
    else {
	cout << "Failed to successfully update gas state due to thermal energy exchange.\n";
#       if PERSISTENT_EE_UPDATE == 1
	cout << "Freezing the gas state and attempting to continue (dt_suggest = " << dt_suggest << ").\n";
	Q.copy_values_from(*Q_save_);
	return SUCCESS;
#       else
	cout << "The initial condition was: \n";
	Q_save_->print_values(false);
	cout << "Bailing Out!\n";
	exit(NUMERICAL_ERROR);
#       endif
    }
#   endif
}

int
Energy_exchange_ODE_update::
s_rate_of_change(Gas_data &Q, vector<double> &dedt)
{
    convert_massf2molef( Q.massf, g_->M(), molef_ );
    ees_->set_gas_data_ptr(Q);
    ees_->set_molef_ptr(molef_);
    // NOTES: - starting at itm=1 to skip translational energy
    //        - scaling by modal_massf to convert J/modal-kg to J/total-kg
    for ( size_t itm = 1; itm < Q.e.size(); ++itm )
        yin_[itm-1] = Q.e[itm];
    ees_->called_at_least_once = false;
    ees_->eval(yin_, ydot_);

    // We don't set dedt[0] -- we fill this in from changes in other modes
    dedt[0] = 0.0;    
    for ( size_t itm = 1; itm < dedt.size(); ++itm ) {
	dedt[itm] = ydot_[itm-1];
    }

    return SUCCESS;
}

int
Energy_exchange_ODE_update::
perform_increment(Gas_data &Q, double t_interval, double &dt_suggest)
{
    // SETUP
    ees_->set_gas_data_ptr(Q);
    ees_->set_molef_ptr(molef_);
    // Store e_total so we now how much so put in translational energy later.
    double e_total = accumulate(Q.e.begin(), Q.e.end(), 0.0);
    // NOTES: - starting at itm=1 to skip translational energy
    for ( size_t itm = 1; itm < Q.e.size(); ++itm ) {
        yin_[itm-1] = Q.e[itm];
    }
    ees_->called_at_least_once = false;
    
    double h = dt_suggest;
    bool flag = false;
    if ( h > 0.0 ) {  // then we have a guess for the timestep
	flag = ode_solver_->solve_over_interval(*ees_, 0.0, t_interval, &h,
						yin_, yout_);
	if ( ! flag ) {
	    // then we retry with the timestep selected by our function
	    h = ees_->stepsize_select(yin_);
	    if ( h > (0.5 * dt_suggest) ) {
		// If we're going to reduce the timestep, it's probably
		// best to do so drastically.  Anything less than half
		// what we started with is not really worth it.
		// Let's just try 10% of first timestep size
		h = 0.1 * dt_suggest;
	    }
	    flag = ode_solver_->solve_over_interval(*ees_, 0.0, t_interval, &h,
						    yin_, yout_);
	    if( ! flag ) {
		return NUMERICAL_ERROR;
	    }
	}
    }
    else { // it's probably our first step (or after T_trigger invocation)	
	h = ees_->stepsize_select(yin_);
	flag = ode_solver_->solve_over_interval(*ees_, 0.0, t_interval, &h,
						yin_, yout_);

	if ( ! flag ) {
	    return NUMERICAL_ERROR;
	}
    }

    // 2. If we've made it this far than we're doing well.
    //    Let's update the gas state and leave.
    double e_other = 0.0; // all other energy components except e[0]
    for ( size_t itm = 1; itm < Q.e.size(); ++itm ) {
    	Q.e[itm] = yout_[itm-1];
	e_other += Q.e[itm];
    }
    Q.e[0] = e_total - e_other;

    // 3. But first apply equilibriation mechanisms
    for ( size_t ieq=0; ieq<eq_mechs_.size(); ++ieq ) {
    	eq_mechs_[ieq]->apply( Q );
    }

    dt_suggest = h;

    return SUCCESS;
}

int
Energy_exchange_ODE_update::
estimate_appropriate_subcycle(double t_interval, double dt_suggest,
			      double &dt_sub, int &no_substeps)
{
    const double allowable_t_factor = 1.0;
    const int default_substeps = 50;
    const int max_substeps = 100;

    if ( dt_suggest < 0.0 ) {
	dt_sub = t_interval/default_substeps;
	no_substeps = default_substeps;
    }
    else {
	double target_t = allowable_t_factor * dt_suggest;
	no_substeps = int(t_interval/target_t) + 1;
	if ( no_substeps > max_substeps ) no_substeps = max_substeps;
	dt_sub = t_interval/no_substeps;
    }
    return SUCCESS;
}

Energy_exchange_update* create_Energy_exchange_ODE_update(lua_State *L, Gas_model &g)
{
    return new Energy_exchange_ODE_update(L, g);
}

