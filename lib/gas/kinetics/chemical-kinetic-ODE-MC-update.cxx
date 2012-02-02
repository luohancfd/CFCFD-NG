// Author: Rowan J. Gollan
// Date: 12-Sep-2008
// Place: NIA, Hampton, Virginia, USA
// History:
//   23-Mar-2009  Revised to accommodate new direct-from-Lua input.
//   07-Apr-2010  Slightly modified to create this 'mass-conserved' version

#include <iostream>
#include <sstream>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "ode_setup.hh"
#include "../models/gas-model.hh"
#include "chemical-kinetic-ODE-MC-update.hh"

#define PERSISTENT_CHEM_UPDATE 0
#define ALLOW_CHEM_SUBCYCLING 1

using namespace std;

Chemical_kinetic_ODE_MC_update::
Chemical_kinetic_ODE_MC_update(lua_State *L, Gas_model &g)
    : Reaction_update()
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
	ost << "Chemical_kinetic_ODE_MC_update::Chemical_kinetic_ODE_MC_update():\n";
	ost << "Error in input: T_upper_limit must be greater than T_lower_limit.\n";
	ost << "T_upper_limit= " << T_upper_limit_ << " T_lower_limit= " << T_lower_limit_ << endl;
	input_error(ost);
    }
    
    // need to determine the number of reactions to initialise the cks
    lua_getglobal(L, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_kinetic_ODE_MC_update::Chemical_kinetic_ODE_MC_update()\n";
	ost << "Error interpreting 'reactions'; a table of reactions is expected.\n";
	input_error(ost);
    }
    int nreac = lua_objlen(L, -1);
    lua_pop(L, 1); // pop reactions
    
    cks_ = new Chemical_kinetic_MC_system(L, g, nreac, error_tol);

    lua_getglobal(L, "ode_t");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_kinetic_ODE_MC_update::Chemical_kinetic_ODE_MC_update():\n";
	ost << "Error in input: ode_solver table is missing.\n";
	input_error(ost);
    }
    ode_solver_ = create_ode_solver(L, 2*nreac, "chemical kinetic ODE MC update system");
    lua_pop(L, 1);
    
    yin_.resize(2*nreac, 0.0);
    yout_.resize(2*nreac, 0.0);
    ydot_.resize(2*nreac, 0.0);

    int nsp = g.get_number_of_species();
    
    cdot_.resize(nsp, 0.0);
    c_.resize(nsp, 0.0);
    M_.resize(nsp, 0.0);
    for ( size_t isp = 0; isp < M_.size(); ++isp ) M_[isp] = g.molecular_weight(isp);
    Q_save_ = new Gas_data(&g);
}

Chemical_kinetic_ODE_MC_update::
~Chemical_kinetic_ODE_MC_update()
{
    delete ode_solver_;
    delete cks_;
    delete Q_save_;
}

int
Chemical_kinetic_ODE_MC_update::
s_update_state(Gas_data &Q, double t_interval, double &dt_suggest, Gas_model *gm)
{
#   if ALLOW_CHEM_SUBCYCLING
    int flag = SUCCESS;
#   endif

    if ( Q.T[0] <= T_lower_limit_  || Q.T[0] >= T_upper_limit_ ) {
	dt_suggest = -1.0;
	return SUCCESS;
    }

    // Keep a copy in case something goes wrong
    // and we need to retry
    Q_save_->copy_values_from(Q);
#   if ALLOW_CHEM_SUBCYCLING
    double dt_suggest_save = dt_suggest;
#   endif

    // 1. Attempt to solve normally.
    if ( perform_increment(Q, t_interval, dt_suggest) == SUCCESS ) {
	// all is well
	return SUCCESS;
    }
#   if ALLOW_CHEM_SUBCYCLING == 1
    else { // Repeat the attempt by subcycling.
	Q.copy_values_from(*Q_save_);
	double dt_sub;
	int no_substeps;
	estimate_appropriate_subcycle(t_interval, dt_suggest_save,
				      dt_sub, no_substeps);
	cout << "Chemical_kinetic_ODE_MC_update::s_update_state()" << endl
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
	    cout << "Failed to successfully update gas state due to chemistry.\n";
	    cks_->print_limiting_species_and_reaction();
#           if PERSISTENT_CHEM_UPDATE == 1
            cout << "Freezing the gas state and attempting to continue.\n";
            Q.copy_values_from(*Q_save_);
            return SUCCESS;
#           else
	    cout << "The initial condition was: \n";
	    Q_save_->print_values(false);
	    cout << "Bailing Out!\n";
	    exit(NUMERICAL_ERROR);
#           endif
	}
    }
#   else
    else {
	cout << "Failed to successfully update gas state due to chemistry.\n";
	cks_->print_limiting_species_and_reaction();
#       if PERSISTENT_CHEM_UPDATE == 1
	cout << "Freezing the gas state and attempting to continue (dt_suggest = " << dt_suggest << ").\n";
	Q.copy_values_from(Q_save_);
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
Chemical_kinetic_ODE_MC_update::
s_rate_of_change(Gas_data &Q, vector<double> &dcdt)
{
    for ( size_t ir=0; ir<yin_.size(); ++ir ) yin_[ir] = 0.0;
    convert_massf2conc(Q.rho, Q.massf, M_, c_);
    cks_->set_gas_data_ptr_and_initial_concs(Q, c_);
    cks_->called_at_least_once = false;
    
    cks_->eval_species_rates(yin_, cdot_);

    for ( size_t isp = 0; isp < cdot_.size(); ++isp ) {
	dcdt[isp] = cdot_[isp];
    }

    return SUCCESS;
}

int
Chemical_kinetic_ODE_MC_update::
s_eval_chemistry_energy_coupling_source_terms( Gas_data &Q, vector<double> &dedt )
{
    // Quick exit if no mechanisms are present
    if ( cks_->cecs_size()==0 ) return 0;

    // 0. clear the yin_ array
    for ( size_t ir=0; ir<yin_.size(); ++ir ) yin_[ir] = 0.0;
    
    // 1. Set the e_old_ values in the coupling components
    // NOTE: c_ should just be a zero vector, we won't be using N_old_
    cks_->initialise_chemistry_energy_coupling(Q, c_);
    
    // 2. Evaluate the source terms
    return cks_->eval_chemistry_energy_coupling_source_terms( Q, yin_, dedt );
}

int
Chemical_kinetic_ODE_MC_update::
perform_increment(Gas_data &Q, double t_interval, double &dt_suggest )
{
    // SETUP
    for ( size_t ir=0; ir<yin_.size(); ++ir ) yin_[ir] = 0.0;
    convert_massf2conc(Q.rho, Q.massf, M_, c_);
    cks_->set_gas_data_ptr_and_initial_concs(Q, c_);
    cks_->initialise_chemistry_energy_coupling(Q, c_);
    cks_->called_at_least_once = false;
    
    double h = dt_suggest;
    bool flag = false;

    if ( h > 0.0 ) {  // then we have a guess for the timestep
	flag = ode_solver_->solve_over_interval(*cks_, 0.0, t_interval, &h,
						yin_, yout_);
	if ( ! flag ) {
	    //printf("step failed! retrying...");
	    // then we retry with the timestep selected by our function
	    h = cks_->stepsize_select(yin_);
	    //printf("subsequent step h=%6.5e\n", h);
	    if ( h > (0.5 * dt_suggest) ) {
		// If we're going to reduce the timestep, it's probably
		// best to do so drastically.  Anything less than half
		// what we started with is not really worth it.
		// Let's just try 10% of first timestep size
		h = 0.1 * dt_suggest;
	    }
	    flag = ode_solver_->solve_over_interval(*cks_, 0.0, t_interval, &h,
						    yin_, yout_);
	    if( ! flag ) {
		return NUMERICAL_ERROR;
	    }
	}
    }
    else { // it's probably our first step (or after T_trigger invocation)	
	h = cks_->stepsize_select(yin_);
	//printf("first step h=%6.5e\n", h);
	flag = ode_solver_->solve_over_interval(*cks_, 0.0, t_interval, &h,
						yin_, yout_);

	//printf("yout_ delta_moles = [");
	//for (size_t ir = 0; ir < yout_.size(); ++ir) {
	//    printf("%4.3e, ", yout_[ir]);
	//}
	//printf("]\n");

	if ( ! flag ) {
	    return NUMERICAL_ERROR;
	}
    }
    // 2. If we've made it this far than we're doing well.
    //    Let's assemble the solution, update the gas state and leave.
    
    cks_->eval_new_concentrations( yout_, c_ );
    
    convert_conc2massf(Q.rho, c_, M_, Q.massf);
    
    // 4. Apply chemistry-energy coupling model(s)
    if ( cks_->apply_chemistry_energy_coupling( Q, yout_, c_ ) ) {
    	cout << "Chemical_kinetic_ODE_MC_update::perform_increment()" << endl
    	     << "Failed when applying chemistry-energy coupling terms." << endl
    	     << "Exiting." << endl;
    	exit( FAILURE );
    }

    dt_suggest = h;

    return SUCCESS;
}

int
Chemical_kinetic_ODE_MC_update::
estimate_appropriate_subcycle(double t_interval, double dt_suggest,
			      double &dt_sub, int &no_substeps)
{
    const double allowable_t_factor = 2.0;
    const int default_substeps = 50;
    const int max_substeps = 1000;

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

Reaction_update* create_Chemical_kinetic_ODE_MC_update(lua_State *L, Gas_model &g)
{
    return new Chemical_kinetic_ODE_MC_update(L, g);
}

Chemical_kinetic_MC_system * create_chemical_kinetic_MC_system( std::string cfile, Gas_model &g )
{
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    // Set up species table
    lua_newtable(L);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	// This species table maps to C++ indices, because
	// it is used to setup the integer maps for
	// the reaction coefficients.
	lua_pushinteger(L, isp);
	lua_setfield(L, -2, g.species_name(isp).c_str());
    }
    // Plus add a field 'size': no of species
    lua_pushinteger(L, g.get_number_of_species());
    lua_setfield(L, -2, "size");
    lua_setglobal(L, "species");

    // Path to reaction parsing script
    string home(getenv("HOME"));
    string script_file(home);
    script_file.append("/e3bin/reaction_parser.lua");

    if ( luaL_dofile(L, script_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error in loading script file: " << script_file << endl;
	input_error(ost);
    }
    
    // Parse the input file...
    lua_getglobal(L, "main");
    lua_pushstring(L, cfile.c_str());
    if ( lua_pcall(L, 1, 0, 0) != 0 ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error trying to load reaction scheme file: " << cfile << endl;
	ost << "Lua error message: " << lua_tostring(L, -1) << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "scheme_t");
    if ( lua_isnil(L, -1) ) {
	cout << "scheme_t not set globally.";
	exit(1);
    }

    lua_getfield(L, -1, "update");
    string update(luaL_checkstring(L, -1));
    lua_pop(L, 1);

    Chemical_kinetic_ODE_MC_update * cku = dynamic_cast<Chemical_kinetic_ODE_MC_update*>(create_Chemical_kinetic_ODE_MC_update(L, g));

    if ( cku == 0 ) {
	ostringstream ost;
	ost << "create_chemical_kinetic_MC_system():\n";
	ost << "Error trying to reaction update of type: " << update << endl;
	input_error(ost);
    }
    
    lua_close(L);
    
    // retrive the chemical kinetic system pointer from the Chemical_kinetic_ODE_MC_update

    return cku->get_cks_pointer();
}
