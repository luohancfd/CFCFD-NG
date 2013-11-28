// Author: Rowan J. Gollan
// Date: 12-Sep-2008
// Place: NIA, Hampton, Virginia, USA
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <numeric>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "ode_setup.hh"
#include "../models/gas-model.hh"
#include "chemical-kinetic-ODE-update.hh"

using namespace std;

Chemical_kinetic_ODE_update::
Chemical_kinetic_ODE_update(lua_State *L, Gas_model &g)
    : Reaction_update()
{
    lua_getglobal(L, "scheme_t");

    lua_getfield(L, -1, "temperature_limits");
    double T_lower = get_positive_number(L, -1, "lower");
    double T_upper = get_positive_number(L, -1, "upper");
    lua_pop(L, 1);

    double error_tol = get_positive_number(L, -1, "error_tolerance");

    lua_pop(L, 1); // pop scheme_t
  
    if ( T_upper <= T_lower ) {
	ostringstream ost;
	ost << "Chemical_kinetic_ODE_update::Chemical_kinetic_ODE_update():\n";
	ost << "Error in input: T_upper must be greater than T_lower.\n";
	ost << "T_upper= " << T_upper << " T_lower= " << T_lower << endl;
	input_error(ost);
    }
    
    int nsp = g.get_number_of_species();

    lua_getglobal(L, "ode_t");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_kinetic_ODE_update::Chemical_kinetic_ODE_update():\n";
	ost << "Error in input: ode_solver table is missing.\n";
	input_error(ost);
    }
    ode_solver_ = create_ode_solver(L, nsp, "chemical kinetic ODE update system");
    lua_pop(L, 1);

    cks_ = new Chemical_kinetic_system(L, g, error_tol, T_upper, T_lower);

    yin_.resize(nsp, 0.0);
    yout_.resize(nsp, 0.0);
    ydot_.resize(nsp, 0.0);
    M_.resize(nsp, 0.0);
    for ( size_t isp = 0; isp < M_.size(); ++isp ) M_[isp] = g.molecular_weight(isp);
    Q_save_ = new Gas_data(&g);
}

Chemical_kinetic_ODE_update::
~Chemical_kinetic_ODE_update()
{
    delete ode_solver_;
    delete cks_;
    delete Q_save_;
}

int
Chemical_kinetic_ODE_update::
s_update_state(Gas_data &Q, double t_interval, double &dt_suggest, Gas_model *gm)
{
    const double DT_SUB_FRAC = 0.1; // fraction of subinterval to use
                                    // as first dt_suggest
    int flag = SUCCESS;
    int nmodes = gm->get_number_of_modes();

    // Keep a copy in case something goes wrong
    // and we need to retry
    Q_save_->copy_values_from(Q);
    double e_total = accumulate(Q.e.begin(), Q.e.end(), 0.0);
    double dt_suggest_save = dt_suggest;

    // 1. Attempt to solve normally.
    if ( perform_increment(Q, t_interval, dt_suggest) == SUCCESS ) {
	if ( nmodes > 1 ) {
	    // Changing mass fractions does not change internal energies
	    // but the amount in kg/mixture has changed, so update e[] values
	    if ( gm->eval_thermo_state_rhoT(Q) != SUCCESS ) {
		return FAILURE;
	    }
	    double e_other = accumulate(Q.e.begin()+1, Q.e.end(), 0.0);
	    Q.e[0] = e_total - e_other;
	}
	// all is well
	return SUCCESS;
    }
    else { // Repeat the attempt by subcycling.
	Q.copy_values_from(*Q_save_);
	double dt_sub;
	int no_substeps;
	estimate_appropriate_subcycle(t_interval, dt_suggest_save,
				      dt_sub, no_substeps);
	// Now that we've split the t_interval into smaller
	// increments, let's use a smaller timestep within
	// the increments.
	dt_suggest = DT_SUB_FRAC*dt_sub;
	for( int i = 0; i < no_substeps; ++i ) {
	    // Update the gas-state assuming constant density and energy
	    if ( i > 0 && gm!=0 ) {
		if ( nmodes > 1 ) {
		    gm->eval_thermo_state_rhoT(Q);
		    double e_other = accumulate(Q.e.begin()+1, Q.e.end(), 0.0);
		    Q.e[0] = e_total - e_other;
		}
	    	gm->eval_thermo_state_rhoe(Q);
	    }
	    if ( perform_increment(Q, dt_sub, dt_suggest) != SUCCESS ) {
		flag = FAILURE;
		break;
	    }
	}
	if ( flag == SUCCESS ) {
	    // Changing mass fractions does not change internal energies
	    // but the amount in kg/mixture has changed, so update e[] values
	    if ( gm->eval_thermo_state_rhoT(Q) != SUCCESS ) {
		return FAILURE;
	    }
	    double e_other = accumulate(Q.e.begin()+1, Q.e.end(), 0.0);
	    Q.e[0] = e_total - e_other;
	    return SUCCESS;
	}
	else {
	    cout << "Failed to successfully update gas state due to chemistry.\n";
	    cout << "The initial condition was: \n";
	    Q_save_->print_values(false);
	    cout << "Bailing Out!\n";
	    exit(NUMERICAL_ERROR);
	}
    }
}	
	
int
Chemical_kinetic_ODE_update::
s_rate_of_change(Gas_data &Q, vector<double> &dcdt)
{
    cks_->set_gas_data_ptr(Q);
    convert_massf2conc(Q.rho, Q.massf, M_, yin_);
    cks_->called_at_least_once = false;
    cks_->eval(yin_, ydot_);

    for ( size_t isp = 0; isp < ydot_.size(); ++isp ) {
	dcdt[isp] = ydot_[isp];
    }

    return SUCCESS;
}

int
Chemical_kinetic_ODE_update::
s_eval_chemistry_energy_coupling_source_terms( Gas_data &Q, vector<double> &dedt )
{
    // No chemistry-energy coupling model currently applied for this chemistry update scheme
    return 0;
}

int
Chemical_kinetic_ODE_update::
perform_increment(Gas_data &Q, double t_interval, double &dt_suggest)
{
    // SETUP
    cks_->set_gas_data_ptr(Q);
    convert_massf2conc(Q.rho, Q.massf, M_, yin_);
    cks_->called_at_least_once = false;
    
    double h = dt_suggest;
    bool flag = false;

    if ( h > 0.0 ) {  // then we have a guess for the timestep
	flag = ode_solver_->solve_over_interval(*cks_, 0.0, t_interval, &h,
						yin_, yout_);
	if ( ! flag ) {
	    // then we retry with the timestep selected by our function
	    h = cks_->stepsize_select(yin_);
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
	flag = ode_solver_->solve_over_interval(*cks_, 0.0, t_interval, &h,
						yin_, yout_);
	if ( ! flag ) {
	    return NUMERICAL_ERROR;
	}
    }
    // 2. If we've made it this far than we're doing well.
    //    Let's update the gas state and leave.
    
    convert_conc2massf(Q.rho, yout_, M_, Q.massf);

    // 3. Normalise: this should just 'tweak' the values.
    //    If they had gotten really out of whack, then it should
    //    not have passed the ODE system test.  That's why
    //    we don't test how bad things are here; they've already
    //    been tested.
    double massf_sum = 0.0;
    for ( size_t isp = 0; isp < Q.massf.size(); ++isp ) {
	Q.massf[isp] = Q.massf[isp] >= 0.0 ? Q.massf[isp] : 0.0;
	massf_sum += Q.massf[isp];
    }

    for ( size_t isp = 0; isp < Q.massf.size(); ++isp ) {
	Q.massf[isp] /= massf_sum;
    }

    dt_suggest = h;

    return SUCCESS;
}

int
Chemical_kinetic_ODE_update::
estimate_appropriate_subcycle(double t_interval, double dt_suggest,
			      double &dt_sub, int &no_substeps)
{
    const double allowable_t_factor = 2.0;
    const int default_substeps = 50;

    if ( dt_suggest < 0.0 ) {
	dt_sub = t_interval/default_substeps;
	no_substeps = default_substeps;
    }
    else {
	double target_t = allowable_t_factor * dt_suggest;
	no_substeps = int(t_interval/target_t) + 1;
	dt_sub = t_interval/no_substeps;
    }
    return SUCCESS;
}

Reaction_update* create_Chemical_kinetic_ODE_update(lua_State *L, Gas_model &g)
{
    return new Chemical_kinetic_ODE_update(L, g);
}

Chemical_kinetic_system * create_chemical_kinetic_system( std::string cfile, Gas_model &g )
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
    char *e3bin = getenv("E3BIN");
    string home;
    if ( e3bin == NULL ) {
	// Assume default location of $HOME/e3bin
	home.append(getenv("HOME")); home.append("/e3bin");
    }
    else {
	home.append(e3bin);
    }
    string script_file(home);
    script_file.append("/reaction_parser.lua");

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

    Chemical_kinetic_ODE_update * cku = dynamic_cast<Chemical_kinetic_ODE_update*>(create_Chemical_kinetic_ODE_update(L, g));

    if ( cku == 0 ) {
	ostringstream ost;
	ost << "create_chemical_kinetic_system():\n";
	ost << "Error trying to reaction update of type: " << update << endl;
	input_error(ost);
    }
    
    lua_close(L);
    
    // retrive the chemical kinetic system pointer from the Chemical_kinetic_ODE_update

    return cku->get_cks_pointer();
}
