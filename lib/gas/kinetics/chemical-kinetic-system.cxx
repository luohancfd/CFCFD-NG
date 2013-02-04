// Author: Rowan J. Gollan
// Date: 12-Sep-2008

#include <iostream>
#include <sstream>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"
#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "chemical-kinetic-system.hh"

using namespace std;

Chemical_kinetic_system::
Chemical_kinetic_system(lua_State *L, Gas_model &g, double error_tol)
    : OdeSystem(g.get_number_of_species(), true), err_tol_(error_tol)
{
    lua_getglobal(L, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_kinetic_system::Chemical_kinetic_system()\n";
	ost << "Error interpreting 'reactions'; a table of reactions is expected.\n";
	input_error(ost);
    }

    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i);
	reaction_.push_back(create_Reaction(L, g));
	lua_pop(L, 1);
    }
    lua_pop(L, 1);

    int nsp = g.get_number_of_species();
    participation_.resize(nsp);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
	    if ( reaction_[ir]->get_nu(isp) != 0 ) {
		participation_[isp].push_back(ir);
	    }
	}
    }

    q_.resize(nsp, 0.0);
    L_.resize(nsp, 0.0);
    ydot_.resize(nsp, 0.0);
    massf_.resize(nsp, 0.0);
    M_.resize(nsp, 0.0);
    for ( size_t isp = 0; isp < M_.size(); ++isp ) M_[isp] = g.molecular_weight(isp);
}

Chemical_kinetic_system::
Chemical_kinetic_system(string cfile, Gas_model &g, double error_tol)
    : OdeSystem(g.get_number_of_species(), true), err_tol_(error_tol)
{
    // Do some pre-work on the cfile to massage
    // it into a state to be parsed for the
    // individual reactions.

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
	ost << "Chemical_kinetic_system():\n";
	ost << "Error in loading script file: " << script_file << endl;
	input_error(ost);
    }
    
    // Parse the input file...
    lua_getglobal(L, "main");
    lua_pushstring(L, cfile.c_str());
    if ( lua_pcall(L, 1, 0, 0) != 0 ) {
	ostringstream ost;
	ost << "Chemical_kinetic_system():\n";
	ost << "Error trying to load reaction scheme file: " << cfile << endl;
	ost << "Lua error message: " << lua_tostring(L, -1) << endl;
	input_error(ost);
    }
    

    lua_getglobal(L, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_kinetic_system::Chemical_kinetic_system()\n";
	ost << "Error interpreting 'reactions'; a table of reactions is expected.\n";
	input_error(ost);
    }

    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i);
	reaction_.push_back(create_Reaction(L, g));
	lua_pop(L, 1);
    }
    lua_pop(L, 1);

    int nsp = g.get_number_of_species();
    participation_.resize(nsp);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
	    if ( reaction_[ir]->get_nu(isp) != 0 ) {
		participation_[isp].push_back(ir);
	    }
	}
    }

    q_.resize(nsp, 0.0);
    L_.resize(nsp, 0.0);
    ydot_.resize(nsp, 0.0);
    massf_.resize(nsp, 0.0);
    M_.resize(nsp, 0.0);
    for ( size_t isp = 0; isp < M_.size(); ++isp ) M_[isp] = g.molecular_weight(isp);
    
    lua_close(L);
}

Chemical_kinetic_system::
~Chemical_kinetic_system()
{
    for ( size_t i = 0; i < reaction_.size(); ++i ) {
	delete reaction_[i];
    }
}

int
Chemical_kinetic_system::
eval(const valarray<double> &y, valarray<double> &ydot)
{
    eval_split(y, q_, L_);
    //printf("ydot=[");
    for ( size_t isp = 0; isp < ydot.size(); ++isp ) {
	ydot[isp] = q_[isp] - L_[isp];
	//printf("%4.3e, ", ydot[isp]);
    }
    //printf("]\n");
    return SUCCESS;
}

int
Chemical_kinetic_system::
eval_split(const valarray<double> &y,
	   valarray<double> &q, valarray<double> &L)
{
    if ( ! called_at_least_once ) {
	for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
	    if ( reaction_[ir]->compute_kf_first() ) {
		reaction_[ir]->compute_k_f(*Q_);
		reaction_[ir]->compute_k_b(*Q_);
	    }
	    else {
		reaction_[ir]->compute_k_b(*Q_);
		reaction_[ir]->compute_k_f(*Q_);
	    }
	    //printf("rxn[%i]: kf=%16.15e, kb=%16.15e\n", ir, reaction_[ir]->k_f(), reaction_[ir]->k_b());
	}
	called_at_least_once = true;
    }
    //printf("conc = [");
    //for (size_t sp = 0; sp < y.size(); ++sp)
    //	printf("%4.3e, ", y[sp]);
    //printf("]\n");

    for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
	reaction_[ir]->compute_forward_rate(y);
	reaction_[ir]->compute_backward_rate(y);
	//printf("rxn[%i]: w_f_=%16.15e, w_b_=%16.15e\n", ir, reaction_[ir]->w_f(), reaction_[ir]->w_b());
    }

    for ( size_t isp = 0; isp < participation_.size(); ++isp ) {
	q[isp] = 0.0;
	L[isp] = 0.0;
	for ( size_t j = 0; j < participation_[isp].size(); ++j ) {
	    int ir = participation_[isp][j];
	    q[isp] += reaction_[ir]->production(isp);
	    L[isp] += reaction_[ir]->loss(isp);
	}
    }
    return SUCCESS;
}

const double eps1 = 0.001;
const double chem_step_upper_limit = 1.0e-3;
const double chem_step_lower_limit = 1.0e-20;
const double zero_tol = 1.0e-30;

double
Chemical_kinetic_system::
stepsize_select(const valarray<double> &y)
{
    eval(y, ydot_);
    double min_dt = chem_step_upper_limit; // to get us started
    double old_dt = 0.0;

    for( int isp = 0; isp < ndim_; ++isp ) {
	if( (y[isp] > 0.0) && (fabs(ydot_[isp]) > zero_tol) ) {
	    old_dt = fabs( y[isp] / ydot_[isp] );
	    if( old_dt < min_dt ) {
		min_dt = old_dt;
	    }
	}
    }

    double dt_chem = eps1 * min_dt;

    // Impose upper and lower chem_step limits
    if( dt_chem > chem_step_upper_limit )
	dt_chem = chem_step_upper_limit;

    if ( dt_chem < chem_step_lower_limit )
	dt_chem = chem_step_lower_limit;

    return dt_chem;
}

const double min_conc = 1.0e-30;

bool
Chemical_kinetic_system::
passes_system_test(valarray<double> &y)
{
    for ( size_t isp = 0; isp < y.size(); ++isp ) {
	y[isp] = y[isp] < min_conc ? 0.0 : y[isp];
    }

    convert_conc2massf(Q_->rho, y, M_, massf_);

    double massf_sum = 0.0;
    for ( size_t isp = 0; isp < massf_.size(); ++isp ) {
	massf_sum += massf_[isp];
    }

    double lower_lim = 1.0 - ( err_tol_ / 100.0 );
    double upper_lim = 1.0 + ( err_tol_ / 100.0 );

    if( (massf_sum < lower_lim) || (massf_sum > upper_lim) )
	return false;
    else {
	return true;
    }
}

int
Chemical_kinetic_system::
get_directional_rates( vector<double> &w_f, vector<double> &w_b )
{
    for ( size_t ir = 0; ir < reaction_.size(); ++ir ) {
        w_f.push_back( reaction_[ir]->w_f() );
        w_b.push_back( reaction_[ir]->w_b() );
    }

    return SUCCESS;
}
