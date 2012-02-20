// Author: Rowan J. Gollan
// Date: 03-Jul-2008
//
// History:
//   20-Feb-2012
//   Changed all indexing of Lua arrays to begin at 0
//   for consistency for use throughout eilmer3 app code.
//


#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "lservice_gas_data.hh"
#include "user-defined-gas-model.hh"

using namespace std;

int call_user_function(lua_State *L, const char *func, Gas_data &Q)
{
    lua_getglobal(L, func);
    push_gas_data_as_table(L, Q);
    // 1 arg, 1 return value: gas data as table
    int number_args = 1;
    int number_results = 1;
    if( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "Error running user function %s:\n%s\n",
			 func, lua_tostring(L, -1));
    }
    get_table_as_gas_data(L, Q);
    // Clear stack
    lua_settop(L, 0);
    return SUCCESS;
}

double call_user_deriv_function(lua_State *L, const char *func, const Gas_data &Q)
{
    lua_getglobal(L, func);
    push_gas_data_as_table(L, Q);
    // 1 arg, 1 return value: gas data as table
    int number_args = 1;
    int number_results = 1;
    if( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L, "Error running user function %s:\n%s\n",
			 func, lua_tostring(L, -1));
    }
    
    double deriv = luaL_checknumber(L, -1);
    // Clear stack
    lua_settop(L, 0);
    return deriv;
}

UD_gas_model::
UD_gas_model(string cfile)
    : Gas_model()
{
    //    cout << "User-defined gas model NOT presently available in elmer3.\n";
    // cout << "Bailing out!\n";
    //exit(1);

    L_ = luaL_newstate();
    luaL_openlibs(L_);

    if( luaL_dofile(L_, cfile.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Error in user-defined gas model file: " << cfile << endl;
	input_error(ost);
    }
    
    // Now check that the user has specified the minimum functions and values.
    int nsp = check_for_integer(L_, cfile, "nsp", "number of species");
    set_number_of_species(nsp);

    int nmodes = check_for_integer(L_, cfile, "nmodes", "number of thermal modes");
    set_number_of_modes(nmodes);

    s_names_.resize(nsp);
    // Set some default names.
    for ( int isp = 0; isp < nsp; ++isp ) {
	ostringstream ost;
	ost << "sp-" << isp;
	s_names_[isp] = ost.str();
    }

    lua_getglobal(L_, "species");
    if ( lua_istable(L_, -1) ) {
	int n = lua_objlen(L_, -1);
	if ( n != nsp ) {
	    cout << "The number of entries in the 'species' table: " << n << endl;
	    cout << "does not match the number of species specified: " << nsp << endl;
	    cout << "Check the following user-defined gas model file for the error: " << cfile << endl;
	    cout << "Bailing out!\n";
	    exit(1);
	}
	
	for ( int isp = 0; isp < nsp; ++isp ) {
	    lua_rawgeti(L_, -1, isp);
	    s_names_[isp] = luaL_checkstring(L_, -1);
	    lua_pop(L_, 1);
	}
    }
    lua_pop(L_, 1);

    if ( !check_for_function(L_, "eval_thermo_state_rhoe") ) {
	ostringstream ost;
	ost << "Error in user-defined gas model file: " << cfile << endl;
	ost << "A user-specified function 'eval_thermo_state_rhoe' \n";
	ost << "must be specified at a minimum.\n";
	input_error(ost);
    }

    if ( !check_for_function(L_, "eval_transport_coefficients") ) {
	ostringstream ost;
	ost << "Error in user-defined gas model file: " << cfile << endl;
	ost << "A user-specified function 'eval_transport_coefficients' \n";
	ost << "must be specified at a minimum.\n";
	input_error(ost);
    }

    if ( !check_for_function(L_, "molecular_weight") ) {
	ostringstream ost;
	ost << "Error in user-defined gas model file: " << cfile << endl;
	ost << "A user-specified function 'molecular_weight' \n";
	ost << "must be specified.\n";
	input_error(ost);
    }

    // Set the molecular weights: M_
    M_.resize(nsp);
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_getglobal(L_, "molecular_weight");
	lua_pushinteger(L_, isp);
	int number_args = 1;
	int number_results = 1;
	if ( lua_pcall(L_, number_args, number_results, 0) != 0 ) {
	    handle_lua_error(L_, "Error running user function molecular_weight:\n%s\n",
			     lua_tostring(L_, -1));
	}
	M_[isp] = luaL_checknumber(L_, -1);
	// Clear stack
	lua_settop(L_, 0);
    }

}

UD_gas_model::
~UD_gas_model()
{
    lua_close(L_);
}

int
UD_gas_model::
s_decode_conserved_energy(Gas_data &Q, const vector<double> &rhoe)
{
    if ( !check_for_function(L_, "decode_conserved_energy") )
	return Gas_model::s_decode_conserved_energy(Q, rhoe);
    // else go on with calling user's function...
    lua_getglobal(L_, "decode_conserved_energy");
    push_gas_data_as_table(L_, Q);
    push_vector_as_table(L_, rhoe);
    // 2 args, 1 return value: gas data as table
    int number_args = 2;
    int number_results = 1;
    lua_pcall(L_, number_args, number_results, 0);
    // Return a table of gas_data...
    lua_getfield(L_, -1, "e");
    get_table_as_vector(L_, Q.e);
    return SUCCESS;
}


int
UD_gas_model::
s_encode_conserved_energy(const Gas_data &Q, vector<double> &rhoe)
{
    if ( !check_for_function(L_, "encode_conserved_energy") )
	return Gas_model::s_encode_conserved_energy(Q, rhoe);
    // else go on with calling user's function...
    lua_getglobal(L_, "encode_conserved_energy");
    push_gas_data_as_table(L_, Q);
    // 2 args, 1 return value: rhoe as array
    int number_args = 1;
    int number_results = 1;
    lua_pcall(L_, number_args, number_results, 0);
    // Return an array of rhoe values
    get_table_as_vector(L_, rhoe);
    return SUCCESS;
}

int
UD_gas_model::
s_eval_thermo_state_rhoe(Gas_data &Q)
{
    return call_user_function(L_, "eval_thermo_state_rhoe", Q);
}


int
UD_gas_model::
s_eval_thermo_state_pT(Gas_data &Q)
{
    if ( !check_for_function(L_, "eval_thermo_state_pT") )
	return Gas_model::s_eval_thermo_state_pT(Q);
    // else proceed to call user-defined function...
    return call_user_function(L_, "eval_thermo_state_pT", Q);
}

int
UD_gas_model::
s_eval_thermo_state_rhoT(Gas_data &Q)
{
    if ( !check_for_function(L_, "eval_thermo_state_rhoT") )
	return Gas_model::s_eval_thermo_state_rhoT(Q);
    // else proceed to call user-defined function...
    return call_user_function(L_, "eval_thermo_state_rhoT", Q);
}

int
UD_gas_model::
s_eval_thermo_state_rhop(Gas_data &Q)
{
    if ( !check_for_function(L_, "eval_thermo_state_rhop") )
	return Gas_model::s_eval_thermo_state_rhop(Q);
    // else proceed to call user-defined function...
    return call_user_function(L_, "eval_thermo_state_rhop", Q);
}

int
UD_gas_model::
s_eval_transport_coefficients(Gas_data &Q)
{
    return call_user_function(L_, "eval_transport_coefficients", Q);
}

int
UD_gas_model::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    return call_user_function(L_, "eval_diffusion_coefficients", Q);
}

double
UD_gas_model::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    if ( !check_for_function(L_, "dTdp_const_rho") )
	return Gas_model::s_dTdp_const_rho(Q, status);
    // else call users function.
    status = SUCCESS;
    return call_user_deriv_function(L_, "dTdp_const_rho", Q);
}


double
UD_gas_model::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    if ( !check_for_function(L_, "dTdrho_const_p") )
	return Gas_model::s_dTdrho_const_p(Q, status);
    // else call users function.
    status = SUCCESS;
    return call_user_deriv_function(L_, "dTdrho_const_p", Q);
}

double
UD_gas_model::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    if ( !check_for_function(L_, "dpdrho_const_T") )
	return Gas_model::s_dpdrho_const_T(Q, status);
    // else call users function.
    status = SUCCESS;
    return call_user_deriv_function(L_, "dpdrho_const_T", Q);
}

double
UD_gas_model::
s_dedT_const_v(const Gas_data &Q, int &status)
{
    if ( !check_for_function(L_, "dedT_const_v") )
	return Gas_model::s_dedT_const_v(Q, status);
    // else call users function.
    status = SUCCESS;
    return call_user_deriv_function(L_, "dedT_const_v", Q);
}

double
UD_gas_model::
s_dhdT_const_p(const Gas_data &Q, int &status)
{
    if ( !check_for_function(L_, "dhdT_const_p") )
	return Gas_model::s_dhdT_const_p(Q, status);
    // else call users function.
    status = SUCCESS;
    return call_user_deriv_function(L_, "dhdT_const_p", Q);
}

double
UD_gas_model::
s_gas_constant(const Gas_data &Q, int &status)
{
    if ( !check_for_function(L_, "gas_constant") )
	return Gas_model::s_gas_constant(Q, status);
    // else call users function.
    status = SUCCESS;
    // Not really a deriv function, but the mechanics
    // of calling it are the same.
    return call_user_deriv_function(L_, "gas_constant", Q);
}

double
UD_gas_model::
s_molecular_weight(int isp)
{
    lua_getglobal(L_, "molecular_weight");
    lua_pushinteger(L_, isp);
    int number_args = 1;
    int number_results = 1;
    if ( lua_pcall(L_, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L_, "Error running user function molecular_weight:\n%s\n",
			 lua_tostring(L_, -1));
    }
    double mw = luaL_checknumber(L_, -1);
    // Clear stack
    lua_settop(L_, 0);
    return mw;
}

double
UD_gas_model::
s_internal_energy(const Gas_data &Q, int isp)
{
    if ( !check_for_function(L_, "energy") )
	return 0.0;
    // else call users function.
    int number_args = 2;
    int number_results = 1;
    lua_getglobal(L_, "energy");
    push_gas_data_as_table(L_, Q);
    lua_pushinteger(L_, isp);
    if ( lua_pcall(L_, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L_, "Error running user function energy:\n%s\n",
			 lua_tostring(L_, -1));
    }
    double e = luaL_checknumber(L_, -1);
    // Clear stack
    lua_settop(L_, 0);
    return e;
}

double
UD_gas_model::
s_enthalpy(const Gas_data &Q, int isp)
{
    // This is only necessary if doing chemistry AND
    // the equilibrium constant is required.
    if ( !check_for_function(L_, "enthalpy") )
	return 0.0;
    // else call users function.
    int number_args = 2;
    int number_results = 1;
    lua_getglobal(L_, "enthalpy");
    push_gas_data_as_table(L_, Q);
    lua_pushinteger(L_, isp);
    if ( lua_pcall(L_, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L_, "Error running user function enthalpy:\n%s\n",
			 lua_tostring(L_, -1));
    }
    double h = luaL_checknumber(L_, -1);
    // Clear stack
    lua_settop(L_, 0);
    return h;
}

double
UD_gas_model::
s_entropy(const Gas_data &Q, int isp)
{
    // This is only necessary if doing chemistry AND
    // the equilibrium constant is required.
    if ( !check_for_function(L_, "entropy") )
	return 0.0;
    // else call users function.
    int number_args = 2;
    int number_results = 1;
    lua_getglobal(L_, "entropy");
    push_gas_data_as_table(L_, Q);
    lua_pushinteger(L_, isp);
    if ( lua_pcall(L_, number_args, number_results, 0) != 0 ) {
	handle_lua_error(L_, "Error running user function entropy:\n%s\n",
			 lua_tostring(L_, -1));
    }
    double s = luaL_checknumber(L_, -1);
    // Clear stack
    lua_settop(L_, 0);
    return s;
}

Gas_model* create_user_defined_gas_model(string cfile)
{
    return new UD_gas_model(cfile);
}

