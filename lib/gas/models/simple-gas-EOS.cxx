// Author: Kan Qin
// Based on perfect gas EOS of Rowan Gollan
// Version: 29-Apirl-2014

#include <iostream>
#include <sstream>
#include <cstdlib>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "simple-gas-EOS.hh"

using namespace std;

Simple_gas::
Simple_gas(lua_State *L)
    : Equation_of_state()
{
    lua_getglobal(L, "species");
    int nsp = lua_objlen(L, -1);
    if ( nsp > 1 ) {
	cout << "ERROR: The simple gas model for very dense gases is only\n";
	cout << "implemented for a single species. However, the number of\n";
	cout << "species was set to: " << nsp << endl;
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Simple_gas::Simple_gas():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    // initialise electron index to -1
    iel_ = -1;

    lua_rawgeti(L, -1, 1); // A Lua list is offset one from the C++ vector index
    const char* sp = luaL_checkstring(L, -1);
    lua_pop(L, 1);

    // Now bring the specific species table to TOS
    lua_getglobal(L, sp);
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Simple_gas::Simple_gas()\n";
	ost << "Error locating information table for species: " << sp << endl;
	input_error(ost);
    }
    double M;
    M = get_positive_value(L, -1, "M");
    M_.push_back(M);
    R_ = PC_R_u/M;
    p_a_ = get_positive_value(L, -1, "p_a");
    rho_a_ = get_positive_value(L, -1, "rho_a");
    k_s_ = get_positive_value(L, -1, "k_s");
        
    // check if this is an electron
    if ( string(sp)=="e_minus" ) iel_ = 0;
    
    lua_pop(L, 1); // pop "sp" off stack
    lua_pop(L, 1); // pop "species" off stack
}

Simple_gas::
~Simple_gas() {}

int
Simple_gas::
s_eval_pressure(Gas_data &Q)
{
    Q.p = simple_pressure(Q.rho, p_a_, rho_a_, k_s_);
    return SUCCESS;
}

int
Simple_gas::
s_eval_temperature(Gas_data &Q)
{
    Q.T[0] = simple_temperature(Q.rho, Q.p, R_);
    return SUCCESS;
}

int
Simple_gas::
s_eval_density(Gas_data &Q)
{
    Q.rho = simple_density(Q.p, p_a_, rho_a_, k_s_);
    return SUCCESS;
}

double
Simple_gas::
s_gas_constant(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return calculate_gas_constant(Q.massf, M_);
}

double 
Simple_gas::
s_prho_ratio(const Gas_data &Q, int isp)
{
    return R_*Q.T[0];
}

double
Simple_gas::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    return 1.0/(Q.rho*R_);
}

double
Simple_gas::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    return (-1.0*Q.p)/(R_*Q.rho*Q.rho);
}

double
Simple_gas::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    return R_*Q.T[0];
}

double
Simple_gas::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    return PC_R_u/M_[isp]*Q.T[0];
}

double
Simple_gas::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    return Q.rho*R_;
}

double simple_pressure(double rho, double p_a, double rho_a, double k_s)
{
    return rho*k_s*p_a/rho_a;
}

double simple_temperature(double rho, double p, double R)
{
    return p/(rho*R);
}

double simple_density(double p, double p_a, double rho_a, double k_s)
{
    return p/p_a*rho_a*k_s;
}
