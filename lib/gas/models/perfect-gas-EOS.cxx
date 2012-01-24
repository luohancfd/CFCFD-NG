// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding.
//

#include <iostream>
#include <sstream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "perfect-gas-EOS.hh"

using namespace std;

Perfect_gas::
Perfect_gas(lua_State *L)
    : Equation_of_state()
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Perfect_gas::Perfect_gas():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    // initialise electron index to -1
    iel_ = -1;
    
    int nsp = lua_objlen(L, -1);
    double M;
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Perfect_gas::Perfect_gas()\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}
	M = get_positive_value(L, -1, "M");
	M_.push_back(M);
	R_.push_back(PC_R_u/M);
	
	// check if this is an electron
	if ( string(sp)=="e_minus" ) iel_ = isp;

	lua_pop(L, 1); // pop "sp" off stack
    }
    lua_pop(L, 1); // pop "species" off stack
}

Perfect_gas::
~Perfect_gas() {}

int
Perfect_gas::
s_eval_pressure(Gas_data &Q)
{
    int status;
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.p = pg_pressure(Q.rho, Q.T[0], R);
    
    // fill in the electron pressure if electrons exist
    if ( iel_ > -1 ) {
    	double R_e = Q.massf[iel_] * R_[iel_];
    	Q.p_e = Q.rho * Q.T[0] * R_e;
    }
    
    return SUCCESS;
}

int
Perfect_gas::
s_eval_temperature(Gas_data &Q)
{
    int status;
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.T[0] = pg_temperature(Q.rho, Q.p, R);
    return SUCCESS;
}

int
Perfect_gas::
s_eval_density(Gas_data &Q)
{
    int status;
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.rho = pg_density(Q.T[0], Q.p, R);
    return SUCCESS;
}

double
Perfect_gas::
s_gas_constant(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return calculate_gas_constant(Q.massf, M_);
}

double 
Perfect_gas::
s_prho_ratio(const Gas_data &Q, int isp)
{
    return R_[isp]*Q.T[0];
}

double
Perfect_gas::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    return 1.0/(Q.rho*R);
}

double
Perfect_gas::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    return (-1.0*Q.p)/(R*Q.rho*Q.rho);
}

double
Perfect_gas::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    return R*Q.T[0];
}

double
Perfect_gas::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    return PC_R_u/M_[isp]*Q.T[0];
}

double
Perfect_gas::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    double R = s_gas_constant(Q, status);
    return Q.rho*R;
}

double pg_pressure(double rho, double T, double R)
{
    return rho*R*T;
}

double pg_temperature(double rho, double p, double R)
{
    return p/(rho*R);
}

double pg_density(double T, double p, double R)
{
    return p/(R*T);
}
