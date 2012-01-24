// Author: Brendan T. O'Flaherty
//  based on perfect-gas-EOS code by Rowan J. Gollan
// Version: 24-May-2008: Initial coding.
//          08-Feb-2011: PJ added pieces needed to recreate gas models for casbar.
//

#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "noble-abel-gas-EOS.hh"

using namespace std;

Noble_Abel_gas::
Noble_Abel_gas(lua_State *L)
    : Equation_of_state()
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Noble_Abel_gas::Noble_Abel_gas():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, -1);
    
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Noble_Abel_gas::Noble_Abel_gas()\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	double M = get_positive_value(L, -1, "M"); // molecular mass kg/mole
        M_.push_back(M);
        R_.push_back(PC_R_u/M);

	// Allow direct setting of the covolume, dropping back to using
	// the critical temperature and pressure if covolume is not found.
	lua_getfield(L, -1, "b");
	int covolume_is_present = !lua_isnil(L, -1);
	lua_pop(L, 1); // discard b because get_value will get it again.
	
	double nu_0 = 0.0;
	if ( covolume_is_present ) {
	    nu_0 = get_value(L, -1, "b");
	    printf("Noble-Abel gas model: setting covolume directly to %g\n", nu_0);
	} else {
	    // Brendan has a reference for this connection between parameters.
	    double T_c = get_value(L, -1, "T_c");
	    double p_c = get_positive_value(L, -1, "p_c");
	    nu_0 = 0.125*PC_R_u*T_c/p_c;
	}
	nu_0_.push_back(nu_0);
	lua_pop(L, 1); // pop "sp" off stack
    }
    lua_pop(L, 1); // pop "species" off stack
}

Noble_Abel_gas::
~Noble_Abel_gas() {}

int
Noble_Abel_gas::
s_eval_pressure(Gas_data &Q)
{
    int status;
    // printf("Noble_Abel_gas::s_eval_pressure()\n");
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    double nu_0 = s_covolume(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.p = nag_pressure(Q.rho, Q.T[0], R, nu_0);
    return SUCCESS;
}

int
Noble_Abel_gas::
s_eval_temperature(Gas_data &Q)
{
    int status;
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    double nu_0 = s_covolume(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.T[0] = nag_temperature(Q.rho, Q.p, R, nu_0);
    return SUCCESS;
}

int
Noble_Abel_gas::
s_eval_density(Gas_data &Q)
{
    int status;
    // printf("Noble_Abel_gas::s_eval_density()\n");
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    double nu_0 = s_covolume(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.rho = nag_density(Q.T[0], Q.p, R, nu_0);
    return SUCCESS;
}

double
Noble_Abel_gas::
s_gas_constant(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return calculate_gas_constant(Q.massf, M_);
}

double 
Noble_Abel_gas::
s_prho_ratio(const Gas_data &Q, int isp)
{ 
    double M = calculate_molecular_weight(Q.massf, M_);
    return Q.p/Q.rho*M/M_[isp];
}

double
Noble_Abel_gas::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    double nu_0 = s_covolume(Q, status);
    double dTdp_const_nu = (1.0/Q.rho - nu_0)/R;
    return dTdp_const_nu;
}

double
Noble_Abel_gas::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    double dTdnu_const_p = Q.p/R;
    return -1.0/(Q.rho*Q.rho)*dTdnu_const_p;
}

double
Noble_Abel_gas::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    double nu_0 = s_covolume(Q, status);
    double dpdnu_const_T = -R*Q.T[0]*pow((1/Q.rho - nu_0), -2);
    return -1.0/(Q.rho*Q.rho)*dpdnu_const_T;
}

double
Noble_Abel_gas::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    // CHECKME - this is just Dan's guess at what dpdrho_const_T for species isp is!
    double R_i = PC_R_u / M_[isp];
    double nu_0_i = nu_0_[isp];
    double rho_i = Q.rho * Q.massf[isp];
    double dpdnu_i_const_T = -R_i*Q.T[0]*pow((1/rho_i - nu_0_i), -2);
    return -1.0/(rho_i*rho_i)*dpdnu_i_const_T;
}

double
Noble_Abel_gas::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    double R = s_gas_constant(Q, status);
    double nu_0 = s_covolume(Q, status);
    double dTdp_const_nu = (1.0/Q.rho - nu_0)/R;
    return 1.0/dTdp_const_nu;
}

double
Noble_Abel_gas::
s_covolume(const Gas_data &Q, int &status)
{
    vector<double> molef;
    molef.resize(M_.size());
    convert_massf2molef(Q.massf, M_, molef);
    double nu_0 = 0.0;
    for ( size_t isp = 0; isp < Q.massf.size(); ++isp ) {
	nu_0 += molef[isp]*nu_0_[isp];
    }
    status = SUCCESS;
    return nu_0;
}

double nag_pressure(double rho, double T, double R, double nu_0)
{
    double nu_1 = 1.0/rho - nu_0;
    if (nu_1 <= 0) {
	cout << "A Noble-Abel gas cannot have a specific-volume smaller than the co-volume.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    // printf("nag_pressure: rho=%g nu_0=%g, R=%g, T=%g\n", rho, nu_0, R, T);
    return R*T/nu_1;
}

double nag_temperature(double rho, double p, double R, double nu_0)
{
    double nu_1 = 1/rho - nu_0;
    if (nu_1 <= 0) {
	cout << "A Noble-Abel gas cannot have a specific-volume smaller than the co-volume.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    return p*nu_1/R;
}

double nag_density(double T, double p, double R, double nu_0)
{
    // printf("nag_desnity: p=%g nu_0=%g, R=%g, T=%g\n", p, nu_0, R, T);
    double nu = R*T/p + nu_0;
    return 1/nu;
}
