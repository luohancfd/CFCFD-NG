// Author: Brendan T. O'Flaherty
//  based on noble-abel-gas-EOS code
// Version: 24-May-2008
//            Initial coding.
//

#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "van-der-waals-gas-EOS.hh"

using namespace std;

van_der_Waals_gas::
van_der_Waals_gas(lua_State *L)
    : Equation_of_state()
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "van_der_Waals_gas::van_der_Waals_gas():\n";
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
	    ost << "van_der_Waals_gas::van_der_Waals_gas()\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	double M = get_positive_value(L, -1, "M");
	M_.push_back(M);
	R_.push_back(PC_R_u/M);

	double T_c = get_value(L, -1, "T_c");
	double p_c = get_positive_value(L, -1, "p_c");
	
	double nu_0 = 0.125*PC_R_u*T_c/p_c;
	nu_0_.push_back(nu_0);

	double a = (27.0/64.0)*((PC_R_u*T_c)*(PC_R_u*T_c))/p_c;
	a_.push_back(a);
	
	lua_pop(L, 1); // pop "sp" off stack
    }
    lua_pop(L, 1); // pop "species" off stack
}

van_der_Waals_gas::
~van_der_Waals_gas() {}

int
van_der_Waals_gas::
s_eval_pressure(Gas_data &Q)
{
    int status;
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    double nu_0 = s_covolume(Q, status);
    if ( status != SUCCESS )
	return status;
    double a = s_a(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.p = vdwg_pressure(Q.rho, Q.T[0], R, nu_0, a);
    return SUCCESS;
}

int
van_der_Waals_gas::
s_eval_temperature(Gas_data &Q)
{
    int status;
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    double nu_0 = s_covolume(Q, status);
    if ( status != SUCCESS )
	return status;
    double a = s_a(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.T[0] = vdwg_temperature(Q.rho, Q.p, R, nu_0, a);
    return SUCCESS;
}

int
van_der_Waals_gas::
s_eval_density(Gas_data &Q)
{
    int status;
    double R = s_gas_constant(Q, status);
    if ( status != SUCCESS )
	return status;
    double nu_0 = s_covolume(Q, status);
    if ( status != SUCCESS )
	return status;
    double a = s_a(Q, status);
    if ( status != SUCCESS )
	return status;
    // else proceed
    Q.rho = vdwg_density(Q.T[0], Q.p, R, nu_0, a);
    return SUCCESS;
}

double
van_der_Waals_gas::
s_gas_constant(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return calculate_gas_constant(Q.massf, M_);
}

double 
van_der_Waals_gas::
s_prho_ratio(const Gas_data &Q, int isp)
{
    double M = calculate_molecular_weight(Q.massf, M_);
    return Q.p/Q.rho*M/M_[isp];
}

double
van_der_Waals_gas::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    double nu = 1/Q.rho;
    double nu_0 = s_covolume(Q, status);
    return (nu - nu_0)/R;
}

double
van_der_Waals_gas::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    double a = s_a(Q, status);
    double nu = 1/Q.rho;
    double nu_0 = s_covolume(Q, status);
    double dTdnu_const_p = (1/R)*(Q.p - a/(nu*nu*nu)*(nu + 2*nu_0));
    return -nu*nu*dTdnu_const_p;
}

double
van_der_Waals_gas::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q, status);
    double a = s_a(Q, status);
    double nu = 1/Q.rho;
    double nu_0 = s_covolume(Q, status);
    double dpdnu_const_T = -R*Q.T[0]/((nu - nu_0)*(nu - nu_0)) - 2*a/nu;
    return -nu*nu*dpdnu_const_T;
}

double
van_der_Waals_gas::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    // CHECKME - this is just Dan's guess at what dpdrho_const_T for species isp is!
    double R_i = PC_R_u / M_[isp];
    double a_i = a_[isp];
    double nu_i = 1/(Q.rho*Q.massf[isp]);
    double nu_0_i = nu_0_[isp];
    double dpdnu_i_const_T = -R_i*Q.T[0]/((nu_i - nu_0_i)*(nu_i - nu_0_i)) - 2*a_i/nu_i;
    return -nu_i*nu_i*dpdnu_i_const_T;
}


double
van_der_Waals_gas::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    double R = s_gas_constant(Q, status);
    double nu = 1/Q.rho;
    double nu_0 = s_covolume(Q, status);
    return R/(nu - nu_0);
}

double
van_der_Waals_gas::
s_covolume(const Gas_data &Q, int &status)
{
    vector<double> molef;
    molef.resize(M_.size());
    convert_massf2molef(Q.massf, M_, molef);

    double nu_0 = 0.0;
    for ( size_t isp = 0; isp < Q.massf.size(); ++isp ) {
	nu_0 += molef[isp]*nu_0_[isp];
    }
    
    double M = calculate_molecular_weight(Q.massf, M_);

    status = SUCCESS;
    return nu_0/M;
}

double
van_der_Waals_gas::
s_a(const Gas_data &Q, int &status)
{
    vector<double> molef;
    molef.resize(M_.size());
    convert_massf2molef(Q.massf, M_, molef);

    double a = 0.0;
    for ( size_t isp = 0; isp < Q.massf.size(); ++isp ) {
	a += molef[isp]*sqrt(a_[isp]);
    }
    double M = calculate_molecular_weight(Q.massf, M_);

    status = SUCCESS;
    return (a*a)/M;
}

double vdwg_pressure(double rho, double T, double R, double nu_0, double a)
{
    double nu_1 = 1/rho - nu_0;
    if (nu_1 <= 0) {
	cout << "Error in vdwg_pressure\n";
	cout << "nu = " << 1/rho << ", nu_0 = " << nu_0 << endl;;
	cout << "A van der Waals gas cannot have a specific-volume smaller than the co-volume.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    return R*T/nu_1 - a*rho*rho;
}

double vdwg_temperature(double rho, double p, double R, double nu_0, double a)
{
    double nu_1 = 1/rho - nu_0;
    if (nu_1 <= 0) {
	cout << "Error in vdwg_temperature\n";
	cout << "A van der Waals gas cannot have a specific-volume smaller than the co-volume.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    return (p + a*rho*rho)*nu_1/R;
}

double vdwg_density(double T, double p, double R, double nu_0, double a)
{
    double nu_old, nu, tol;
    nu = R*T/p + nu_0;
    tol = 1e-8;
    do {
	nu_old = nu;
	nu = R*T/(p + a/(nu_old*nu_old)) + nu_0;
    } while (fabs(nu_old - nu) > tol);
    
    return 1.0/nu;
}
