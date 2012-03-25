/** \file Bender-real-gas-EOS.cxx
 *  \ingroup gas
 *
 *  \author Peter Blyton
 *  \version 13-Mar-2012
 */

#include <iostream>
#include <sstream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "Bender-real-gas-EOS.hh"

using namespace std;

const double BRG_PROPS_CONVERGE_TOL = 0.001;
const double BRG_MAX_ITERATIONS = 10;

Bender_real_gas::
Bender_real_gas(lua_State *L)
    : Equation_of_state()
{
    lua_getglobal(L, "species"); // Put the list of species table to TOS
    
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "Bender_real_gas::Bender_real_gas():\n";
        ost << "Error in the declaration of species: a table is expected.\n";
        input_error(ost);
    }

    if ( lua_objlen(L, -1) != 1 ) {
        ostringstream ost;
        ost << "Bender_real_gas::Bender_real_gas():\n";
        ost << "Error in the declaration of species: single species expected.\n";
        input_error(ost);
    }

    lua_rawgeti(L, -1, 1); // Put the specific species name to TOS
    const char* sp = luaL_checkstring(L, -1);
    lua_pop(L, 1); // pop specific species name off stack
    lua_pop(L, 1); // pop list of species table off stack
    

    lua_getglobal(L, sp); // Now bring the specific species table to TOS
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "Bender_real_gas::Bender_real_gas()\n";
        ost << "Error locating information table for species: " << sp << endl;
        input_error(ost);
    }

    double M;
    M = get_positive_value(L, -1, "M");
    M_.push_back(M);
    R_ = PC_R_u/M;

    lua_getfield(L, -1, "Bender_EOS_coeffs"); // Bring EOS coeffs to TOS
    A_ = get_vector(L, -1, "A");
    lua_pop(L, 1); // pop EOS coeffs off stack
    lua_pop(L, 1); // pop specific species table off stack
}

Bender_real_gas::
~Bender_real_gas() {}

int
Bender_real_gas::
s_eval_pressure(Gas_data &Q)
{
    Q.p = brg_pressure(Q.rho, Q.T[0], R_, A_);
    return SUCCESS;
}

int
Bender_real_gas::
s_eval_temperature(Gas_data &Q)
{
    int status;
    Q.T[0] = brg_temperature(Q.rho, Q.p, R_, A_, status);
    return status;
}

int
Bender_real_gas::
s_eval_density(Gas_data &Q)
{
    int status;
    Q.rho = brg_density(Q.T[0], Q.p, R_, A_, status);
    return status;
}

double
Bender_real_gas::
s_gas_constant(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    // Return a gas constant for use by the CFD code, which assumes
    // ideal gas. For this equation of state, it is not simply R_.
    return Q.p / (Q.rho * Q.T[0]);
}

double 
Bender_real_gas::
s_prho_ratio(const Gas_data &Q, int isp)
{
    return brg_pressure(Q.rho, Q.T[0], R_, A_)/Q.rho;
}

double
Bender_real_gas::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    return 1.0/brg_dpdT(Q.rho, Q.T[0], R_, A_);
}

double
Bender_real_gas::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    // Should be able to derive this derivative analytically, or do it numerically.
    // For now, just returning derivative of the ideal gas term.
    return (-1.0*Q.p)/(R_*Q.rho*Q.rho);
}

double
Bender_real_gas::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return brg_dpdrho(Q.T[0], Q.p, R_, A_);
}

double
Bender_real_gas::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    status = SUCCESS;
    return brg_dpdrho(Q.T[0], Q.p, R_, A_);
}

double
Bender_real_gas::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    status = SUCCESS;
    return brg_dpdT(Q.rho, Q.T[0], R_, A_);
}

double brg_pressure(double rho, double T, double R, const std::vector<double> &A)
{
    // Bender real gas equation for pressure, given density and temperature.
    // As transcribed from Reynold's text and manipulated a little by PJ.
    double Tinv = 1.0/T;
    double p = rho * R * T;
    double rhos = rho * rho;
    p += rhos * (A[1]*T + A[2] + Tinv*(A[3] + Tinv*(A[4] + Tinv*A[5])));
    rhos *= rho;
    p += rhos * (T*A[6] + A[7] + Tinv*A[8]);
    rhos *= rho;
    p += rhos * (T*A[9] + A[10]);
    rhos *= rho;
    p += rhos * (T*A[11] + A[12]);
    rhos *= rho;
    p += rhos * A[13];
    // last two terms
    rhos = rho * rho * rho;
    double term1 = rhos * Tinv*Tinv*(A[14] + Tinv*(A[15] + Tinv*A[16]));
    rhos *= rho * rho;
    double term2 = rhos * Tinv*Tinv*(A[17] + Tinv*(A[18] + Tinv*A[19]));
    p += (term1 + term2) * exp(-A[20]*rho*rho);
    return p;
}

double brg_dpdT(double rho, double T, double R, const std::vector<double> &A)
{
    // Derivative (with respect to T) of brg_pressure.
    // This is needed for the integrals for internal energy and entropy.
    // Also used to iteratively solve the EOS for temperature.
    double Tinv = 1.0/T;
    double dpdT = rho * R;
    double rhos = rho * rho;
    dpdT += rhos * (A[1] - Tinv*Tinv*(A[3] + Tinv*(2.0*A[4] + Tinv*3.0*A[5])));
    rhos *= rho;
    dpdT += rhos * (A[6] - Tinv*Tinv*A[8]);
    rhos *= rho;
    dpdT += rhos * A[9];
    rhos *= rho;
    dpdT += rhos * A[11];
    // last two terms
    rhos = rho * rho * rho;
    double term1 = (-1.0)*rhos*Tinv*Tinv*Tinv*(2.0*A[14]+Tinv*(3.0*A[15]+Tinv*4.0*A[16]));
    rhos *= rho * rho;
    double term2 = (-1.0)*rhos*Tinv*Tinv*Tinv*(2.0*A[17]+Tinv*(3.0*A[18]+Tinv*4.0*A[19]));
    dpdT += (term1 + term2) * exp(-A[20]*rho*rho);
    return dpdT;
}

double brg_dpdrho(double rho, double T, double R, const std::vector<double> &A)
{
    // Derivative (with respect to T) of brg_pressure.
    // This is needed for the integrals for internal energy and entropy.
    // Also used to iteratively solve the EOS for density.
    double Tinv = 1.0/T;
    double dpdrho = R * T;
    dpdrho += 2.0 * rho * (A[1]*T + A[2] + Tinv*(A[3] + Tinv*(A[4] + Tinv*A[5])));
    double rhos = rho * rho;
    dpdrho += 3.0 * rhos * (T*A[6] + A[7] + Tinv*A[8]);
    rhos *= rho;
    dpdrho += 4.0 * rhos * (T*A[9] + A[10]);
    rhos *= rho;
    dpdrho += 5.0 * rhos * (T*A[11] + A[12]);
    rhos *= rho;
    dpdrho += 6.0 * rhos * A[13];
    // last two terms, un-differentiated
    rhos = rho * rho * rho;
    double term1 = rhos * Tinv*Tinv*(A[14] + Tinv*(A[15] + Tinv*A[16]));
    rhos *= rho * rho;
    double term2 = rhos * Tinv*Tinv*(A[17] + Tinv*(A[18] + Tinv*A[19]));
    // last two terms, differentiated
    rhos = rho * rho;
    double dterm1_drho = 3.0 * rhos * Tinv*Tinv*(A[14] + Tinv*(A[15] + Tinv*A[16]));
    rhos *= rho * rho;
    double dterm2_drho = 5.0 * rhos * Tinv*Tinv*(A[17] + Tinv*(A[18] + Tinv*A[19]));
    // Differentiate exponential and sub in terms:
    dpdrho += exp(-A[20]*rho*rho) * (-2.0*A[20]*rho*(term1 + term2) + (dterm1_drho + dterm2_drho));
    return dpdrho;
}

double brg_temperature(double rho, double p, double R, const std::vector<double> &A, int &status)
{
    // Evaluate temperature from density and pressure. Must solve EOS iteratively.
    // Temperature from ideal gas EOS is a good first guess.
    double T_guess = p / (rho * R);
    int i = 0;
    double T_temp;
    status = SUCCESS;
    do {
        if (i >= BRG_MAX_ITERATIONS) {
            status = ITERATION_ERROR;
            cout << "brg_temperature():\n";
            cout << "    Iterations did not converge.\n";
            return T_guess;
        }
        T_temp = T_guess;
        T_guess -= (brg_pressure(rho, T_guess, R, A) - p)/brg_dpdT(rho, T_guess, R, A);
        i++;
    } while (fabs(T_temp - T_guess) > BRG_PROPS_CONVERGE_TOL);
    return T_guess;
}

double brg_density(double T, double p, double R, const std::vector<double> &A, int &status)
{
    // Evaluate density from temperature and pressure. Must solve EOS iteratively.
    // Density from ideal gas EOS is a good first guess.
    double rho_guess = p / (R * T);
    int i = 0;
    double rho_temp;
    status = SUCCESS;
    do {
        if (i >= BRG_MAX_ITERATIONS) {
            status = ITERATION_ERROR;
            cout << "brg_density():\n";
            cout << "    Iterations did not converge.\n";
            return rho_guess;
        }
        rho_temp = rho_guess;
        rho_guess -= (brg_pressure(rho_guess, T, R, A) - p)/brg_dpdrho(rho_guess, T, R, A);
        i++;
    } while (fabs(rho_temp - rho_guess) > BRG_PROPS_CONVERGE_TOL);
    return rho_guess;
}

