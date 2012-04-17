/** \file Bender-EOS.cxx
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
#include "Bender-EOS.hh"

using namespace std;

const double BENDER_PROPS_CONVERGE_TOL = 0.001;
const int BENDER_MAX_ITERATIONS = 10;

Bender_EOS::
Bender_EOS(lua_State *L)
    : Equation_of_state()
{
    lua_getglobal(L, "species"); // Put the list of species table to TOS

    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "Bender_EOS::Bender_EOS():\n";
        ost << "Error in the declaration of species: a table is expected.\n";
        input_error(ost);
    }

    if ( lua_objlen(L, -1) != 1 ) {
        ostringstream ost;
        ost << "Bender_EOS::Bender_EOS():\n";
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
        ost << "Bender_EOS::Bender_EOS()\n";
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

Bender_EOS::
~Bender_EOS() {}

int
Bender_EOS::
s_eval_pressure(Gas_data &Q)
{
    Q.p = Bender_pressure(Q.rho, Q.T[0], R_, A_);
    return SUCCESS;
}

int
Bender_EOS::
s_eval_temperature(Gas_data &Q)
{
    int status;
    Q.T[0] = Bender_temperature(Q.rho, Q.p, R_, A_, status);
    return status;
}

int
Bender_EOS::
s_eval_density(Gas_data &Q)
{
    int status;
    Q.rho = Bender_density(Q.T[0], Q.p, R_, A_, status);
    return status;
}

double
Bender_EOS::
s_gas_constant(const Gas_data &Q, int &status)
{   // Return the "equivalent" R for the ideal gas equation, needed by CFD code.
    status = SUCCESS;
    return Q.p / (Q.rho * Q.T[0]);
}

double 
Bender_EOS::
s_prho_ratio(const Gas_data &Q, int isp)
{
    return Bender_pressure(Q.rho, Q.T[0], R_, A_)/Q.rho;
}

double
Bender_EOS::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    return 1.0/Bender_dpdT(Q.rho, Q.T[0], R_, A_);
}

double
Bender_EOS::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{   // Use the partial derivative cyclic relation.
    return -s_dTdp_const_rho(Q, status)*s_dpdrho_const_T(Q, status);
}

double
Bender_EOS::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    status = SUCCESS;
    return Bender_dpdrho(Q.rho, Q.T[0], R_, A_);
}

double
Bender_EOS::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    status = SUCCESS;
    return Bender_dpdrho(Q.rho, Q.T[0], R_, A_);
}

double
Bender_EOS::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    status = SUCCESS;
    return Bender_dpdT(Q.rho, Q.T[0], R_, A_);
}

double
Bender_EOS::
s_integral_const_T_energy(const Gas_data &Q)
{
    return Bender_integral_const_T_energy(Q.rho, Q.T[0], A_);
}

double
Bender_EOS::
s_dintegral_const_T_energy_dT(const Gas_data &Q)
{
    return Bender_dintegral_const_T_energy_dT(Q.rho, Q.T[0], A_);
}

double
Bender_EOS::
s_integral_const_T_entropy(const Gas_data &Q)
{
    return Bender_integral_const_T_entropy(Q.rho, Q.T[0], A_);
}

double Bender_pressure(double rho, double T, double R, const std::vector<double> &A)
{
    // Bender p-rho-T equation of state.
    double Tinv = 1/T;
    double calc_exp = exp(-A[20]*pow(rho, 2));
    return rho*R*T
        + pow(rho, 2)*(A[1]*T + A[2] + Tinv*(A[3] + Tinv*(A[4] + A[5]*Tinv)))
        + pow(rho, 3)*(A[6]*T + A[7] + A[8]*Tinv)
        + pow(rho, 4)*(A[9]*T + A[10])
        + pow(rho, 5)*(A[11]*T + A[12])
        + pow(rho, 6)*A[13]
        + pow(rho, 3)*pow(T, -2)*(A[14] + Tinv*(A[15] + A[16]*Tinv))*calc_exp
        + pow(rho, 5)*pow(T, -2)*(A[17] + Tinv*(A[18] + A[19]*Tinv))*calc_exp;
}

double Bender_dpdT(double rho, double T, double R, const std::vector<double> &A)
{
    // Derivative of pressure p-rho-T equation w.r.t temperature.
    // This is needed for the integrals for internal energy and entropy.
    // Also used to iteratively solve the EOS for temperature using Newtons method.
    double Tinv = 1/T;
    double calc_exp = exp(-A[20]*pow(rho, 2));
    return rho*R
        + pow(rho, 2)*(A[1] - pow(T, -2)*(A[3] + Tinv*(2*A[4] + 3*A[5]*Tinv)))
        + pow(rho, 3)*(A[6] - A[8]*pow(T, -2))
        + pow(rho, 4)*A[9] + pow(rho, 5)*A[11]
        - pow(rho, 3)*pow(T, -3)*(2*A[14] + Tinv*(3*A[15] + 4*A[16]*Tinv))*calc_exp
        - pow(rho, 5)*pow(T, -3)*(2*A[17] + Tinv*(3*A[18] + 4*A[19]*Tinv))*calc_exp;
}

double Bender_dpdrho(double rho, double T, double R, const std::vector<double> &A)
{
    // Derivative of pressure p-rho-T equation w.r.t density.
    // Used to iteratively solve the EOS for density using Newtons method.
    double Tinv = 1/T;
    double calc_exp = exp(-A[20]*pow(rho, 2));
    return R*T
        + 2*rho*(A[1]*T + A[2] + Tinv*(A[3] + Tinv*(A[4] + A[5]*Tinv)))
        + 3*pow(rho, 2)*(A[6]*T + A[7] + A[8]*Tinv)
        + 4*pow(rho, 3)*(A[9]*T + A[10])
        + 5*pow(rho, 4)*(A[11]*T + A[12])
        + 6*pow(rho, 5)*A[13]
        + (3*pow(rho, 2) - 2*rho*A[20]*pow(rho, 3))*pow(T, -2)*(A[14] + Tinv*(A[15] + A[16]*Tinv))*calc_exp
        + (5*pow(rho, 4) - 2*rho*A[20]*pow(rho, 5))*pow(T, -2)*(A[17] + Tinv*(A[18] + A[19]*Tinv))*calc_exp;
}

double Bender_temperature(double rho, double p, double R, const std::vector<double> &A, int &status)
{
    // Evaluate temperature from density and pressure. Must solve EOS iteratively.
    // Temperature from ideal gas EOS is a good first guess.
    double T_guess = p / (rho * R);
    int i = 0;
    double T_temp;
    status = SUCCESS;
    do { // Newtons method to solve for temperature.
        if (i >= BENDER_MAX_ITERATIONS) {
            status = ITERATION_ERROR;
            cout << "Bender_temperature():\n";
            cout << "    Iterations did not converge.\n";
            return T_guess;
        }
        T_temp = T_guess;
        T_guess -= (Bender_pressure(rho, T_guess, R, A) - p)/Bender_dpdT(rho, T_guess, R, A);
        i++;
    } while (fabs(T_temp - T_guess) > BENDER_PROPS_CONVERGE_TOL);
    return T_guess;
}

double Bender_density(double T, double p, double R, const std::vector<double> &A, int &status)
{
    // Evaluate density from temperature and pressure. Must solve EOS iteratively.
    // Density from ideal gas EOS is a good first guess.
    double rho_guess = p / (R * T);
    int i = 0;
    double rho_temp;
    status = SUCCESS;
    do { // Newtons method to solve for density.
        if (i >= BENDER_MAX_ITERATIONS) {
            status = ITERATION_ERROR;
            cout << "Bender_density():\n";
            cout << "    Iterations did not converge.\n";
            return rho_guess;
        }
        rho_temp = rho_guess;
        rho_guess -= (Bender_pressure(rho_guess, T, R, A) - p)/Bender_dpdrho(rho_guess, T, R, A);
        i++;
    } while (fabs(rho_temp - rho_guess) > BENDER_PROPS_CONVERGE_TOL);
    return rho_guess;
}

double Bender_integral_const_T_energy(double rho, double T, const std::vector<double> &A)
{
    // Integral of (rho^-2)*(P - T(dP/dT)) from zero to current density, w.r.t density.
    // For use in the evaluation of internal energy.
    double Tinv = 1/T;
    double calc_exp = exp(-A[20]*pow(rho, 2));
    return rho*(A[2] + Tinv*(2*A[3] + Tinv*(3*A[4] + 4*Tinv*A[5])))
        + pow(rho, 2)*(A[7] + 2*A[8]*Tinv)/2
        + pow(rho, 3)*A[10]/3
        + pow(rho, 4)*A[12]/4
        + pow(rho, 5)*A[13]/5
        + pow(T, -2)*(3*A[14] + Tinv*(4*A[15] + 5*A[16]*Tinv))*(1 - calc_exp)/(2*A[20])
        + pow(T, -2)*(3*A[17] + Tinv*(4*A[18] + 5*A[19]*Tinv))*(1 - (A[20]*pow(rho, 2) + 1)*calc_exp)/(2*pow(A[20], 2));
}

double Bender_dintegral_const_T_energy_dT(double rho, double T, const std::vector<double> &A)
{
    // Derivative of Bender_integral_const_T_energy w.r.t temperature.
    // For use in the evaluation of Cv.
    double Tinv = 1/T;
    double calc_exp = exp(-A[20]*pow(rho, 2));
    return -rho*pow(T, -2)*(2*A[3] + Tinv*(6*A[4] + 12*Tinv*A[5]))
        - pow(rho, 2)*A[8]*pow(T, -2)
        - pow(T, -3)*(6*A[14] + Tinv*(12*A[15] + 20*Tinv*A[16]))*(1 - calc_exp)/(2*A[20])
        - pow(T, -3)*(6*A[17] + Tinv*(12*A[18] + 20*Tinv*A[19]))*(1 - (A[20]*pow(rho, 2) + 1)*calc_exp)/(2*pow(A[20], 2));
}

double Bender_integral_const_T_entropy(double rho, double T, const std::vector<double> &A)
{
    // Integral of (rho^-2)*(rho*R - (dP/dT)) from zero to current density, w.r.t density.
    // For use in the evaluation of entropy.
    double Tinv = 1/T;
    double calc_exp = exp(-A[20]*pow(rho, 2));
    return rho*(-A[1] + pow(T, -2)*(A[3] + Tinv*(2*A[4] + 3*Tinv*A[5])))
        + pow(rho, 2)*(-A[6] + A[8]*pow(T, -2))/2
        - pow(rho, 3)*A[9]/3
        - pow(rho, 4)*A[11]/4
        + pow(T, -3)*(2*A[14] + Tinv*(3*A[15] + 4*A[16]*Tinv))*(1 - calc_exp)/(2*A[20])
        + pow(T, -3)*(2*A[17] + Tinv*(3*A[18] + 4*A[19]*Tinv))*(1 - (A[20]*pow(rho, 2) + 1)*calc_exp)/(2*pow(A[20], 2));
}
