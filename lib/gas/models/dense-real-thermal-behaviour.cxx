/** \file dense-real-thermal-behaviour.cxx
 *  \ingroup gas
 *
 *  \author Peter Blyton
 *  \version 26-Mar-2012
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "dense-real-thermal-behaviour.hh"

using namespace std;

const double DRTB_PROPS_CONVERGE_TOL = 0.001;
const int DRTB_MAX_ITERATIONS = 10;

Dense_real_thermal_behaviour::
Dense_real_thermal_behaviour(lua_State *L)
    : Thermal_behaviour_model()
{
    lua_getglobal(L, "species"); // Put the list of species table to TOS

    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "Dense_real_thermal_behaviour::Dense_real_thermal_behaviour():\n";
        ost << "Error in the declaration of species: a table is expected.\n";
        input_error(ost);
    }

    if ( lua_objlen(L, -1) != 1 ) {
        ostringstream ost;
        ost << "Dense_real_thermal_behaviour::Dense_real_thermal_behaviour():\n";
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
        ost << "Dense_real_thermal_behaviour::Dense_real_thermal_behaviour()\n";
        ost << "Error locating information table for species: " << sp << endl;
        input_error(ost);
    }

    double M = get_positive_value(L, -1, "M");
    M_.push_back(M);
    R_ = PC_R_u/M;

    lua_getfield(L, -1, "Cp0_coeffs"); // Bring Cp0 coeffs to TOS
    G_ = get_vector(L, -1, "G");
    lua_getfield(L, -1, "T_low");
    T_low_ = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop T_low off stack
    lua_getfield(L, -1, "T_high");
    T_high_ = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop T_high off stack
    lua_pop(L, 1); // pop Cp0 coeffs off stack

    lua_getfield(L, -1, "reference_state"); // Bring reference states to TOS
    lua_getfield(L, -1, "T0");
    T0_ = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop T0 off stack
    lua_getfield(L, -1, "u0");
    u0_ = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop u0 off stack
    lua_getfield(L, -1, "s0");
    s0_ = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop s0 off stack
    lua_pop(L, 1); // pop reference states off stack

    lua_pop(L, 1); // pop specific species table off stack
}

Dense_real_thermal_behaviour::
~Dense_real_thermal_behaviour() {}

int
Dense_real_thermal_behaviour::
check_temperature_valid(const double &T)
{
    // Check that temperature is within range of validity of Cv0 data.
    if (T < T_low_ || T > T_high_) {
        cout << "Dense_real_thermal_behaviour::check_temperature_valid():\n"
             << "   Temperature = " << T << " out of validity range!\n";
        return FAILURE;
    } else {
        return SUCCESS;
    }
}

double
Dense_real_thermal_behaviour::
s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    double Cp = s_dedT_const_v(Q, EOS_, status);
    // Use chain rule to find (dv/dT)p
    double dvdT_const_p = -pow(Q.rho, -2)/EOS_->dTdrho_const_p(Q, status);
    // Apply the Mayer relation, equation 12-45 from Cengel and Boles Thermodynamics, 6e.
    Cp += Q.T[0]*dvdT_const_p/EOS_->dTdp_const_rho(Q, status);
    return Cp;
}

double
Dense_real_thermal_behaviour::
s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    // Evaluation derivative analytically.
    check_temperature_valid(Q.T[0]);
    double Cv = drtb_Cp0(Q.T[0], G_) - R_;
    Cv += EOS_->eval_dintegral_const_T_energy_dT(Q);
    return Cv;
}

int
Dense_real_thermal_behaviour::
s_eval_energy(Gas_data &Q, Equation_of_state *EOS_)
{
    int status = SUCCESS;
    Q.e[0] = s_eval_energy_isp(Q, EOS_, 0);
    return status;
}

int
Dense_real_thermal_behaviour::
s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_)
{
    // From the given energy and density, calculate temperature.
    int i = 0;
    // Cannot simply use T_guess = e/Cv0 as the integral of EOS has large effect.
    double T_guess = 300.0;
    int status;
    do { // Newtons method to solve for temperature.
        if (i >= DRTB_MAX_ITERATIONS) {
            cout << "Dense_real_thermal_behaviour::s_eval_temperature():\n";
            cout << "    Iterations did not converge.\n";
            exit(ITERATION_ERROR);
        }
        Q.T[0] = T_guess;
        T_guess -= (s_eval_energy_isp(Q, EOS_, 0) - Q.e[0])/s_dedT_const_v(Q, EOS_, status);
        i++;
    } while (fabs(Q.T[0] - T_guess) > DRTB_PROPS_CONVERGE_TOL);
    return SUCCESS;
}

double
Dense_real_thermal_behaviour::
s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    // Evaluate the internal energy of a real gas using integration.
    // See Equation 4, Section 2 of Reynolds, WC (1979). Thermodynamic Properties in SI.
    check_temperature_valid(Q.T[0]);
    double e = drtb_integral_Cp0(T0_, Q.T[0], G_) - R_*(Q.T[0] - T0_); // First integral, ideal gas at zero density
    e += EOS_->eval_integral_const_T_energy(Q); // Second integral, constant temperature
    e += u0_;
    return e;
}

double
Dense_real_thermal_behaviour::
s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    return s_eval_energy_isp(Q, EOS_, isp) + Q.p/Q.rho;
}

double
Dense_real_thermal_behaviour::
s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp)
{
    // Evaluate the entropy of a real gas using integration.
    // See Equation 8, Section 2 of Reynolds, WC (1979). Thermodynamic Properties in SI.
    check_temperature_valid(Q.T[0]);
    double s = drtb_integral_Cp0_over_T(T0_, Q.T[0], G_) - R_*log(Q.T[0]/T0_); // First integral, ideal gas at zero density
    s += - R_*log(Q.rho);
    s += EOS_->eval_integral_const_T_entropy(Q); // Second integral, constant temperature
    s += s0_;
    return s;
}

double drtb_Cp0(double T, const std::vector<double> &G)
{
    // Isobaric ideal gas specific heat capacity equation.
    // For use in the evaluation of Cv.
    return G[1]/T + G[2] + G[3]*T + G[4]*pow(T, 2) + G[5]*pow(T, 3) + G[6]*pow(T, 4);
}

double drtb_integral_Cp0(double T0, double T, const std::vector<double> &G)
{
    // Integral of heat capacity equation from the reference temperature T0 to T,
    // with respect to temperature.
    // For use in the evaluation of internal energy.
    return G[1]*log(T/T0) + G[2]*(T - T0)
        + G[3]*(pow(T, 2) - pow(T0, 2))/2
        + G[4]*(pow(T, 3) - pow(T0, 3))/3
        + G[5]*(pow(T, 4) - pow(T0, 4))/4
        + G[6]*(pow(T, 5) - pow(T0, 5))/5;
}

double drtb_integral_Cp0_over_T(double T0, double T, const std::vector<double> &G)
{
    // Integral of heat capacity equation divided by T from the reference temperature T0 to T,
    // with respect to temperature.
    // For use in the evaluation of entropy.
    return G[1]*(1/T0 - 1/T)
        + G[2]*log(T/T0) + G[3]*(T - T0)
        + G[4]*(pow(T, 2) - pow(T0, 2))/2
        + G[5]*(pow(T, 3) - pow(T0, 3))/3
        + G[6]*(pow(T, 4) - pow(T0, 4))/4;
}
