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

    lua_getfield(L, -1, "Cv0_coeffs"); // Bring Cv0 coeffs to TOS
    G_ = get_vector(L, -1, "G");
    lua_pop(L, 1); // pop Cv0 coeffs off stack

    lua_getfield(L, -1, "reference_state"); // Bring reference states to TOS
    lua_getfield(L, -1, "T0");
    T0 = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop T0 off stack
    lua_getfield(L, -1, "u0");
    u0 = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop u0 off stack
    lua_getfield(L, -1, "s0");
    s0 = lua_tonumber(L, -1);
    lua_pop(L, 1); // pop s0 off stack
    lua_pop(L, 1); // pop reference states off stack

    lua_pop(L, 1); // pop specific species table off stack
}

Dense_real_thermal_behaviour::
~Dense_real_thermal_behaviour() {}

int
Dense_real_thermal_behaviour::
s_decode_conserved_energy(Gas_data &Q, const vector<double> &rhoe)
{
    return tbm_decode_conserved_energy(Q.e, rhoe, Q.rho);
}

int
Dense_real_thermal_behaviour::
s_encode_conserved_energy(const Gas_data &Q, vector<double> &rhoe)
{
    return tbm_encode_conserved_energy(rhoe, Q.e, Q.rho);
}

double
Dense_real_thermal_behaviour::
s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    // This could be improved by using Univariate_functor
    // and Richardson extrapolation from the nm library.
    double dT = 0.01;
    double h1 = s_eval_enthalpy_isp(Q, EOS_, 0);
    Gas_data Q_peturb(Q);
    Q_peturb.T[0] += dT;
    // Recalculate the thermodynamic state at the perturbed temp and original pressure.
    EOS_->eval_density(Q_peturb);
    double h2 = s_eval_enthalpy_isp(Q_peturb, EOS_, 0);
    return (h2 - h1)/dT;
}

double
Dense_real_thermal_behaviour::
s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status)
{
    // This could be improved by using Univariate_functor
    // and Richardson extrapolation from the nm library.
    double dT = 0.01;
    double e1 = s_eval_energy_isp(Q, EOS_, 0);
    Gas_data Q_peturb(Q);
    Q_peturb.T[0] += dT;
    // Recalculate the thermodynamic state at the perturbed temp and original density.
    EOS_->eval_pressure(Q_peturb);
    double e2 = s_eval_energy_isp(Q_peturb, EOS_, 0);
    return (e2 - e1)/dT;
}

int
Dense_real_thermal_behaviour::
s_eval_energy(Gas_data &Q, Equation_of_state *EOS_)
{
    int status = SUCCESS;
    double e = s_eval_energy_isp(Q, EOS_, 0);
    Q.e[0] = e;
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
    double e = drtb_integral_Cv0(T0, Q.T[0], G_); // First integral, ideal gas at zero density
    e += EOS_->eval_integral_const_T_energy(Q); // Second integral, constant temperature
    e += u0;
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
    double s = drtb_integral_Cv0_over_T(T0, Q.T[0], G_); // First integral, ideal gas at zero density
    s += - R_*log(Q.rho);
    s += EOS_->eval_integral_const_T_entropy(Q); // Second integral, constant temperature
    s += s0;
    return s;
}

double drtb_integral_Cv0(double T0, double T, const std::vector<double> &G)
{
    // Integral of the classic fourth order polynomial ideal gas specific
    // heat capacity equation from the reference temperature T0 to T,
    // with respect to temperature.
    // Generated by sympy from the utility script Bender-gas-equations.py.
    // For use in the evaluation of internal energy.
    return G[1]*log(T) - G[1]*log(T0) + G[2]*T - G[2]*T0 + G[3]*pow(T, 2)/2
        - G[3]*pow(T0, 2)/2 + G[4]*pow(T, 3)/3 - G[4]*pow(T0, 3)/3
        + G[5]*pow(T, 4)/4 - G[5]*pow(T0, 4)/4 + G[6]*pow(T, 5)/5
        - G[6]*pow(T0, 5)/5;
}

double drtb_integral_Cv0_over_T(double T0, double T, const std::vector<double> &G)
{
    // Integral of the classic fourth order polynomial ideal gas specific
    // heat capacity equation divided by T from the reference temperature T0 to T,
    // with respect to temperature.
    // Generated by sympy from the utility script Bender-gas-equations.py.
    // For use in the evaluation of entropy.
    return G[1]/T0 - G[1]/T + G[2]*log(T) - G[2]*log(T0) + G[3]*T - G[3]*T0
        + G[4]*pow(T, 2)/2 - G[4]*pow(T0, 2)/2 + G[5]*pow(T, 3)/3
        - G[5]*pow(T0, 3)/3 + G[6]*pow(T, 4)/4 - G[6]*pow(T0, 4)/4;
}
