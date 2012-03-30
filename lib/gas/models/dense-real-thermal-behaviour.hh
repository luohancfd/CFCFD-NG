/** \file dense-real-thermal-behaviour.hh
 *  \ingroup gas
 *
 *  \author Peter Blyton
 *  \version 26-Mar-2012
 */

#ifndef DENSE_REAL_THERMAL_BEHAVIOUR_HH
#define DENSE_REAL_THERMAL_BEHAVIOUR_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "thermal-behaviour-model.hh"

/**
 *  \brief Class for implementing a real gas thermal behaviour model.
 *  Evaluates energy and entropy by integration of the equation of state
 *  and the ideal gas specific heat capacity, as described by:
 *
 *  Reynolds, WC (1979). Thermodynamic Properties in SI:
 *  Graphs, tables, and computational equations for forty substances.
 *  Department of Mechanical Engineering, Stanford University.
 *
 *  The classic fourth order polynomial ideal gas specific heat equation is
 *  used, this is equation C-6 published by Reynolds.
 *
 *  As there are no mixing rules implemented for this thermal behaviour model,
 *  it should only be used with single species gases.
 */

class Dense_real_thermal_behaviour : public Thermal_behaviour_model {
public:
    Dense_real_thermal_behaviour(lua_State *L);
    ~Dense_real_thermal_behaviour();

private:
    std::vector<double> G_; // Coefficients for the Cv0 equation.
    double T0, u0, s0;      // Thermodynamic reference state parameters.
    double R_;
    std::vector<double> M_;

    int s_decode_conserved_energy(Gas_data &Q, const std::vector<double> &rhoe);
    int s_encode_conserved_energy(const Gas_data &Q, std::vector<double> &rhoe);
    double s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status);
    double s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status);
    int s_eval_energy(Gas_data &Q, Equation_of_state *EOS_);
    int s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_);
    double s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
};

double drtb_integral_Cv0(double T0, double T, const std::vector<double> &G);
double drtb_integral_Cv0_over_T(double T0, double T, const std::vector<double> &G);

#endif
