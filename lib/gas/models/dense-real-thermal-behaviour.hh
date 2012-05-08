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
 *  Evaluates energy, entropy, Cv and Cp by integration of an equation of state
 *  and an isobaric ideal gas heat capacity equation, Cp0, of the form:
 *
 *  Cp0 = G2/T + G3 + G4*T + G5*T^2 + G6*T^3 + G7*T^4 [J/kg.K]
 *
 *  This form is used as it fits with the coefficients provided by
 *  McLinden (1989, 1992, 1995, 1997), Reynolds (1979), Platzer (1990), Polt (1992).
 */

class Dense_real_thermal_behaviour : public Thermal_behaviour_model {
public:
    Dense_real_thermal_behaviour(lua_State *L);
    ~Dense_real_thermal_behaviour();

private:
    std::vector<double> G_; // Coefficients for the Cp0 equation.
    double T_low_, T_high_;
    double T0_, u0_, s0_;   // Thermodynamic reference state parameters.
    double R_;
    std::vector<double> M_;

    int check_temperature_valid(const double &T);
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

double drtb_Cp0(double T, const std::vector<double> &G);
double drtb_integral_Cp0(double T0, double T, const std::vector<double> &G);
double drtb_integral_Cp0_over_T(double T0, double T, const std::vector<double> &G);

#endif
