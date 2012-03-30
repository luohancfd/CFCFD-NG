/** \file Bender-real-gas-EOS.hh
 *  \ingroup gas
 *
 *  \author Peter Blyton
 *  \version 13-Mar-2012
 */

#ifndef BENDER_REAL_GAS_EOS_HH
#define BENDER_REAL_GAS_EOS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "equation-of-state.hh"

/**
 *  \brief Class for implementing the Bender equation of state described by:
 *
 *  Reynolds, WC (1979). Thermodynamic Properties in SI:
 *  Graphs, tables, and computational equations for forty substances.
 *  Department of Mechanical Engineering, Stanford University.
 *
 *  The Bender equation of state is equation P-3 published by Reynolds.
 *
 *  As there are no mixing rules implemented for this equation of state,
 *  it should only be used with single species gases.
 */

class Bender_real_gas : public Equation_of_state {
public:
    Bender_real_gas(lua_State *L);
    ~Bender_real_gas();

private:
    std::vector<double> A_; // Coefficients for Bender EOS.
    double R_;
    int s_eval_pressure(Gas_data &Q);
    int s_eval_temperature(Gas_data &Q);
    int s_eval_density(Gas_data &Q);
    double s_gas_constant(const Gas_data &Q, int &status);
    double s_prho_ratio(const Gas_data &Q, int isp);
    double s_dTdp_const_rho(const Gas_data &Q, int &status);
    double s_dTdrho_const_p(const Gas_data &Q, int &status);
    double s_dpdrho_const_T(const Gas_data &Q, int &status);
    double s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status);
    double s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status);
    double s_integral_const_T_energy(const Gas_data &Q);
    double s_integral_const_T_entropy(const Gas_data &Q);
};

double brg_pressure(double rho, double T, double R, const std::vector<double> &A);
double brg_dpdT(double rho, double T, double R, const std::vector<double> &A);
double brg_dpdrho(double rho, double T, double R, const std::vector<double> &A);
double brg_temperature(double rho, double p, double R, const std::vector<double> &A, int &status);
double brg_density(double T, double p, double R, const std::vector<double> &A, int &status);
double brg_integral_const_T_energy(double rho, double T, double R, const std::vector<double> &A);
double brg_integral_const_T_entropy(double rho, double T, double R, const std::vector<double> &A);

#endif
