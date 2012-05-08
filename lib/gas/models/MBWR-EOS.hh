/** \file MBWR-EOS.hh
 *  \ingroup gas
 *
 *  \author Peter Blyton
 *  \version 18-Apr-2012
 */

#ifndef MBWR_EOS_HH
#define MBWR_EOS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "equation-of-state.hh"

/**
 *  \brief Class for implementing the modified Benedict-Webb-Rubin equation of state.
 *
 *  The MBWR equation of state is the modification of the BWR equation proposed by:
 *
 *  Jacobsen, RT and RB Stewart (1973). "Thermodynamic properties of nitrogen including
 *  liquid and vapor phases from 63K to 2000K with pressures to 10,000 bar".
 *
 *  As there are no mixing rules implemented for this equation of state,
 *  it should only be used with single species gases.
 *
 *  The base gas model and gas data classes make us of SI units for all properties,
 *  however it is conventional for the MBWR equation and published coefficients to
 *  use the following units (convention of Jacobsen and Stewart, McLinden, DuPont):
 *
 *  pressure (bar)
 *  density (mol/L, or mol/(dm^3), 1L = 1 cubic decimetre)
 *  temperature (K)
 *  molecular mass (kg/mol)
 *  universal gas constant (L.bar/(mol.K) = J/(100 mol.K))
 *
 *  Unit conversion is performed in-function so that all interfaces remain in
 *  SI units. One condition is that the MBWR equation of state coefficients
 *  provided by the species file must be suitable for these conventional units.
 */

class MBWR_EOS : public Equation_of_state {
public:
    MBWR_EOS(lua_State *L);
    ~MBWR_EOS();

private:
    std::vector<double> b_; // Coefficients for MBWR EOS
    double rho_c_; // Critical density
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
    double s_dintegral_const_T_energy_dT(const Gas_data &Q);
    double s_integral_const_T_entropy(const Gas_data &Q);
};

double MBWR_pressure(double rho, double T, double rho_c, double M, const std::vector<double> &b);
double MBWR_dpdT(double rho, double T, double rho_c, double M, const std::vector<double> &b);
double MBWR_dpdrho(double rho, double T, double rho_c, double M, const std::vector<double> &b);
double MBWR_temperature(double rho, double p, double rho_c, double M, const std::vector<double> &b, int &status);
double MBWR_density(double T, double p, double rho_c, double M, const std::vector<double> &b, int &status);
double MBWR_integral_const_T_energy(double rho, double T, double rho_c, double M, const std::vector<double> &b);
double MBWR_dintegral_const_T_energy_dT(double rho, double T, double rho_c, double M, const std::vector<double> &b);
double MBWR_integral_const_T_entropy(double rho, double T, double rho_c, double M, const std::vector<double> &b);

#endif
