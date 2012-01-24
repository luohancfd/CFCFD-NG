// Author: Brendan T. O'Flaherty
//  based on perfect-gas-EOS code by Rowan J. Gollan
// Version: 14-May-2009
//            Initial coding.
//

#ifndef NOBLE_ABEL_GAS_EOS_HH
#define NOBLE_ABEL_GAS_EOS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "equation-of-state.hh"

class Noble_Abel_gas : public Equation_of_state {
public:
    Noble_Abel_gas(lua_State *L);
    ~Noble_Abel_gas();

    // python function
    double covolume(Gas_data Q) { int status; return s_covolume(Q, status); }

private:
    std::vector<double> R_;
    std::vector<double> nu_0_;

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
    double s_covolume(const Gas_data &Q, int &status);
};

double nag_pressure(double rho, double T, double R, double nu_0);
double nag_temperature(double rho, double p, double R, double nu_0);
double nag_density(double T, double p, double R, double nu_0);

#endif
