// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding.
//

#ifndef PERFECT_GAS_EOS_HH
#define PERFECT_GAS_EOS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "equation-of-state.hh"

class Perfect_gas : public Equation_of_state {
public:
    Perfect_gas(lua_State *L);
    ~Perfect_gas();
    
private:
    std::vector<double> R_;
    int iel_;	// electron isp index
    
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
};

double pg_pressure(double rho, double T, double R);
double pg_temperature(double rho, double p, double R);
double pg_density(double T, double p, double R);

#endif
