// Author: Kan Qin
// Based on perfect gas EOS of Rowan Gollan
// Version: 29-Apirl-2014
// A simple EOS where pressure is only the function of density
// used to create a very dense air

#ifndef SIMPLE_GAS_EOS_HH
#define SIMPLE_GAS_EOS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "equation-of-state.hh"

class Simple_gas : public Equation_of_state {
public:
    Simple_gas(lua_State *L);
    ~Simple_gas();
    
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

double simple_pressure(double rho);
double simple_temperature(double rho, double p, double R);
double simple_density(double p);

#endif
