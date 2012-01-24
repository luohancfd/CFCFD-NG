// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#ifndef NONEQ_GAS_EOS_HH
#define NONEQ_GAS_EOS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "equation-of-state.hh"

class Noneq_gas : public Equation_of_state {
public:
    Noneq_gas(lua_State *L);
    ~Noneq_gas();
    
private:
    // heavy-particle / free electron indices
    int iT_;
    int iTe_;
    int iel_;
    
    // gas-constants
    std::vector<double> R_;
    
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
    double hp_gas_constant(const Gas_data &Q);
};


#endif
