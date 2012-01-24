// Author: Rowan J. Gollan
// Version: 24-May-2008
//            Initial coding.
//

#ifndef CONSTANT_SPECIFIC_HEATS_HH
#define CONSTANT_SPECIFIC_HEATS_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "thermal-behaviour-model.hh"

class Constant_specific_heats : public Thermal_behaviour_model {
public:
    Constant_specific_heats(lua_State *L);
    ~Constant_specific_heats();

private:
    std::vector<double> Cv_;
    std::vector<double> Cp_;
    std::vector<double> e_zero_;
    std::vector<double> q_;
    
    int s_decode_conserved_energy(Gas_data &Q, const std::vector<double> &rhoe);
    int s_encode_conserved_energy(const Gas_data &Q, std::vector<double> &rhoe);
    double s_dhdT_const_p(const Gas_data &Q, int &status);
    double s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status);
    int s_eval_energy(Gas_data &Q, Equation_of_state *EOS_);
    int s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_);
    double s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
};

#endif
