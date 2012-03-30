// Author: Brendan T. O'Flaherty
//  based on perfect-thermal-behaviour.hh
// Date: 10-Aug-2009
// Place: Brisbane, Qld., A.

#ifndef REAL_THERMAL_BEHAVIOUR_HH
#define REAL_THERMAL_BEHAVIOUR_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/segmented-functor.hh"
#include "gas_data.hh"
#include "thermal-behaviour-model.hh"

class Real_thermal_behaviour : public Thermal_behaviour_model {
public:
    Real_thermal_behaviour(lua_State *L);
    ~Real_thermal_behaviour();

private:
    double T_COLD_;
    std::vector<double> M_;
    std::vector<double> R_;
    std::vector<Segmented_functor *> Cp_;
    std::vector<Segmented_functor *> h_;
    std::vector<Segmented_functor *> s_;

    int s_decode_conserved_energy(Gas_data &Q, const std::vector<double> &rhoe);
    int s_encode_conserved_energy(const Gas_data &Q, std::vector<double> &rhoe);
    int s_eval_energy(Gas_data &Q, Equation_of_state *EOS_);
    int s_eval_temperature(Gas_data &Q, Equation_of_state *EOS_);
    double s_dhdT_const_p(const Gas_data &Q, Equation_of_state *EOS_, int &status);
    double s_dedT_const_v(const Gas_data &Q, Equation_of_state *EOS_, int &status);
    double s_eval_energy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_enthalpy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double s_eval_entropy_isp(const Gas_data &Q, Equation_of_state *EOS_, int isp);
    double zero_function(Gas_data &Q, Equation_of_state *EOS_, double e_given, double T);
    double deriv_function(Gas_data &Q, Equation_of_state *EOS_, double T);
    bool test_T_for_polynomial_breaks(double T);
};


#endif
