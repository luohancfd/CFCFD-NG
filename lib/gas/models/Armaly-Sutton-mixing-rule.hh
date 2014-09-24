// Author: Elise J. Fahy
// Date: 24-Sep-2014

#ifndef ARMALY_SUTTON_MIXING_RULE_HH
#define ARMALY_SUTTON_MIXING_RULE_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "transport-coefficients-model.hh"
#include "viscosity-model.hh"
#include "thermal-conductivity-model.hh"

class Armaly_Sutton_mixing_rule : public Transport_coefficients_model {
public:
    Armaly_Sutton_mixing_rule(lua_State *L);
    ~Armaly_Sutton_mixing_rule();

private:
    std::vector<Viscosity_model*> VM_;

    double ignore_mole_fraction_;
    std::vector<double> M_;
    std::vector<double> x_;
    std::vector<double> mu_;
    matrix phi_;
    double A_;
    double B_;
    double F_;

    int s_eval_transport_coefficients(Gas_data &Q, double A_, double B_, double F_);
};

#endif
