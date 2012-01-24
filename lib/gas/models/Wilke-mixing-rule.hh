// Author: Rowan J. Gollan
// Date: 10-Jul-2008

#ifndef WILKE_MIXING_RULE_HH
#define WILKE_MIXING_RULE_HH

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

class Wilke_mixing_rule : public Transport_coefficients_model {
public:
    Wilke_mixing_rule(lua_State *L);
    ~Wilke_mixing_rule();

private:
    std::vector<Viscosity_model*> VM_;
    std::vector<Thermal_conductivity_model*> TCM_;

    double ignore_mole_fraction_;
    std::vector<double> M_;
    std::vector<double> x_;
    std::vector<double> mu_;
    std::vector<double> k_;
    matrix phi_;
    matrix psi_;

    int s_eval_transport_coefficients(Gas_data &Q);
};

#endif
