// Author: Elise J. Fahy
// Date: 24-Sep-2014

#ifndef CONSTANT_PRANDTL_THERMAL_CONDUCTIVITY_HH
#define CONSTANT_PRANDTL_THERMAL_CONDUCTIVITY_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "thermal-conductivity-model.hh"

class constant_Prandtl_thermal_conductivity : public Thermal_conductivity_model {
public:
    constant_Prandtl_thermal_conductivity(lua_State *L);
    ~constant_Prandtl_thermal_conductivity();

private:
    double Pr_;
    
    double s_eval_thermal_conductivity(const Gas_data &Q);
};

double S_thermal_conductivity(double mu, double Cp, double Pr_);

#endif
