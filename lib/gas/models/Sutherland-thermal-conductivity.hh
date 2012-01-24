// Author: Rowan J. Gollan
// Date: 10-Jul-2008

#ifndef SUTHERLAND_THERMAL_CONDUCTIVITY_HH
#define SUTHERLAND_THERMAL_CONDUCTIVITY_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "thermal-conductivity-model.hh"

class Sutherland_thermal_conductivity : public Thermal_conductivity_model {
public:
    Sutherland_thermal_conductivity(lua_State *L);
    ~Sutherland_thermal_conductivity();

private:
    double k0_;
    double T0_;
    double S_;
    
    double s_eval_thermal_conductivity(const Gas_data &Q);
};

double S_thermal_conductivity(double T, double k0, double T0, double S);

#endif
