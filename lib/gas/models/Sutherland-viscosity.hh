// Author: Rowan J. Gollan
// Date: 10-Jul-2008

#ifndef SUTHERLAND_VISCOSITY_HH
#define SUTHERLAND_VISCOSITY_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "viscosity-model.hh"

class Sutherland_viscosity : public Viscosity_model {
public:
    Sutherland_viscosity(lua_State *L);
    ~Sutherland_viscosity();

private:
    double mu0_;
    double T0_;
    double S_;

    double s_eval_viscosity(const Gas_data &Q);
};

double S_viscosity(double T, double mu0, double T0, double S);

#endif
