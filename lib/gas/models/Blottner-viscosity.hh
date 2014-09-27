// Author: Elise J. Fahy
// Date: 24-Sep-2014

#ifndef BLOTTNER_VISCOSITY_HH
#define BLOTTNER_VISCOSITY_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "viscosity-model.hh"

class Blottner_viscosity : public Viscosity_model {
public:
    Blottner_viscosity(lua_State *L);
    ~Blottner_viscosity();

private:
    double A_mu_;
    double B_mu_;
    double C_mu_;

    double s_eval_viscosity(const Gas_data &Q);
};

#endif
