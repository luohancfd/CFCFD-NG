// Author: Daniel F Potter
// Date: 07-Dec-2009

#ifndef PARK_NONEQUILIBRIUM_HH
#define PARK_NONEQUILIBRIUM_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "reaction-rate-coeff.hh"
#include "../models/gas_data.hh"
#include "generalised-Arrhenius.hh"

class Park_nonequilibrium : public Generalised_Arrhenius {
public:
    Park_nonequilibrium(lua_State *L, Gas_model &g);
    Park_nonequilibrium(double A, double n, double E_a, double s_p, int iTp, int iTq);
    ~Park_nonequilibrium();

private:
    int s_eval(const Gas_data &Q);
    int iTp_;
    int iTq_;
    double s_p_;
};

Reaction_rate_coefficient* create_Park_nonequilibrium_coefficient(lua_State *L, Gas_model &g);

#endif
