// Author: Rowan J. Gollan
// Date: 29-Mar-2009
// Place: Poquoson, Virginia, USA

#ifndef THIRD_BODY_REACTION_HH
#define THIRD_BODY_REACTION_HH
#include <valarray>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "normal-reaction.hh"

class Third_body_reaction : public Normal_reaction {
public:
    Third_body_reaction(lua_State *L, Gas_model &g);
    virtual ~Third_body_reaction();

private:
    double s_compute_forward_rate(const std::valarray<double> &y);
    double s_compute_backward_rate(const std::valarray<double> &y);
    void compute_third_body_concentration(const std::valarray<double> &y);

    std::map<int, double> efficiencies_;
    double conc_;
    bool conc_just_computed_;
};

Reaction* create_Third_body_reaction(lua_State *L, Gas_model &g);

#endif
