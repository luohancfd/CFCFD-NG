// Author: Rowan J. Gollan
// Version: 17-Oct-2008
// Place: Hampton, Virginia, USA

#ifndef NORMAL_REACTION_HH
#define NORMAL_REACTION_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "reaction.hh"

class Normal_reaction : public Reaction {
public:
    Normal_reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower);
    virtual ~Normal_reaction();

protected:
    virtual double s_compute_forward_rate(const std::vector<double> &y);
    virtual double s_compute_backward_rate(const std::vector<double> &y);

private:
    std::map<int, int> f_coeffs_;
    std::map<int, int> b_coeffs_;
};

Reaction* create_Normal_reaction(lua_State *L, Gas_model &g, double T_upper, double T_lower);


#endif
