// Author: Rowan J. Gollan
// Date: 02-May-2009
// Place: Poquoson, Virginia, USA

#ifndef LREACTION_RATE_COEFF_HH
#define LREACTION_RATE_COEFF_HH

#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"
#include "../models/gas_data.hh"
#include "reaction-rate-coeff.hh"

class luaReaction_rate_coefficient {
public:
    static const char className[];
    static Lunar<luaReaction_rate_coefficient>::RegType member_data[];
    static Lunar<luaReaction_rate_coefficient>::RegType methods[];
    static Lunar<luaReaction_rate_coefficient>::MetaType metamethods[];

    std::string str() const;

    luaReaction_rate_coefficient(lua_State *L);
    ~luaReaction_rate_coefficient();

    int k(lua_State *L);
    int eval(lua_State *L);

private:
    Reaction_rate_coefficient *rrc_;
    Gas_data *Q_;
};

int open_reaction_rate_coefficient(lua_State *L, int table);

#endif


