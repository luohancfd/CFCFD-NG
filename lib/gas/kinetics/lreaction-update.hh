// Author: Rowan J. Gollan
// Date: 27-Mar-2009
// Place: Poquoson, Virginia, USA

#ifndef LREACTION_UPDATE_HH
#define LREACTION_UPDATE_HH

#include <vector>
#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"
#include "../models/gas_data.hh"
#include "reaction-update.hh"

class luaReaction_update {
public:
    static const char className[];
    static Lunar<luaReaction_update>::RegType member_data[];
    static Lunar<luaReaction_update>::RegType methods[];
    static Lunar<luaReaction_update>::MetaType metamethods[];

    std::string str() const;

    luaReaction_update(lua_State *L);
    ~luaReaction_update();

    int update_state(lua_State *L);
    int rate_of_change(lua_State *L);

private:
    Reaction_update *r_;
    Gas_data *Q_;
    std::vector<double> dfdt_;
};

int open_reaction_update(lua_State *L, int table);

#endif


