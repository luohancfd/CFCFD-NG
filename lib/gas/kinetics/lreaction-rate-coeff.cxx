// Author: Rowan J. Gollan
// Date: 02-May-2009
// Place: Poquoson, Virginia, USA

#include <string>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "../models/gas_data.hh"
#include "lreaction-rate-coeff.hh"
#include "../models/lgas-model.hh"
#include "../models/lservice_gas_data.hh"

using namespace std;

luaReaction_rate_coefficient::
luaReaction_rate_coefficient(lua_State *L)
{
    luaGas_model *lg = Lunar<luaGas_model>::check(L, 2);
    lua_pushvalue(L, 1);
    rrc_ = create_Reaction_rate_coefficient(L, *(lg->r_pointer()));
    lua_pop(L, 1);
    Q_ = new Gas_data(lg->r_pointer());
}

luaReaction_rate_coefficient::
~luaReaction_rate_coefficient()
{
    delete rrc_;
    delete Q_;
}

const char luaReaction_rate_coefficient::className[] = "Reaction_rate_coefficient";

#define member_data(class, name) {#name, &class::name}

Lunar<luaReaction_rate_coefficient>::RegType luaReaction_rate_coefficient::member_data[] = {
    {0, 0}
};

#define method(class, name) {#name, &class::name}

Lunar<luaReaction_rate_coefficient>::RegType luaReaction_rate_coefficient::methods[] = {
    method(luaReaction_rate_coefficient, k),
    method(luaReaction_rate_coefficient, eval),
    {0, 0}
};

Lunar<luaReaction_rate_coefficient>::MetaType luaReaction_rate_coefficient::metamethods[] = {
    {0, 0}
};

string
luaReaction_rate_coefficient::
str() const
{
    return string("");
}

int
luaReaction_rate_coefficient::
k(lua_State *L)
{
    lua_pushnumber(L, rrc_->k());
    return 1;
}

int
luaReaction_rate_coefficient::
eval(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to eval():\n";
	ost << "1 argument expected.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int flag = rrc_->eval(*Q_);

    lua_pushinteger(L, flag);

    return 1;
}
    
int open_reaction_rate_coefficient(lua_State *L, int table)
{
    Lunar<luaReaction_rate_coefficient>::Register(L, table);
    return 0;
}
