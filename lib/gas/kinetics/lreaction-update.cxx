// Author: Rowan J. Gollan
// Date: 27-Mar-2009
// Place: Poquoson, Virginia, USA

#include <sstream>
#include <string>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "lreaction-update.hh"
#include "../models/lgas-model.hh"
#include "../models/lservice_gas_data.hh"

using namespace std;

luaReaction_update::
luaReaction_update(lua_State *L)
{
    string cfile(luaL_checkstring(L, 1));
    luaGas_model *lg = Lunar<luaGas_model>::check(L, 2);
    r_ = create_Reaction_update(cfile, *(lg->r_pointer()));
    Q_ =  new Gas_data(lg->r_pointer());
    dfdt_.resize(Q_->massf.size());
}

luaReaction_update::
~luaReaction_update()
{
    delete r_;
    delete Q_;
}

const char luaReaction_update::className[] = "Reaction_update";

#define member_data(class, name) {#name, &class::name}

Lunar<luaReaction_update>::RegType luaReaction_update::member_data[] = {
    {0, 0}
};

#define method(class, name) {#name, &class::name}

Lunar<luaReaction_update>::RegType luaReaction_update::methods[] = {
    method(luaReaction_update, update_state),
    method(luaReaction_update, rate_of_change),
    {0, 0}
};

Lunar<luaReaction_update>::MetaType luaReaction_update::metamethods[] = {
    {0, 0}
};

string
luaReaction_update::
str() const
{
    return string("");
}

int
luaReaction_update::
update_state(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 2 and narg != 3) {
	ostringstream ost;
	ost << "Error in call to update_state():\n";
	ost << "2 or 3 arguments expected: gas_data table, dt and dt_suggest.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    double dt = luaL_checknumber(L, 2);
    double dt_suggest = -1.0;
    if ( narg == 3 ) 
	dt_suggest = luaL_checknumber(L, 3);

    // NOTE: passing a null pointer as Gas_model * gmodel (DFP 26/02/09)
    int flag = r_->update_state(*Q_, dt, dt_suggest);

    set_gas_data_at_table(L, 1, *Q_);
    
    lua_pushnumber(L, dt_suggest);

    if ( flag == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}

int
luaReaction_update::
rate_of_change(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 1 ) {
	ostringstream ost;
	ost << "Error in call to rate_of_change():\n";
	ost << "1 argument expected.\n";
	ost << narg << " argument(s) received.\n";
	input_error(ost);
    }

    lua_pushvalue(L, 1);
    get_table_as_gas_data(L, *Q_);
    lua_pop(L, 1);

    int flag = r_->rate_of_change(*Q_, dfdt_);

    push_vector_as_table(L, dfdt_);

    if ( flag == SUCCESS )
	lua_pushstring(L, "success");
    else
	lua_pushstring(L, "fail");

    return 2;
}
    
int open_reaction_update(lua_State *L, int table)
{
    Lunar<luaReaction_update>::Register(L, table);
    return 0;
}
