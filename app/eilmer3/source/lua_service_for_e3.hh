// lua_service_for_e3.hh
// Author: Rowan J. Gollan
// Date: 13-Mar-2013
// Place: The University of Queensland
//
// History: 13-Mar-2013
//          Refactored some service functions from bc_user_defined.hh
//

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

int luafn_sample_flow(lua_State *L);
int luafn_locate_cell(lua_State *L);


int register_luafns(lua_State *L);
