// Author: Rowan J. Gollan
// Date: 26-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#ifndef L_NURBS_UTILS_HH
#define L_NURBS_UTILS_HH

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"

int l_quad_curve_interp(lua_State *L);
int l_quad_curve_approx(lua_State *L);
int l_make_one_arc(lua_State *L);
int l_skinned_surface(lua_State *L);
int open_nurbs_utils(lua_State *L, int table=LUA_GLOBALSINDEX);

#endif

