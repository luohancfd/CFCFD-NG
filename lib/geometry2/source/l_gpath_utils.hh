#ifndef L_GPATH_UTILS_HH
#define L_GPATH_UTILS_HH

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/lunar.hh"
#include "geom.hh"
#include "gpath.hh"
#include "gpath_utils.hh"

int l_create_line(lua_State *L);
int l_create_best_fit_XBezier(lua_State *L);
int l_create_optimised_Bezier(lua_State *L);
int l_create_optimised_Bezier_YZ(lua_State *L);
int l_create_optimised_Bezier3D(lua_State *L);
int l_create_best_fit_poly(lua_State *L);
int l_create_quintic_spline(lua_State *L);
int open_gpath_utils(lua_State *L, int table=LUA_GLOBALSINDEX);
#endif
