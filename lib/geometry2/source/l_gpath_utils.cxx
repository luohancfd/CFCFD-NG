#include <cstdlib>
#include <iostream>
#include <sstream>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/useful.h"
#include "../../util/source/lunar.hh"
#include "geom.hh"
#include "l_geom.hh"
#include "gpath.hh"
#include "l_gpath.hh"
#include "gpath_utils.hh"
#include "l_gpath_utils.hh"

using namespace std;

int l_create_line(lua_State *L)
{
    // Accept a call from Lua of the form:
    // create_line(P, m, x0, x1)
    // where:
    // P     -- point (lua_Vector3)
    // m     -- gradient (number)
    // x0    -- start of line (number)
    // x1    -- end of line (number)
    
    int narg = lua_gettop(L);

    if( narg != 4 ) {
	ostringstream ost;
	ost << "error in call create_line():\n"
	    << "   4 arguments expected. " << narg << " argument(s) received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    luaVector3* P = Lunar<luaVector3>::check(L, 1);
    double m = luaL_checknumber(L, 2);
    double x0 = luaL_checknumber(L, 3);
    double x1 = luaL_checknumber(L, 4);

    XPoly a = create_line(P->r_object(), m, x0, x1);
    Lunar<luaXPoly>::push(L, new luaXPoly(a), true);
    return 1;
}

int l_create_best_fit_XBezier(lua_State *L)
{
    // Accept call from Lua of the form:
    // create_best_fit_XBezier(P, n, lslope, rslope,
    //                         {lcurv=0.0, rcurv=0.0})
    // where:
    // P      -- points (Lua array)
    // n      -- order of Bezier (int)
    // lslope -- slope at left end (number)
    // rslope -- slope at right end (number)
    // Optional arguments:
    // lcurv  -- curvature at left end (number) [default = 0.0]
    // rcurv  -- curvature at right end (number) [default = 0.0]

    int narg = lua_gettop(L);

    if( narg != 4 && narg != 5 ) {
	ostringstream ost;
	ost << "error in call create_best_fit_XBezier():\n"
	    << "   4 or 5 arguments expected. " << narg << " argument(s) received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
		ostringstream ost;
	ost << "error in call create_best_fit_XBezier():\n"
	    << "   Table of Vector3 objects expected.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    int nVectors = lua_objlen(L, 1);
    for( int i = 1; i <= nVectors; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    int n = luaL_checkinteger(L, 2);
    double lslope = luaL_checknumber(L, 3);
    double rslope = luaL_checknumber(L, 4);
    
    double lcurv = 0.0;
    double rcurv = 0.0;

    if( narg == 5 ) {
	if( !lua_istable(L, 5) ) {
		ostringstream ost;
	ost << "error in call create_best_fit_XBezier():\n"
	    << "   Table with optional arguments expected.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
	}
	// Extract optional arguments
	lua_getfield(L, 5, "lcurv");
	if( lua_isnumber(L, -1) ) {
	    lcurv = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);

	lua_getfield(L, 5, "rcurv");
	if( lua_isnumber(L, -1) ) {
	    rcurv = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);
    }
    
    XBezier xbez = create_best_fit_XBezier(tmp, n, lslope, rslope,
					   lcurv, rcurv);
    Lunar<luaXBezier>::push(L, new luaXBezier(xbez), true);
    return 1;
}

int l_create_optimised_Bezier(lua_State *L)
{
    // Accept call from Lua of the form:
    // create_optimised_Bezier(P, n, lslope, rslope)
    // where:
    // P      -- points (Lua array)
    // n      -- order of Bezier (int)
    // dydxA -- slope at start (number)
    // dydxB -- slope at end (number)

    int narg = lua_gettop(L);

    if ( narg != 4 && narg != 5) {
	ostringstream ost;
	ost << "error in call create_best_fit_Bezier():\n"
	    << "   4 or 5 arguments expected. " << narg << " argument(s) received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
		ostringstream ost;
	ost << "error in call create_best_fit_Bezier():\n"
	    << "   Table of Vector3 objects expected.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    int nVectors = lua_objlen(L, 1);
    for( int i = 1; i <= nVectors; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    int n = luaL_checkinteger(L, 2);
    double dydxA = luaL_checknumber(L, 3);
    double dydxB = luaL_checknumber(L, 4);
    double z_pos = 0.0;
    if ( narg == 5 ) {
	z_pos = luaL_checknumber(L, 5);
    }
    Bezier bez = create_optimised_Bezier(tmp, n, dydxA, dydxB, z_pos);
    Lunar<luaBezier>::push(L, new luaBezier(bez), true);
    return 1;
}

int l_create_optimised_Bezier_YZ(lua_State *L)
{
    // Accept call from Lua of the form:
    // create_optimised_Bezier_YZ(P, n, lslope, rslope)
    // where:
    // P      -- points (Lua array)
    // n      -- order of Bezier (int)
    // dydzA -- slope at start (number)
    // dydzB -- slope at end (number)

    int narg = lua_gettop(L);

    if ( narg != 4 && narg != 5) {
	ostringstream ost;
	ost << "error in call create_best_fit_Bezier():\n"
	    << "   4 or 5 arguments expected. " << narg << " argument(s) received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
		ostringstream ost;
	ost << "error in call create_best_fit_Bezier():\n"
	    << "   Table of Vector3 objects expected.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    int nVectors = lua_objlen(L, 1);
    for( int i = 1; i <= nVectors; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    int n = luaL_checkinteger(L, 2);
    double dydzA = luaL_checknumber(L, 3);
    double dydzB = luaL_checknumber(L, 4);
    double x_pos = 0.0;
    if ( narg == 5 ) {
	x_pos = luaL_checknumber(L, 5);
    }
    Bezier bez = create_optimised_Bezier_YZ(tmp, n, dydzA, dydzB, x_pos);
    Lunar<luaBezier>::push(L, new luaBezier(bez), true);
    return 1;
}


int l_create_optimised_Bezier3D(lua_State *L)
{
    // Accept call from Lua of the form:
    // create_optimised_Bezier3D(P, n, dydxA, dydxB, dzdxA, dzdxB)
    // where:
    // P      -- points (Lua array)
    // n      -- order of Bezier (int)
    // dydxA -- dydx slope at start (number)
    // dydxB -- dydx slope at end (number)
    // dzdxA -- dzdx slope at start
    // dzdxB -- dzdx slope at end

    int narg = lua_gettop(L);

    if ( narg != 6 ) {
	ostringstream ost;
	ost << "error in call create_best_fit_Bezier3D():\n"
	    << "   6 arguments expected. " << narg << " argument(s) received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
		ostringstream ost;
	ost << "error in call create_best_fit_Bezier():\n"
	    << "   Table of Vector3 objects expected.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    int nVectors = lua_objlen(L, 1);
    for( int i = 1; i <= nVectors; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    int n = luaL_checkinteger(L, 2);
    double dydxA = luaL_checknumber(L, 3);
    double dydxB = luaL_checknumber(L, 4);
    double dzdxA = luaL_checknumber(L, 5);
    double dzdxB = luaL_checknumber(L, 6);
    
    Bezier bez = create_optimised_Bezier3D(tmp, n, dydxA, dydxB, dzdxA, dzdxB);
    Lunar<luaBezier>::push(L, new luaBezier(bez), true);
    return 1;
}

int l_create_best_fit_poly(lua_State *L)
{
    // Accept call from Lua of the form:
    // create_best_fit_poly(P, n)
    // where:
    // P      -- points (Lua array)
    // n      -- order of polynomial (int)

    int narg = lua_gettop(L);

    if( narg != 2 ) {
	ostringstream ost;
	ost << "error in call create_best_fit_poly():\n"
	    << "   2 arguments expected. " << narg << " argument(s) received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
		ostringstream ost;
	ost << "error in call create_best_fit_poly():\n"
	    << "   Table of Vector3 objects expected.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    int nVectors = lua_objlen(L, 1);
    for( int i = 1; i <= nVectors; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    int n = luaL_checkinteger(L, 2);
    
    XPoly xpoly = create_best_fit_poly(tmp, n);
    Lunar<luaXPoly>::push(L, new luaXPoly(xpoly), true);
    return 1;
}

int l_create_quintic_spline(lua_State *L)
{
    // Accept call from Lua of the form:
    // create_quintic_spline(points, slopes,
    //                      {lcurv=0.0, rcurv=0.0})
    // where:
    // points -- points (Lua array of Vector3)
    // slopes -- slopes (lua array of numbers)
    // Optional arguments:
    // lcurv  -- curvature at left end (number) [default = 0.0]
    // rcurv  -- curvature at right end (number) [default = 0.0]

    int narg = lua_gettop(L);

    if( narg != 2 && narg != 3 ) {
	ostringstream ost;
	ost << "error in call create_quintic_spline():\n"
	    << "   2 or 3 arguments expected. " << narg << " argument(s) received.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    vector<Vector3> tmp;
    if( !lua_istable(L, 1) ) {
		ostringstream ost;
	ost << "error in call create_quintic_spline():\n"
	    << "   Table of Vector3 objects expected as argument 1.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    int nVectors = lua_objlen(L, 1);
    for( int i = 1; i <= nVectors; ++i ) {
	lua_rawgeti(L, 1, i);
	luaVector3* v = Lunar<luaVector3>::check(L, -1);
	tmp.push_back(v->r_object());
	lua_pop(L, 1);
    }

    vector<double> tmp2;
    if( !lua_istable(L, 2) ) {
		ostringstream ost;
	ost << "error in call create_quintic_spline():\n"
	    << "   Table of numbers objects expected as argument 2.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    int ndoubles = lua_objlen(L, 2);
    for( int i = 1; i <= ndoubles; ++i ) {
	lua_rawgeti(L, 2, i);
	tmp2.push_back(luaL_checknumber(L, -1));
	lua_pop(L, 1);
    }

    double lcurv = 0.0;
    double rcurv = 0.0;

    if( narg == 3 ) {
	if( !lua_istable(L, 3) ) {
		ostringstream ost;
	ost << "error in call create_quintic_spline():\n"
	    << "   Table with optional arguments expected.\n"
	    << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
	}
	// Extract optional arguments
	lua_getfield(L, 3, "lcurv");
	if( lua_isnumber(L, -1) ) {
	    lcurv = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);

	lua_getfield(L, 3, "rcurv");
	if( lua_isnumber(L, -1) ) {
	    rcurv = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);
    }
    
    XSpline xspline = create_quintic_spline(tmp, tmp2, lcurv, rcurv);

    Lunar<luaXSpline>::push(L, new luaXSpline(xspline), true);
    return 1;
}

int open_gpath_utils(lua_State *L, int table)
{
    lua_pushcfunction(L, l_create_line);
    lua_setfield(L, table, "create_line");
    lua_pushcfunction(L, l_create_best_fit_XBezier);
    lua_setfield(L, table, "create_best_fit_XBezier");
    lua_pushcfunction(L, l_create_optimised_Bezier);
    lua_setfield(L, table, "create_optimised_Bezier");
    lua_pushcfunction(L, l_create_optimised_Bezier_YZ);
    lua_setfield(L, table, "create_optimised_Bezier_YZ");
    lua_pushcfunction(L, l_create_optimised_Bezier3D);
    lua_setfield(L, table, "create_optimised_Bezier3D");
    lua_pushcfunction(L, l_create_best_fit_poly);
    lua_setfield(L, table, "create_best_fit_poly");
    lua_pushcfunction(L, l_create_quintic_spline);
    lua_setfield(L, table, "create_quintic_spline");

    return 0;
}
