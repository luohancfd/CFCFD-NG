// Author: Rowan J. Gollan
// Date: 26-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

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
#include "nurbs_utils.hh"
#include "l_nurbs_utils.hh"
#include "l_surface.hh"

using namespace std;

int l_quad_curve_interp(lua_State *L)
{
    int narg = lua_gettop(L);

    if ( narg < 4 ) {
	ostringstream ost;
	ost << "quad_curve_interp():\n";
	ost << "   Bad number of arguments, 4 expected but " << narg << " received.\n";
	ost << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    luaVector3* Q0 = Lunar<luaVector3>::check(L, 1);
    luaVector3* T0 = Lunar<luaVector3>::check(L, 2);
    luaVector3* Q1 = Lunar<luaVector3>::check(L, 3);
    luaVector3* T1 = Lunar<luaVector3>::check(L, 4);

    vector<Vector3> R;
    vector<double> w;

    int result = quad_curve_interp(Q0->r_object(), T0->r_object(),
				   Q1->r_object(), T1->r_object(),
				   R, w);

    lua_pushinteger(L, result);

    lua_newtable(L);
    for ( size_t i = 0; i < R.size(); ++i ) {
	Lunar<luaVector3>::push(L, new luaVector3(R[i]), true);
	lua_rawseti(L, -2, i+1);
    }

    lua_newtable(L);
    for ( size_t i = 0; i < w.size(); ++i ) {
	lua_pushnumber(L, w[i]);
	lua_rawseti(L, -2, i+1);
    }

    return 3;
}

int l_quad_curve_approx(lua_State *L)
{
    // Accept a call from Lua of the form:
    // quad_curve_approx(Q, E [, Kmax, pc])
    // where:
    // Q     -- collection of points to approximate
    // E     -- tolerance for the fit
    // Kmax  -- upper limit on segment size (optional)
    // pc    -- Boolean indicating if corners should be
    //          preserved or not

    int narg = lua_gettop(L);
    vector<Vector3> Q;
    luaVector3* P;
    double E;
    int Kmax = -1;
    bool pc = false;

    if ( narg <= 1 ) {
	ostringstream ost;
	ost << "error in call quad_curve_approx():\n";
	ost << "   2, 3 or 4 arguments expected.\n";
	ost << "   Number of arguments received: " << narg << endl;
	ost << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
	
    if ( narg >= 2 ) {
	if ( ! lua_istable(L, 1) ) {
	    ostringstream ost;
	    ost << "error in call quad_curve_approx():\n";
	    ost << "   A table of points is expected as the first argument.\n";
	    ost << "Bailing out!\n";
	    luaL_error(L, ost.str().c_str());
	    exit(LUA_ERROR);
	}

	for ( size_t i = 1; i <= lua_objlen(L, 1); ++i ) {
	    lua_rawgeti(L, 1, i);
	    P = Lunar<luaVector3>::check(L, -1);
	    Q.push_back(P->r_object());
	    lua_pop(L, 1);
	}

	E = luaL_checknumber(L, 2);
    }

    if ( narg >= 3 ) {
	Kmax = luaL_checkint(L, 3);
    }

    if ( narg >= 4 ) {
	pc = (bool) lua_toboolean(L, 4);
    }


    Nurbs n = quad_curve_approx(Q, E, Kmax, pc);

    Lunar<luaNurbs>::push(L, new luaNurbs(n), true);

    return 1;

}

int l_make_one_arc(lua_State *L)
{
    // Accept a call from Lua of the form:
    // make_one_arc(P0, T0, P2, T2, P)
    // where:
    // P0    -- starting point on arc
    // T0    -- tangent at starting point
    // P2    -- finishing point on arc
    // T2    -- tangent at finishing point
    // P     -- an interior point on arc
    //
    // Return --
    // P1    -- the mid control point for arc
    // w1    -- weight at mid control point
    //

    int narg = lua_gettop(L);

    if ( narg < 5 ) {
	ostringstream ost;
	ost << "error in call make_one_arc():\n";
	ost << "   5 arguments are expected.\n";
	ost << "   Number of arguments received: " << narg << endl;
	ost << "Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    luaVector3 *P0 = Lunar<luaVector3>::check(L, 1);
    luaVector3 *T0 = Lunar<luaVector3>::check(L, 2);
    luaVector3 *P2 = Lunar<luaVector3>::check(L, 3);
    luaVector3 *T2 = Lunar<luaVector3>::check(L, 4);
    luaVector3 *P = Lunar<luaVector3>::check(L, 5);

    Vector3 P1;
    double w1;

    int result = make_one_arc(P0->r_object(), T0->r_object(),
			      P2->r_object(), T2->r_object(),
			      P->r_object(), P1, w1);

    string result_str;
    if ( result == SUCCESS )
	result_str = "success";
    else
	result_str = "fail";

    Lunar<luaVector3>::push(L, new luaVector3(P1));
    lua_pushnumber(L, w1);
    lua_pushstring(L, result_str.c_str());

    return 3;

}

int l_skinned_surface(lua_State *L)
{
    luaL_argcheck(L, lua_istable(L, 1), 1, 
		  "a table of 'Nurbs' is expected");
    int q = luaL_checkint(L, 2);
    luaL_argcheck(L, q >= 1, 2, 
		  "an integere >= 1 is expected");

    size_t n = lua_objlen(L, 1);
    vector<Nurbs> C;
    Nurbs *n_pointer;
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 1, i);
	n_pointer = dynamic_cast<Nurbs*>(Lunar<luaNurbs>::check(L, -1)->r_pointer());
	C.push_back(*n_pointer);
	lua_pop(L, 1);
    }

    NurbsSurface ns = skinned_surface(C, q);

    Lunar<luaNurbsSurface>::push(L, new luaNurbsSurface(ns), true);
    return 1;
}

int open_nurbs_utils(lua_State *L, int table)
{
    lua_pushcfunction(L, l_quad_curve_interp);
    lua_setfield(L, table, "quad_curve_interp");
    lua_pushcfunction(L, l_quad_curve_approx);
    lua_setfield(L, table, "quad_curve_approx");
    lua_pushcfunction(L, l_make_one_arc);
    lua_setfield(L, table, "make_one_arc");
    lua_pushcfunction(L, l_skinned_surface);
    lua_setfield(L, table, "skinned_surface");

    return 0;
}
