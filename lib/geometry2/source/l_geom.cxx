#include <cstdlib>
#include <sstream>
#include <string>

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "../../util/source/useful.h"
#include "../../util/source/lunar.hh"
#include "geom.hh"

#include "l_geom.hh"
#include "l_gpath.hh"
#include "l_gpath_utils.hh"
#include "l_nurbs_utils.hh"
#include "l_surface.hh"

using namespace std;

string
luaVector3::
str() const
{
    return v_.str();
}

luaVector3::
luaVector3(lua_State *L)
    : v_(0.0, 0.0, 0.0)
{
    int narg = lua_gettop(L);

    if( narg >= 1 )
	v_.x = luaL_checknumber(L, 1);

    if( narg >= 2)
	v_.y = luaL_checknumber(L, 2);

    if( narg >= 3 )
	v_.z = luaL_checknumber(L, 3);
}

luaVector3::
luaVector3(const Vector3 &v)
    : v_(v) {}

luaVector3::
~luaVector3() {}

int
luaVector3::
x(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg == 0 ) {
	// This is a getter.
	lua_pushnumber(L, v_.x);
	return 1;
    }
    // else
    // Treat as a setter.
    v_.x = luaL_checknumber(L, 1);
    return 0;
}

int
luaVector3::
y(lua_State *L)
{
    //    cout << "luaVector3::y()  this= " << this << " &(v_)= " << &(v_) << endl;
    int narg = lua_gettop(L);

    if( narg == 0 ) {
	// This is a getter.
	lua_pushnumber(L, v_.y);
	return 1;
    }
    // else
    // Treat as a setter.
    v_.y = luaL_checknumber(L, 1);
    return 0;
}

int
luaVector3::
z(lua_State *L)
{
    int narg = lua_gettop(L);

    if( narg == 0 ) {
	// This is a getter.
	lua_pushnumber(L, v_.z);
	return 1;
    }
    // else
    // Treat as a setter.
    v_.z = luaL_checknumber(L, 1);
    return 0;
}

int
luaVector3::
norm(lua_State *L)
{
    v_.norm();
    return 0;
}

int
luaVector3::
copy(lua_State *L)
{
    Lunar<luaVector3>::push(L, new luaVector3(v_), true);
    return 1;
}

bool istype(lua_State *L, int index, const char *tname)
{
    void *p = lua_touserdata(L, index);
    if( p != NULL ) {
	if( lua_getmetatable(L, index) ) {
	    lua_getfield(L, LUA_REGISTRYINDEX, tname);
	    if( lua_rawequal(L, -1, -2) ) {
		lua_pop(L, 2);
		return true;
	    }
	    // Pop values anyway
	    lua_pop(L, 2);
	}
    }
    return false;
}

int l_equals(lua_State *L)
{
    const double tol = 1.0e-4;
    luaVector3 *a = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    luaVector3 *b = (luaVector3*) Lunar<luaVector3>::check(L, 2);
    lua_pushboolean(L, equal(a->r_object(), b->r_object(), tol));
    return 1;
}

int l_add(lua_State *L)
{
    luaVector3 *a = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    luaVector3 *b = (luaVector3*) Lunar<luaVector3>::check(L, 2);

    Vector3 c = a->r_object() + b->r_object();
    Lunar<luaVector3>::push(L, new luaVector3(c), true);
    return 1;
}

int l_subtract(lua_State *L)
{
    luaVector3 *a = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    luaVector3 *b = (luaVector3*) Lunar<luaVector3>::check(L, 2);

    Vector3 c = a->r_object() - b->r_object();
    Lunar<luaVector3>::push(L, new luaVector3(c), true);
    return 1;
}

int l_mul(lua_State *L)
{
    luaVector3 *a;
    double scalar;
    if( istype(L, 1, "Vector3") ) {
	a = (luaVector3*) Lunar<luaVector3>::check(L, 1);
	scalar = luaL_checknumber(L, 2);
    }
    else if( istype(L, 2, "Vector3") ) {
	a = (luaVector3*) Lunar<luaVector3>::check(L, 2);
	scalar = luaL_checknumber(L, 1);
    }
    else {
	ostringstream ost;
	ost << "error in call __mul():\n"
	    << "  2 arguments expected: 1 Vector3, 1 number.\n"
	    << "  Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }
	
    Vector3 b = scalar * a->r_object();
    Lunar<luaVector3>::push(L, new luaVector3(b), true);
    return 1;
}

int l_unm(lua_State *L)
{
    luaVector3 *a = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    Vector3 b = -(a->r_object());
    Lunar<luaVector3>::push(L, new luaVector3(b), true);
    return 1;
}

int l_vabs(lua_State *L)
{
    luaVector3 *l_v = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    Vector3 v = l_v->r_object();
    lua_pushnumber(L, vabs(v));
    return 1;
}

int l_unit(lua_State *L)
{
    luaVector3 *l_v = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    Vector3 v = l_v->r_object();
    Vector3 u = unit(v);
    Lunar<luaVector3>::push(L, new luaVector3(u), true);
    return 1;
}

int l_dot(lua_State *L)
{
    luaVector3 *a = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    luaVector3 *b = (luaVector3*) Lunar<luaVector3>::check(L, 2);
    lua_pushnumber(L, dot(a->r_object(), b->r_object()));
    return 1;
}

int l_cross(lua_State *L)
{
    luaVector3 *a = (luaVector3*) Lunar<luaVector3>::check(L, 1);
    luaVector3 *b = (luaVector3*) Lunar<luaVector3>::check(L, 2);
    Vector3 c = cross(a->r_object(), b->r_object());
    Lunar<luaVector3>::push(L, new luaVector3(c), true);
    return 1;
}

const char luaVector3::className[] = "Vector3";

#define member_data(class, name) {#name, &class::name}

Lunar<luaVector3>::RegType luaVector3::member_data[] = {
    member_data(luaVector3, x),
    member_data(luaVector3, y),
    member_data(luaVector3, z),
    {0, 0}
};

#define method(class, name) {#name, &class::name}

Lunar<luaVector3>::RegType luaVector3::methods[] = {
    //    method(luaVector3, x),
    //    method(luaVector3, y),
    //    method(luaVector3, z),
    method(luaVector3, norm),
    method(luaVector3, copy),
    {0, 0}
};

Lunar<luaVector3>::MetaType luaVector3::metamethods[] = {
    {"__eq", &l_equals},
    {"__add", &l_add},
    {"__sub", &l_subtract},
    {"__mul", &l_mul},
    {"__unm", &l_unm},
    {0, 0}
};

int open_geom(lua_State *L)
{
    Lunar<luaVector3>::Register(L);
    lua_pushcfunction(L, l_vabs);
    lua_setglobal(L, "vabs");
    lua_pushcfunction(L, l_unit);
    lua_setglobal(L, "unit");
    lua_pushcfunction(L, l_dot);
    lua_setglobal(L, "dot");
    lua_pushcfunction(L, l_cross);
    lua_setglobal(L, "cross");

    open_gpath(L);
    open_gpath_utils(L);
    return 0;
}

int luaopen_geometry(lua_State *L)
{
    // _G["geometry"] = {}
    lua_newtable(L);
    lua_setglobal(L, "geometry");
    
    // Bring "geometry" table to TOS
    lua_getglobal(L, "geometry");
    int table = lua_gettop(L);
    Lunar<luaVector3>::Register(L, table);

    // Set some general Vector3 functions
    lua_pushcfunction(L, l_vabs);
    lua_setfield(L, table, "vabs");
    lua_pushcfunction(L, l_unit);
    lua_setfield(L, table, "unit");
    lua_pushcfunction(L, l_dot);
    lua_setfield(L, table, "dot");
    lua_pushcfunction(L, l_cross);
    lua_setfield(L, table, "cross");

    // And the paths
    open_gpath(L, table);

    // And the surfaces
    open_surface(L, table);

    // And the path utilities
    open_gpath_utils(L, table);

    // And the special NURBS utilities
    open_nurbs_utils(L, table);

    return 1;
}
