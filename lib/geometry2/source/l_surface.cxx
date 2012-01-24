// Author: Rowan J. Gollan
// Date: 15-Jun-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#include <sstream>
#include <cstdlib>

#include "../../util/source/useful.h"
#include "l_geom.hh"
#include "l_surface.hh"

using namespace std;

luaSurface::
luaSurface()
    : surf_(0) {}

luaSurface::
~luaSurface()
{
    delete surf_;
}

int
luaSurface::
r0(lua_State *L)
{
    int narg = lua_gettop(L);

    if ( narg == 0 ) { // This is a getter
	lua_pushnumber(L, surf_->r0);
	return 1;
    }
    // else, a setter
    surf_->r0 = luaL_checknumber(L, 1);
    return 0;
}

int
luaSurface::
r1(lua_State *L)
{
    int narg = lua_gettop(L);

    if ( narg == 0 ) { // This is a getter
	lua_pushnumber(L, surf_->r1);
	return 1;
    }
    // else, a setter
    surf_->r1 = luaL_checknumber(L, 1);
    return 0;
}

int
luaSurface::
s0(lua_State *L)
{
    int narg = lua_gettop(L);

    if ( narg == 0 ) { // This is a getter
	lua_pushnumber(L, surf_->s0);
	return 1;
    }
    // else, a setter
    surf_->s0 = luaL_checknumber(L, 1);
    return 0;
}

int
luaSurface::
s1(lua_State *L)
{
    int narg = lua_gettop(L);

    if ( narg == 0 ) { // This is a getter
	lua_pushnumber(L, surf_->s1);
	return 1;
    }
    // else, a setter
    surf_->s1 = luaL_checknumber(L, 1);
    return 0;
}


int
luaSurface::
eval(lua_State *L)
{
    double r = luaL_checknumber(L, 1);
    double s = luaL_checknumber(L, 2);

    Vector3 p = surf_->eval(r, s);
    Lunar<luaVector3>::push(L, new luaVector3(p), true);
    return 1;
}

int
luaSurface::
dpdr(lua_State *L)
{
    double r = luaL_checknumber(L, 1);
    double s = luaL_checknumber(L, 2);

    Vector3 dpdr = surf_->dpdr(r, s);
    Lunar<luaVector3>::push(L, new luaVector3(dpdr), true);
    return 1;
}

int
luaSurface::
dpds(lua_State *L)
{
    double r = luaL_checknumber(L, 1);
    double s = luaL_checknumber(L, 2);

    Vector3 dpds = surf_->dpds(r, s);
    Lunar<luaVector3>::push(L, new luaVector3(L), true);
    return 1;
}

int
luaSurface::
translate(lua_State *L)
{
    int narg = lua_gettop(L);

    if ( narg == 1 ) {
	luaVector3* tmp = Lunar<luaVector3>::check(L, 1);
	surf_->translate(tmp->r_object());
    }
    else if ( narg == 3 ) {
	double vx = luaL_checknumber(L, 1);
	double vy = luaL_checknumber(L, 2);
	double vz = luaL_checknumber(L, 3);
	surf_->translate(vx, vy, vz);
    }
    else {
	ostringstream ost;
	ost << "error in call to translate():\n"
	    << "   expect 1 arg:  of type Vector3  <OR>\n"
	    << "   expect 3 args: each of type double.\n"
	    << "  " << narg << " argument(s) received.\n"
	    << "   Bailing out!\n";
	luaL_error(L, ost.str().c_str());
	exit(LUA_ERROR);
    }

    return 0;

}


int
luaSurface::
mirror_image(lua_State *L)
{
    luaVector3* point = Lunar<luaVector3>::check(L, 1);
    luaVector3* normal = Lunar<luaVector3>::check(L, 2); 
    surf_->mirror_image(point->r_object(), normal->r_object());

    return 0;
}

luaNurbsSurface::
luaNurbsSurface(lua_State *L)
    : luaSurface()
{
    // Accept a constructor of the form:
    // n = NurbsSurface(Q, w, p, U, q, V,
    //                  {label=.., r0=..., r1=...,
    //                   s0=..., s1=...})
    //
    luaL_argcheck(L, lua_istable(L, 1), 1,
		  "a table of 'Vector3' objects is expected");
    luaL_argcheck(L, lua_istable(L, 2), 2,
		  "a table of weights is expected");
    int p = luaL_checkint(L, 3);
    luaL_argcheck(L, p >= 1, 3,
		  "an integer >= 1 is expected");
    luaL_argcheck(L, lua_istable(L, 4), 4,
		  "a table of knots is expected");
    int q = luaL_checkint(L, 5);
    luaL_argcheck(L, q >= 1, 5,
		  "an integer >= 1 is expected");
    luaL_argcheck(L, lua_istable(L, 6), 6,
		  "a table of knots is expected");

    size_t n = lua_objlen(L, 1);
    vector<vector<Vector3> > Q(n);
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 1, i);
	size_t ncols = lua_objlen(L, -1);
	Q[i-1].resize(ncols);
	for ( size_t j = 1; j <= ncols; ++j ) {
	    lua_rawgeti(L, -1, j);
	    Q[i-1][j-1] = Lunar<luaVector3>::check(L, -1)->r_object();
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);
    }

    n = lua_objlen(L, 2);
    vector<vector<double> > w(n);
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 2, i);
	size_t ncols = lua_objlen(L, -1);
	w[i-1].resize(ncols);
	for ( size_t j = 1; j <= ncols; ++j ) {
	    lua_rawgeti(L, -1, j);
	    w[i-1][j-1] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);
    }

    n = lua_objlen(L, 4);
    vector<double> U;
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 4, i);
	U.push_back(luaL_checknumber(L, -1));
	lua_pop(L, 1);
    }

    n = lua_objlen(L, 6);
    vector<double> V;
    for ( size_t i = 1; i <= n; ++i ) {
	lua_rawgeti(L, 6, i);
	V.push_back(luaL_checknumber(L, -1));
	lua_pop(L, 1);
    }

    string label = "";
    double r0 = 0.0;
    double r1 = 1.0;
    double s0 = 0.0;
    double s1 = 1.0;

    if ( lua_istable(L, 7) ) {
	// optional table present.
	lua_getfield(L, 7, "label");
	if ( !lua_isnil(L, -1) ) {
	    label = luaL_checkstring(L, -1);
	}
	lua_pop(L, 1);

	lua_getfield(L, 7, "r0");
	if ( !lua_isnil(L, -1) ) {
	    r0 = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);

	lua_getfield(L, 7, "r1");
	if ( !lua_isnil(L, -1) ) {
	    r1 = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);

	lua_getfield(L, 7, "s0");
	if ( !lua_isnil(L, -1) ) {
	    s0 = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);

	lua_getfield(L, 7, "s1");
	if ( !lua_isnil(L, -1) ) {
	    s1 = luaL_checknumber(L, -1);
	}
	lua_pop(L, 1);
    }

    surf_ = new NurbsSurface(Q, w, p, U, q, V,
			     label, r0, r1, s0, s1);

}

luaNurbsSurface::
luaNurbsSurface(const luaNurbsSurface &n)
    : luaSurface()
{
    NurbsSurface* ns = dynamic_cast<NurbsSurface*>(n.surf_);
    surf_ = new NurbsSurface(*ns);
}

luaNurbsSurface::
luaNurbsSurface(const NurbsSurface &n)
    : luaSurface()
{
    surf_ = new NurbsSurface(n);
}

luaNurbsSurface::
~luaNurbsSurface() {}

string
luaNurbsSurface::
str() const
{
    ostringstream ost;
    ost << "luaNurbsSurface(... string representation not impelemented ...)\n";
    return ost.str();
}

luaNurbsSurface*
luaNurbsSurface::
clone() const
{
    return new luaNurbsSurface(*this);
}

//int
//luaNurbsSurface::
//write_IGES_file(lua_State *L)
//{
//    string fname = luaL_checkstring(L, 1);
//
//   NurbsSurface* ns = dynamic_cast<NurbsSurface*>(surf_);
//
//    ns->write_IGES_file(fname);
//
//    return 0;
//}

const char luaNurbsSurface::className[] = "NurbsSurface";

#define member_data(class, name) {#name, &class::name}

Lunar<luaNurbsSurface>::RegType luaNurbsSurface::member_data[] = {
    member_data(luaNurbsSurface, r0),
    member_data(luaNurbsSurface, r1),
    member_data(luaNurbsSurface, s0),
    member_data(luaNurbsSurface, s1),
    {0, 0}
};

#define method(class, name) {#name, &class::name}

Lunar<luaNurbsSurface>::RegType luaNurbsSurface::methods[] = {
    method(luaNurbsSurface, eval),
    method(luaNurbsSurface, dpdr),
    method(luaNurbsSurface, dpds),
    method(luaNurbsSurface, translate),
    method(luaNurbsSurface, mirror_image),
//    method(luaNurbsSurface, write_IGES_file),
    {0, 0}
};

Lunar<luaNurbsSurface>::MetaType luaNurbsSurface::metamethods[] = {
    {0, 0}
};

int open_surface(lua_State *L, int table)
{
    Lunar<luaNurbsSurface>::Register(L, table);
    return 0;
}
