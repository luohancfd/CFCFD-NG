/**
 * An Lua interface for the D geom module.
 *
 * This follows and adapts the examples given
 * in PIL in Chapter 28.
 *
 * Reference:
 * Ierusalimschy, R. (2006)
 * Programming in Lua, 2nd Edition
 * Lua.org, Rio de Janeiro 
 *
 * Author: Rowan G. and Peter J.
 * Date: 21-Feb-2014
 */

module luageom2;

// We cheat to get the C Lua headers by using LuaD.
import luad.all;
import luad.c.lua;
import luad.c.lauxlib;
import std.stdio;
import geom;

immutable string Vector3MT = "Vector3"; // Name of Vector3 metatable

/**
 * This creates a new userdata spot on the stack
 * and populates it with a Vector3 struct.
 */
int pushVector3(lua_State *L, Vector3 vec)
{
    auto vPtr = cast(Vector3*) lua_newuserdata(L, vec.sizeof);
    *vPtr = vec;
    luaL_getmetatable(L, Vector3MT.toStringz);
    lua_setmetatable(L, -2);
    return 1;
}

/**
 * This function will serve as our "constructor"
 * in the Lua script.
 *
 * We are free to support any type of construction.
 * Multiple construction types are also possible.
 * It's just a question of how much flexibility
 * we offer the user.
 */
extern(C) int newVector3(lua_State *L)
{
    auto vec = Vector3(0.0, 0.0, 0.0);
    /* This is where we decide how the user will instantiate
     * an object in Lua-land. Let's support positional argument
     * construction like our old Python interface for the moment.
     */
    lua_remove(L, 1); // remove self.
    int narg = lua_gettop(L);
    if ( narg >= 1 ) vec.refx = luaL_checknumber(L, 1);
    if ( narg >= 2 ) vec.refy = luaL_checknumber(L, 2);
    if ( narg >= 3 ) vec.refz = luaL_checknumber(L, 3);
    
    /* Regardless of how we filled in vec. We are now
     * ready to grab a piece of the lua stack and
     * place our new Vector3 there as userdata.
     */
    return pushVector3(L, vec);
}

/**
 * Provides a sanity check that the raw userdata
 * is in fact what we think it is.
 */
Vector3* checkVector3(lua_State *L, int index)
{
    auto vPtr = cast(Vector3*) luaL_checkudata(L, index, Vector3MT.toStringz);
    return vPtr;
}

/*-------- exposed Vector3 methods ------------ */

// The x(), y(), z() methods are a little funny.
// We are faking data access in a sense.
extern(C) int xVector3(lua_State* L)
{
    int narg = lua_gettop(L);
    auto a = checkVector3(L, 1);
    if( narg == 1 ) { // This is a getter
	lua_pushnumber(L, a.x);
	return 1;
    }
    // else: treat as a setter.
    a.refx = luaL_checknumber(L, 2);
    return 0;
}

/**
 * Normalizes a Vector3 object. Exposes geom.Vector3.normalize()
 */
extern(C) int normalizeVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    a.normalize();
    return 0;
}

extern(C) int toStringVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    lua_pushstring(L, a.toString.toStringz);
    return 1;
}

/*--------- exposed Vector3 overloaded operators */

/**
 * Adds two Vector3 objects. Exposes geom.Vector3.opBinary("+")
 */
extern(C) int addVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto b = checkVector3(L, 2);
    auto c = (*a) + (*b);
    return pushVector3(L, c);
}

void registerVector3(LuaState lua)
{
    auto L = lua.state;
    luaL_newmetatable(L, Vector3MT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newVector3);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &xVector3);
    lua_setfield(L, -2, "x");
    lua_pushcfunction(L, &normalizeVector3);
    lua_setfield(L, -2, "normalize");
    lua_pushcfunction(L, &addVector3);
    lua_setfield(L, -2, "__add");
    lua_pushcfunction(L, &toStringVector3);
    lua_setfield(L, -2, "__tostring");
    lua_setglobal(L, Vector3MT.toStringz);

    /* Also attempt to put "add" in the global namespace. */
    lua_pushcfunction(L, &addVector3);
    lua_setglobal(L, "add");
}
    
