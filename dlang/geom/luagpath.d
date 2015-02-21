/**
 * A Lua interface for the D gpath module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-22
 */

module luagpath;

// We cheat to get the C Lua headers by using LuaD.
import luad.all;
import luad.c.lua;
import luad.c.lauxlib;
import std.stdio;
import std.string;
import geom;
import gpath;
import luageom2;

immutable string LineMT = "Line"; // Name of Line metatable

/* ----------------  Parameterized functions -----------------
 * These should work for all path types.
 * ------------------------------------------------------------

/**
 * This creates a new userdata spot on the stack and
 * populates it with a Path object.
 *
 * Note: we explicitly set the metatable name as
 * a string and do NOT try to get the type name
 * using T.stringof (or a more elaborate __traits function).
 * The reason is that we want control of the name that
 * appears in the Lua script, no matter what name the D
 * compiler decides to give you class type.
 */
int pushPath(T, string MTname)(lua_State *L, in T path)
{
    auto pPtr = cast(T*) lua_newuserdata(L, path.sizeof);
    *pPtr = path.dup();
    luaL_getmetatable(L, MTname.toStringz);
    lua_setmetatable(L, -2);
    return 1;
}

/**
 * Provides a sanity check that the raw userdata
 * is in fact what we think it is.
 */
T* checkPath(T, string MTname)(lua_State *L, int index)
{
    auto pPtr = cast(T*) luaL_checkudata(L, index, MTname.toStringz);
    return pPtr;
}

extern(C) int opCallPath(T, string MTname)(lua_State* L)
{
    auto path = checkPath!(T, MTname)(L, 1);
    auto t = luaL_checknumber(L, 2);
    auto pt = (*path)(t);
    return pushVector3(L, pt);
}


extern(C) int toStringPath(T, string MTname)(lua_State* L)
{
    auto path = checkPath!(T, MTname)(L, 1);
    lua_pushstring(L, path.toString.toStringz);
    return 1;
}

/* ----------------- Line specific functions --------------- */

/**
 * The Lua constructor for a Line.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{}
 * b = Vector3:new{1, 1}
 * ab = Line:new(a, b)
 * ---------------------------------
 */
extern(C) int newLine(lua_State* L)
{
    lua_remove(L, 1); // remove first agurment "this"
    int narg = lua_gettop(L);
    if ( narg != 2 ) {
	string errMsg = "Error in number of arguments to Line:new()\n";
	errMsg ~= format("2 arguments expected, but %s recevied.\n", narg);
	luaL_error(L, errMsg.toStringz);
    }
    // else, just assume we can get two good Vector3s.
    auto a = checkVector3(L, 1);
    auto b = checkVector3(L, 2);
    auto ab = new Line(*a, *b);
    return pushPath!(Line, LineMT)(L, ab);
}

void registerPaths(LuaState lua)
{
    auto L = lua.state;

    // Register the Line object
    luaL_newmetatable(L, LineMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newLine);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Line, LineMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Line, LineMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringPath!(Line, LineMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, LineMT.toStringz);

}
    






