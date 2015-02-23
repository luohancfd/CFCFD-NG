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
import luageom;

immutable string LineMT = "Line"; // Name of Line metatable

extern(C) int opCallPath(T, string MTname)(lua_State* L)
{
    auto path = checkObj!(T, MTname)(L, 1);
    auto t = luaL_checknumber(L, 2);
    auto pt = path(t);
    return pushVector3(L, pt);
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
    return pushObj!(Line, LineMT)(L, ab);
    // We could think about deleting line ab,
    // or just rely on the garbage collector
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
    lua_pushcfunction(L, &toStringObj!(Line, LineMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, LineMT.toStringz);

}
    






