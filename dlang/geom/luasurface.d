/**
 * luasurface.d
 * Lua interface to Surface objects.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-24
 */

module luasurface;

import std.stdio;
import std.string;
import luad.all;
import luad.c.lua;
import luad.c.lauxlib;
import util.lua_service;
import geom;
import gpath;
import surface;
import luageom;
import luagpath;

/// Name of CoonsPatch metatable -- this is the Lua access name.
immutable string CoonsPatchMT = "CoonsPatch";
/// Name of AOPatch metatable -- this is the Lua access name.
immutable string AOPatchMT = "AOPatch"; // Name of 


extern(C) int opCallSurface(T, string MTname)(lua_State* L)
{
    auto surface = checkObj!(T, MTname)(L, 1);
    auto r = luaL_checknumber(L, 2);
    auto s = luaL_checknumber(L, 3);
    auto pt = surface(r, s);
    return pushVector3(L, pt);
}

void getPaths(lua_State *L, string ctorName, out Path[string] paths)
{
    // Assume that table is at index 1.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the '%s' path was not of Path type.`; 
    int index = 1;
    int nstack = lua_gettop(L);
    lua_getfield(L, index, "north");
    paths["north"] = checkPath(L, -1);
    if ( paths["north"] is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "north")));
    lua_pop(L, 1);
    lua_getfield(L, index, "east");
    paths["east"] = checkPath(L, -1);
    if ( paths["east"] is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "east")));
    lua_pop(L, 1);
    lua_getfield(L, index, "south");
    paths["south"] = checkPath(L, -1);
    if ( paths["south"] is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "south")));
    lua_pop(L, 1);
    lua_getfield(L, index, "west");
    paths["west"] = checkPath(L, -1);
    if ( paths["west"] is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "west")));
    lua_pop(L, 1);
}

void getVector3s(lua_State *L, string ctorName,
		 out Vector3 p00, Vector3 p10, out Vector3 p11, out Vector3 p01)
{
    // Assume that table is at index 1.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the '%s' corner was not of Vector3 type.`; 
    int index = 1;
    lua_getfield(L, index, "p00");
    auto p00Ptr = checkVector3(L, -1);
    if ( p00Ptr is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "p00")));
    else p00 = *p00Ptr;
    lua_pop(L, 1);
    lua_getfield(L, index, "p10");
    auto p10Ptr = checkVector3(L, -1);
    if ( p10Ptr is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "p10")));
    else p10 = *p10Ptr;
    lua_pop(L, 1);
    lua_getfield(L, index, "p11");
    auto p11Ptr = checkVector3(L, -1);
    if ( p11Ptr is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "p11")));
    else p11 = *p11Ptr;
    lua_pop(L, 1);
    lua_getfield(L, index, "p01");
    auto p01Ptr = checkVector3(L, -1);
    if ( p01Ptr is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "p01")));
    else p01 = *p01Ptr;
    lua_pop(L, 1);
}

void getRandS(lua_State* L, string ctorName,
	      out double r0, out double r1, out double s0, out double s1)
{
    // Table is at index 1
    string errMsgTmplt = format("Error in call to %s:new.", ctorName);
    errMsgTmplt ~= "The value for variable '%s' is not valid. It should be a number value.";
    r0 = getNumberFromTable(L, 1, "r0", false, 0.0, true, format(errMsgTmplt, "r0"));
    r1 = getNumberFromTable(L, 1, "r1", false, 1.0, true, format(errMsgTmplt, "r1"));
    s0 = getNumberFromTable(L, 1, "s0", false, 0.0, true, format(errMsgTmplt, "s0"));
    s1 = getNumberFromTable(L, 1, "s1", false, 1.0, true, format(errMsgTmplt, "s1"));
}

// --------------- Constructor specific to CoonsPatch -------------------

/**
 * This function implements our constructor from the Lua interface.
 *
 * At successful completion of this function, a new CoonsPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = CoonsPatch:new{north=nPath, east=ePath, south=sPath, west=wPath}
 * patch1 = CoonsPatch:new{north=nPath, east=ePath, south=sPath, west=wPath,
 *                         r0=0.0, r1=1.0, s0=0.0, s1=1.0}
 * patch2 = CoonsPatch:new{p00=a, p10=b, p11=c, p01=d}
 * patch3 = CoonsPatch:new{p00=a, p10=b, p11=c, p01=d, r0=0, r1=1, s0=0, s1=1}
 * --------------------------
 * Notes:
 * 1. See PJs diagram at top of geom.surface.d for ordering and labelling of
 *    paths and corners.
 * 2. Any missing r and s parameters default to defaults given in the
 *    CoonsPatch constructor.
 * 3. No mix-n-match of constructors allowed. It is one of: 4 paths OR 4 corner points. If a path is found first, that constructor wins.
 *
 */

extern(C) int newCoonsPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor CoonPatch:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Look for a path. If found, proceed with construction from paths.
    lua_getfield(L, 1, "north");
    if ( ! lua_isnil(L, -1) ) { 
	lua_pop(L, 1);
	Path[string] paths;
	getPaths(L, "CoonsPatch", paths);
	double r0, r1, s0, s1;
	getRandS(L, "CoonsPatch", r0, r1, s0, s1);
	auto cpatch = new CoonsPatch(paths["south"], paths["north"],
				     paths["west"], paths["east"],
				     r0, r1, s0, s1);
	return pushObj!(CoonsPatch, CoonsPatchMT)(L, cpatch);
    }
    else {
	lua_pop(L, 1);
    }
    // Instead, look for Vector3 objects.
    lua_getfield(L, 1, "p00");
    if ( ! lua_isnil(L, -1) ) {
	lua_pop(L, 1);
	Vector3 p00, p10, p11, p01;
	getVector3s(L, "CoonsPath", p00, p10, p11, p01);
	double r0, r1, s0, s1;
	getRandS(L, "CoonsPatch", r0, r1, s0, s1);
	auto cpatch = new CoonsPatch(p00, p10, p11, p01, r0, r1, s0, s1);
	return pushObj!(CoonsPatch, CoonsPatchMT)(L, cpatch);
    }
    lua_pop(L, 1);
    // If we make it here, there's been an error in construction.
    string errMsg = `There's a problem in call to CoonsPath.new.
Neither a list of named paths ('north', 'east', 'south', 'west')
nor a list of named corners ('p00', 'p10', 'p11', 'p01') were found.`;
    luaL_error(L, errMsg.toStringz);
    return 0;
}

void registerSurfaces(LuaState lua)
{
    auto L = lua.state;
    
    // Register the CoonsPatch object
    luaL_newmetatable(L, CoonsPatchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newCoonsPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(CoonsPatch, CoonsPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(CoonsPatch, CoonsPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(CoonsPatch, CoonsPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, CoonsPatchMT.toStringz);
}
