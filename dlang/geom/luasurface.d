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
import std.conv;
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
immutable string AOPatchMT = "AOPatch";

ParametricSurface checkSurface(lua_State* L, int index) {
    // We have to do a brute force test for each object
    // type in turn.
    if ( isObjType(L, index, CoonsPatchMT) )
	return checkObj!(CoonsPatch, CoonsPatchMT)(L, index);
    if ( isObjType(L, index, AOPatchMT ) )
	return checkObj!(AOPatch, AOPatchMT)(L, index);
    // if no match found then
    return null;
}

extern(C) int isSurface(lua_State* L)
{
    if ( checkSurface(L, 1) )
	lua_pushboolean(L, 1);
    else
	lua_pushboolean(L, 0);
    return 1;
}

extern(C) int opCallSurface(T, string MTname)(lua_State* L)
{
    auto surface = checkObj!(T, MTname)(L, 1);
    auto r = luaL_checknumber(L, 2);
    auto s = luaL_checknumber(L, 3);
    auto pt = surface(r, s);
    return pushVector3(L, pt);
}

string getPathFromTable(string pName)
{
    return `lua_getfield(L, index, "`~pName~`");
    paths["`~pName~`"] = checkPath(L, -1);
    if ( paths["`~pName~`"] is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "`~pName~`")));
    lua_pop(L, 1);`;
}

void getPaths(lua_State *L, string ctorName, out Path[string] paths)
{
    // Assume that table is at index 1.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the '%s' path was not of Path type.`; 
    int index = 1;
    mixin(getPathFromTable("north"));
    mixin(getPathFromTable("east"));
    mixin(getPathFromTable("south"));
    mixin(getPathFromTable("west"));
}

string getVecFromTable(string vName)
{
    return `lua_getfield(L, index, "`~vName~`");
auto `~vName~`Ptr = checkVector3(L, -1);
if ( `~vName~`Ptr is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "`~vName~`")));
    else `~vName~` = *`~vName~`Ptr;
    lua_pop(L, 1);`;
}

void getVector3s(lua_State *L, string ctorName,
		 out Vector3 p00, out Vector3 p10, out Vector3 p11, out Vector3 p01)
{
    // Assume that table is at index 1.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the '%s' corner was not of Vector3 type.`; 
    int index = 1;
    mixin(getVecFromTable("p00"));
    mixin(getVecFromTable("p10"));
    mixin(getVecFromTable("p11"));
    mixin(getVecFromTable("p01"));
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

// --------------- Specific constructors for each Surface type -------------------

/**
 * This is the constructor for a CoonsPatch to be used from the Lua interface.
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
	getVector3s(L, "CoonsPatch", p00, p10, p11, p01);
	writeln("p00= ", p00, " p10= ", p10, " p11= ", p11, " p01= ", p01);
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

/**
 * This is constructor for an AOPatch object to be used from the Lua interface.
 *
 * At successful completion of this function, a new AOPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = AOPatch:new{north=nPath, east=ePath, south=sPath, west=wPath}
 * patch1 = AOPatch:new{north=nPath, east=ePath, south=sPath, west=wPath,
 *                     r0=0.0, r1=1.0, s0=0.0, s1=1.0, nx=10, ny=10}
 * patch2 = AOPatch:new{p00=a, p10=b, p11=c, p01=d}
 * patch3 = AOPatch:new{p00=a, p10=b, p11=c, p01=d, r0=0, r1=1, s0=0, s1=1, nx=10, ny=10}
 * --------------------------
 * Notes:
 * 1. See PJs diagram at top of geom.surface.d for ordering and labelling of
 *    paths and corners.
 * 2. Any missing r and s parameters default to defaults given in the
 *    CoonsPatch constructor.
 * 3. No mix-n-match of constructors allowed. It is one of: 4 paths OR 4 corner points. If a path is found first, that constructor wins.
 *
 */
extern(C) int newAOPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor AOPatch:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    string errMsgTmplt = "Error in call to AOPatch:new.\n";
    errMsgTmplt ~= "A valid value for '%s' is not found in arguments.\n";
    errMsgTmplt ~= "The value, if present, should be a number.";
    int nx = to!int(getNumberFromTable(L, 1, "nx", false, 10.0, true, format(errMsgTmplt, "nx")));
    int ny = to!int(getNumberFromTable(L, 1, "ny", false, 10.0, true, format(errMsgTmplt, "ny")));
    // Look for a path. If found, proceed with construction from paths.
    lua_getfield(L, 1, "north");
    if ( ! lua_isnil(L, -1) ) { 
	lua_pop(L, 1);
	Path[string] paths;
	getPaths(L, "AOPatch", paths);
	double r0, r1, s0, s1;
	getRandS(L, "AOPatch", r0, r1, s0, s1);
	

	auto aopatch = new AOPatch(paths["south"], paths["north"],
				   paths["west"], paths["east"],
				   nx, ny, r0, r1, s0, s1);
	return pushObj!(AOPatch, AOPatchMT)(L, aopatch);
    }
    else {
	lua_pop(L, 1);
    }
    // Instead, look for Vector3 objects.
    lua_getfield(L, 1, "p00");
    if ( ! lua_isnil(L, -1) ) {
	lua_pop(L, 1);
	Vector3 p00, p10, p11, p01;
	getVector3s(L, "AOPatch", p00, p10, p11, p01);
	writeln("p00= ", p00, " p10= ", p10, " p11= ", p11, " p01= ", p01);
	double r0, r1, s0, s1;
	getRandS(L, "AOPatch", r0, r1, s0, s1);
	auto aopatch = new AOPatch(p00, p10, p11, p01, nx, ny, r0, r1, s0, s1);
	return pushObj!(AOPatch, AOPatchMT)(L, aopatch);
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

    // Register the AOPatch object
    luaL_newmetatable(L, AOPatchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newAOPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(AOPatch, AOPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(AOPatch, AOPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(AOPatch, AOPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, AOPatchMT.toStringz);

    lua_pushcfunction(L, &isSurface);
    lua_setglobal(L, "isSurface");
}
