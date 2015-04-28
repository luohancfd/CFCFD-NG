/**
 * luasolidprops.d
 * Lua interface to access writeInitialSolidFile
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-04-29
 */

module luasolidprops;

import std.conv;
import luad.all;
import luad.stack;
import luad.c.lua;
import luad.c.lauxlib;
import util.lua_service;
import luasgrid;

import solidprops;

extern(C) int writeInitialSolidFileFromLua(lua_State* L)
{
    auto fname = to!string(luaL_checkstring(L, 1));
    auto grid = checkStructuredGrid(L, 2);
    double T_init = luaL_checknumber(L, 3);
    double rho = getNumberFromTable(L, 4, "rho", true);
    double k = getNumberFromTable(L, 4, "k", true);
    double Cp = getNumberFromTable(L, 4, "Cp", true);
    auto sp = new SolidProps(rho, k, Cp);
    double t0 = luaL_checknumber(L, 5);
    writeInitialSolidFile(fname, grid, T_init, sp, t0);
    return 0;
}

void registerSolidProps(LuaState lua)
{
    auto L = lua.state;
    lua_pushcfunction(L, &writeInitialSolidFileFromLua);
    lua_setglobal(L, "writeInitialSolidFile");
}
