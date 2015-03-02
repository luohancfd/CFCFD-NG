/**
 * luaglobalconfig_demo.d
 * Demonstrate the wrapped GlobalConfig object.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-03-02
 */

import std.stdio;

import luad.all;
import globalconfig;
import luaglobalconfig;

void main()
{
    writeln("Begin demonstration of Lua connection to GlobalConfig object.");
    auto lua = new LuaState;
    lua.openLibs();
    registerGlobalConfig(lua);
    lua.doString(`
-- Gas Model 
nsp, nmodes = setGasModel('sample-data/ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
--
-- Access GlobalConfig via a specially named 'config' table.
print("Testing a different form of config interaction.")
config{title="Testing plan 3", dimensions=3}
config["max_time"] = 1.0e-5
config.dt_init = 1.0e-9
print("title= ", config["title"])
print("dimensions= ", config.dimensions)
print("max_time= ", config.max_time)
print("dt_init= ", config["dt_init"])
    `);
    writeln("Done with demo.");
}
