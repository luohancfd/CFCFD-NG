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
-------------------------- Gas Model -----------------------
nsp, nmodes = setGasModel('sample-data/ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
--
------------------ Plan 1: GlobalConfig via metatable functions -------------------
print("Test access to GlobalConfig via functions.")
GlobalConfig.title("Some grand plan")
print("title=", GlobalConfig.title())
GlobalConfig.dimensions(3)
print("dimensions=", GlobalConfig.dimensions())
GlobalConfig.dt_init(0.001)
print("dt_init=", GlobalConfig.dt_init())
--
print("default value: axisymmetric=", GlobalConfig.axisymmetric())
GlobalConfig.axisymmetric(true)
print("axisymmetric=", GlobalConfig.axisymmetric())
if GlobalConfig.axisymmetric() then
   print("script sees axisymmetric as true")
else
   print("script sees axisymmetric as false")
end
GlobalConfig.axisymmetric(false)
print("axisymmetric=", GlobalConfig.axisymmetric())
if GlobalConfig.axisymmetric() then
   print("script sees axisymmetric as true")
else
   print("script sees axisymmetric as false")
end
--
------------------ Plan 2: GlobalConfig input via a table -------------------------
print("Test setting of GlobalConfig via a table.")
configSet{title="A better plan", dimensions=2, max_time=0.00123}
configSet{axisymmetric=true, dt_init=0.002}
print("dimensions=", configGet("dimensions"))
print("title=", configGet("title"))
print("axisymmetric=", configGet("axisymmetric"))
print("max_time=", configGet("max_time"))
print("dt_init=", configGet("dt_init"))
------------------ Plan 3: GlobalConfig via a specially named 'config' table. ----
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
