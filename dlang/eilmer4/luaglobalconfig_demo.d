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
GlobalConfig.title("Some grand plan")
print("title=", GlobalConfig.title())
GlobalConfig.dimensions(3)
print("dimensions=", GlobalConfig.dimensions())
GlobalConfig.dt_init(0.001)
print("dt_init=", GlobalConfig.dt_init())
print("axisymmetric=", GlobalConfig.axisymmetric())
GlobalConfig.axisymmetric(0) -- doesn't accept true false
print("axisymmetric=", GlobalConfig.axisymmetric())
if GlobalConfig.axisymmetric then
   print("script sees axisymmetric as true")
else
   print("script sees axisymmetric as false")
end
    `);
    writeln("Done with demo.");
}
