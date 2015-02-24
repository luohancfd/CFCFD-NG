/**
 * luaflowstate_demo.d
 * Demonstrate the wrapped FlowState object.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-23
 */

import std.stdio;

import luad.all;
import luageom;
import luaflowstate;

void main()
{
    writeln("Begin demonstration of Lua connection to FlowState object.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    registerFlowState(lua);
    lua.doString(`
nsp, nmodes = setGasModel('sample-data/ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
fs = FlowState:new{p=1.0e5, T=300.0, u=1000.0, v=200.0}
fsTab = fs:toTable{}
for k,v in pairs(fsTab) do
    print(k,v)
    if ( k == 'gas' ) then
       for k1,v1 in pairs(v) do
           print(k1, v1)
           if ( k1 == 'massf' ) then
               print("massf[1]= ", v1[1])
           end
       end
    end
end
    `);
    writeln("Done with demo.");
}
