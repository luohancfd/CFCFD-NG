/**
 * luasurface_demo.d
 * Demonstrate the wrapped Surface objects.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-24
 */

import std.stdio;

import luad.all;
import luageom;
import luagpath;
import luasurface;
import luaunifunction;
import luasgrid;

void main()
{
    writeln("Begin demonstration of Lua connection to StructuredGrid objects.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    registerPaths(lua);
    registerSurfaces(lua);
    registerUnivariateFunctions(lua);
    registerStructuredGrid(lua);
    lua.doString(`
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.0, 1.0}
c = Vector3:new{1.0, 0.0}
d = Vector3:new{1.0, 1.0}
surf = CoonsPatch:new{north=Line:new{b, d}, east=Line:new{c, d},
                      south=Line:new{a, c}, west=Line:new{a, b}}
print("CoonsPatch representation: ", surf)
myrf = RobertsFunction:new{end0=true, end1=false, beta=1.01}
grid = StructuredGrid2D:new{surf=surf, niv=10, njv=10} 
    `);
}

    
