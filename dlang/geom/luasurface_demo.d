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

void main()
{
    writeln("Begin demonstration of Lua connection to Surface objects.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    registerPaths(lua);
    registerSurfaces(lua);
    lua.doString(`
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.0, 1.0}
c = Vector3:new{1.0, 0.0}
d = Vector3:new{1.0, 1.0}
bd = Line:new(b, d)
cd = Line:new(c, d)
ac = Line:new(a, c)
ab = Line:new(a, b)
surf = CoonsPatch:new{north=bd, east=cd, south=ac, west=ab}
--print("CoonsPatch representation: ", surf)
ctr = surf(0.5, 0.5)
print("ctr= ", ctr)
    `);
}

    
