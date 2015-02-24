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
surf = CoonsPatch:new{north=Line:new(b, d), east=Line:new(c, d),
                      south=Line:new(a, c), west=Line:new(a, b)}
print("CoonsPatch representation: ", surf)
ctr = surf(0.5, 0.5)
print("ctr= ", ctr)
-- Try construction using corners and a subset surface 
surf2 = CoonsPatch:new{p00=a, p01=b, p11=c, p10=d,
                       r0=0, r1=0.5, s0=0.5, s1=1}
p = surf2:eval(0.5, 0.5)
print("p= ", p)
-- Copy PJ's AO patch demo
p00 = Vector3:new{0.0, 0.1, 3.0}
p10 = Vector3:new{1.0, 0.4, 3.0}
p11 = Vector3:new{1.0, 1.1, 3.0}
p01 = Vector3:new{0.0, 1.1, 3.0}
my_aopatch = AOPatch:new{p00=p00, p10=p10, p11=p11, p01=p01}
p = my_aopatch(0.1, 0.1);
print("my_aopatch(0.1, 0.1)= ", p)
print("isSurface(my_aopatch)= ", isSurface(my_aopatch))
print("isSurface(surf2)= ", isSurface(surf2));
print("isSurface(a)= ", isSurface(a));

    `);
}

    
