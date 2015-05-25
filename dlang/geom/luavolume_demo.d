/**
 * luavolume_demo.d
 * Demonstrate the wrapped ParametricVolume objects.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2015-04-07
 */

import std.stdio;

import util.lua;
import luageom;
import luagpath;
import luasurface;
import luavolume;

void main()
{
    writeln("Begin demonstration of Lua connection to Volume objects.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    registerSurfaces(L);
    registerVolumes(L);
    luaL_dostring(L, `
p000 = Vector3:new{0.0, 0.1, 0.0}
p100 = Vector3:new{1.0, 0.1, 0.0}
p110 = Vector3:new{1.0, 1.1, 0.0}
p010 = Vector3:new{0.0, 1.1, 0.0}
p001 = Vector3:new{0.0, 0.1, 3.0}
p101 = Vector3:new{1.0, 0.1, 3.0}
p111 = Vector3:new{1.0, 1.1, 3.0}
p011 = Vector3:new{0.0, 1.1, 3.0}
my_volume = TFIVolume:new{vertices={p000,p100,p110,p010,p001,p101,p111,p011}}
print("my_volume=", my_volume)
p = my_volume(0.1, 0.1, 0.5);
print("my_volume(0.1, 0.1, 0.5)= ", p)
print("isVolume(my_volume)= ", isVolume(my_volume))
    `);
}

    
