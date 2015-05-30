/**
 * luagpath_demo.d Demonstrate some of the behaviour of the Path primitives.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-02-22
 */

import std.stdio;
import util.lua;
import geom;
import gpath;
import luageom;
import luagpath;

void main()
{
    writeln("Begin demonstration of LuaD connection to Paths.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    luaL_dostring(L, `
-- Add a couple of points and alter their data.
a = Vector3:new{x=1.0, y=2.0}
b = Vector3:new(0.0, 5.0, 4.0)
ab = Line:new{a, b}
print("ab= ", ab)
print("Try evaluating a point midway on the line.")
pt = ab(0.5)
print("pt= ", pt)
print("Or with an eval.")
pt2 = ab:eval(0.5)
print("pt2= ", pt2)
ab:t1(0.8)
ab2 = ab:copy()
ab2:t0(0.2)
print("ab:t0()= ", ab:t0(), "ab:t1()= ", ab:t1())
print("ab2:t0()= ", ab2:t0(), "ab2:t2()= ", ab2:t1())
--
print("Arc")
a = Vector3:new(2.0, 2.0, 0.0)
b = Vector3:new(1.0, 2.0, 1.0)
c = Vector3:new(1.0, 2.0, 0.0)
abc = Arc:new{a, b, c}
d = abc(0.5)
print("d=", d, "expected approximately Vector3(1.7071068, 2.0, 0.7071068)")
--
print("Arc3")
a = Vector3:new(2.0, 2.0, 0.0)
b = Vector3:new(1.0, 2.0, 1.0)
m = Vector3:new(1.7071068, 2.0, 0.7071068)
amb = Arc3:new{a, m, b}
dd = amb(0.5)
print("dd=", dd, "expected approximately Vector3(1.7071068, 2.0, 0.7071068)")
--
print("Bezier")
adb = Bezier:new{a, d, b}
e = adb(0.5)
print("e=", e, "expected approximately Vector3(1.60355, 2, 0.603553)")
--
print("Polyline")
polyline = Polyline:new{abc, Line:new{b,c}}
print("polyline= ", polyline)
f = polyline(0.5)
print("polyline(0.5)= ", f, "expected approximately Vector3(1.28154, 2, 0.95955)")
--
print("LuaFnPath")
function myLuaFunction(t)
   -- Straight line from 0,0,0 to 1.0,2.0,3.0
   return {t, 2*t, 3*t}
end
myPath = LuaFnPath:new{"myLuaFunction"}
print("myLine= ", myPath)
g = myPath(0.5)
print("myPath(0.5)= ", g)
    `);
    writeln("Done with demo.");
}
