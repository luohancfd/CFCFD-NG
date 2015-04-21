/**
 * luagpath_demo.d Demonstrate some of the behaviour of the Path primitives.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-02-22
 */

import std.stdio;
import luad.all;
import geom;
import gpath;
import luageom;
import luagpath;

void main()
{
    writeln("Begin demonstration of LuaD connection to Paths.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    registerPaths(lua);
    lua.doString(`
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
print("Bezier")
adb = Bezier:new{a, d, b}
e = adb(0.5)
print("e=", e, "expected approximately Vector3(1.60355, 2, 0.603553)")
    `);
    writeln("Done with demo.");
}
