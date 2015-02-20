/**
 * gpath_demo.d Demonstrate some of the behaviour of the geometric primitives.
 *
 * Author: Peter J.
 * Version: 2014-06-16
 */

import std.stdio;
import luad.all;
import geom;
import gpath;

void main()
{
    writeln("Begin demonstration of the geometric Path primitives.");
    auto a = Vector3([1.0, 2.2, 3.0]);
    auto b = Vector3(1.0);
    writeln("a = ", a, ", b = ", b);
    auto ab = new Line(a, b);
    writeln("ab= ", ab);
    auto c = ab(0.5);
    writeln("ab(0.5)= ", c);

    writeln("Try LuaD connection.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    registerPath(lua);
    lua.doString(`
-- Add a couple of points and make a Line.
a = Vector{x=1.0, y=2.0}
b = Vector{x=3.0, y=5.0}
ab = Line(a, b)
c = {}
evalLine(c, ab, 0.5)
print("c=", c)
print("uglyprint c=[")
for k,v in pairs(c) do
   print(k, "=", v, ",")
end
print("]")
ef = Line(VectorA{0.0, 10.0}, VectorA{10.0, 0.0})      
    `);
    writeln("points.length= ", points.length);
    foreach (i; 0 .. points.length) {
	writeln("points[", i, "]= ", points[i]);
    }
    writeln("paths.length= ", paths.length);
    foreach (i; 0 .. paths.length) {
	writeln("paths[", i, "]= ", paths[i]);
    }
    writeln("Done.");
}
