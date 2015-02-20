/**
 * gpath_demo.d Demonstrate some of the behaviour of the geometric primitives.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-21
 */

import std.stdio;
import luad.all;
import geom;
import gpath;
import luageom;

void main()
{
    writeln("Begin demonstration of LuaD connection to geometric elements.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    registerPath(lua);
    lua.doString(`
-- Add a couple of points and alter their data.
a = Vector{x=1.0, y=2.0}
print("a=", a, " remember that it is an index")
print("a.x=", getX(a), "a.y=", getY(a), "a.z=", getZ(a))
b = {}
Vector3Value(b, a)
print("b=", b)
print("uglyprint b=[")
for k,v in pairs(b) do
   print(k, "=", v, ",")
end
print("]")
c = Vector3({x=2.0, z=-4.2})
print("c.x= ", c.x, "c.y= ", c.y, "c.z=", c.z)
c.x = 32.0
print("After modification of c.x...")
print("c.x= ", c.x, "c.y= ", c.y, "c.z=", c.z)      

-- Try a path
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
p = Vector3{x=9.0}
q = Vector3{x=1.0, y=-1.0, z=8}
pq = Line2(p, q)
r = evalLine2(pq, 0.2)
print("r.x=", r.x, "r.y=", r.y, "r.z=", r.z)
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
