/**
 * luageom2_demo.d
 * Shows the wrapped Vector3 in action.
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-02-21
 */

import std.stdio;
import luad.all;
import luageom2;

void main()
{
    writeln("Begin demo of wrapped D Vector3 struct for use in Lua.");
    auto lua = new LuaState;
    lua.openLibs();
    registerVector3(lua);
    lua.doString(`
-- Add some points and manipulate them.
a = Vector3:new()
print("a= ", a)
print("Change a's x value.")
a:x(4.4)
print("a:x= ", a:x())
b = Vector3:new(7.0, 3.0, -2.5)
print("b= ", b)
c = a + b
print("c= ", c)
d = add(a, b)
print("d= ", d)
d:normalize()
print("After normalizing, d=", d)
    `);
    writeln("End demo.");
}

