/**
 * geom_demo.d Demonstrate some of the behaviour of the geometric primitives.
 *
 * Author: Peter J.
 * Version: 2014-06-16
 */

import std.stdio;
import geom;

void main()
{
    writeln("Begin demonstration of the geometric primitives.");
    Vector3 a = new Vector3([1.0, 2.2, 3.0]);
    Vector3 b = new Vector3(1.0);
    writeln("a = ", a.toString(), ", b = ", b.toString());
    writeln("Done.");
    Vector3 c = a + b;
    writeln("c = a + b = ", c.toString());
    Vector3 d = a.dup;
    a.y = 99.0;
    writeln("a = ", a.toString(), ", d = ", d.toString());   
}
