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
    writeln("a = ", a, ", b = ", b);
    writeln("Done.");
    Vector3 c = a + b;
    writeln("c = a + b = ", c);
    Vector3 d = a.dup;
    a.y = 99.0;
    writeln("a = ", a, ", d = ", d); 
    Vector3 e = a * 2.0;
    Vector3 f = 3 * d; // int promoted to double OK
    writeln("e = ", e, " f = ", f);
    Vector3 g = d / 3.0;
    writeln("g = ", g);
    g += f;
    writeln("g += f -> ", g);
    g /= 2.0;
    writeln("g /= 2.0 -> ", g);
}
