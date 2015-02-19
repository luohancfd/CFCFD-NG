/**
 * gpath_demo.d Demonstrate some of the behaviour of the geometric primitives.
 *
 * Author: Peter J.
 * Version: 2014-06-16
 */

import std.stdio;
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
    auto c = ab.eval(0.5);
    writeln("ab.eval(0.5)= ", c);
    writeln("Done.");
}
