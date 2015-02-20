/**
 * path3_demo.d Demonstrate some of the behaviour of the geometric primitives.
 *
 * Author: Peter J.
 * Version: 2014-06-16
 */

import std.stdio;
import geom;
import gpath;
import surface;

void main()
{
    writeln("Begin demonstration of the ParametricSurface elements.");
    auto p00 = Vector3([0.0, 0.1, 3.0]);
    auto p10 = Vector3(1.0, 0.1, 3.0);
    auto p11 = Vector3(1.0, 1.1, 3.0);
    auto p01 = Vector3(0.0, 1.1, 3.0);

    writeln("CoonsPatch demo");
    auto my_patch = new CoonsPatch(p00, p10, p11, p01);
    writeln("my_patch= ", my_patch);
    auto c = my_patch(0.1, 0.1);
    writeln("my_patch(0.1, 0.1)= ", c);

    writeln("AOPatch demo");
    p10 = Vector3(1.0, 0.4, 3.0);
    auto my_AOpatch = new AOPatch(p00, p10, p11, p01);
    writeln("my_AOpatch= ", my_AOpatch);
    c = my_AOpatch(0.1, 0.1);
    writeln("my_AOpatch(0.1, 0.1)= ", c);

    writeln("Done.");
}
