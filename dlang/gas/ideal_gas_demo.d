/**
 * ideal_gas_demo.d
 *
 * Author: Peter J.
 * Version: 2014-06-22
 */

import std.stdio;
import gas_model;
import ideal_gas;

void main() {
    writeln("Begin demonstration of using the Ideal_gas and Gas_data classes...");
    auto gm = new Ideal_gas();
    auto gd = new Gas_data(gm);
    writeln("pressure=", gd.p);
    writeln("Done.");
}
