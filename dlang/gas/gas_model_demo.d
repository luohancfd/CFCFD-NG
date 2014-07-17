/**
 * gas_model_demo.d
 *
 * Author: Peter J.
 * Version: 2014-06-22
 */

import std.stdio;
import gas_model;
import gas_model_util;

void main() {
    writeln("Begin demonstration of using the Gas_model and Gas_data classes...");
    auto gm = init_gas_model("ideal-air-gas-model.json");
    auto gd = new Gas_data(gm, 100.0e3, 300.0);
    writefln("R= %s, pressure= %s, temperature= %s", gm.R(gd), gd.p, gd.T[0]);
    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    writefln("rho= %s, e= %s, a= %s", gd.rho, gd.e[0], gd.a); 
    writeln("gd= ", gd);
    writeln("Done.");
}
