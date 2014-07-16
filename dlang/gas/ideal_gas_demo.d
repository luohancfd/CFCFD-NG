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
    writefln("R= %s, pressure= %s, temperature= %s", gm.R(gd), gd.p, gd.T[0]);
    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    writefln("rho= %s, e= %s, a= %s", gd.rho, gd.e[0], gd.a);
    gm.update_trans_coeffs(gd);
    writefln("mu= %s, k= %s", gd.mu, gd.k[0]);
    writeln("Done.");
}
