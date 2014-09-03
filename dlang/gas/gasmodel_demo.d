/**
 * gasmodel_demo.d
 *
 * Author: Peter J.
 * Version: 2014-06-22
 */

import std.stdio;
import gasmodel;
import gasmodelutil;

void main() {
    writeln("Begin demonstration of using the gasmodel and Gas_data classes...");
    auto gm = init_gas_model("sample-data/ideal-air-gas-model.json");
    foreach(i; 0 .. gm.n_species) {
	writeln("species[", i, "] name=", gm.species_name(i));
    }
    auto gd = GasState(gm, 100.0e3, 300.0);
    writefln("R= %s, pressure= %s, temperature= %s", gm.R(gd), gd.p, gd.T[0]);
    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    writefln("rho= %s, e= %s, a= %s", gd.rho, gd.e[0], gd.a); 
    writeln("gd= ", gd);
    auto gd2 = GasState(gm, 200.0e3, 400.0);
    writeln("gd2=", gd2);
    auto gd3 = GasState(gm, 100.0e3, 300.0);
    gd3.copy_average_values_from([gd2, gd], gm);
    writeln("after average gd3=", gd3);
    gd2.copy_values_from(gd);
    writeln("after copy gd2=", gd2);
    writeln("Done.");
}
