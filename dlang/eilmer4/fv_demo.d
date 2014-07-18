// fv_demo.d
// Exerside the finite-volume classes.
// PJ, 2014-07-17

import std.stdio;
import geom;
import gasmodelutil;
import flowstate;
import conservedquantities;
import fvcore;
import fvinterface;
import fvvertex;

void main()
{
    writeln("test fv modules");
    auto gm = init_gas_model("ideal-air-gas-model.json");
    auto flow = new FlowState(gm, 100.0e3, [300.0,], Vector3(1.0,0.0,0.0));
    writeln("flow=", flow);
    auto flow2 = new FlowState(gm, 120.0e3, [300.0,], Vector3(0.0,1.0,0.0));
    auto flow3 = new FlowState(gm, 120.0e3, [300.0,], Vector3(0.0,1.0,0.0));
    flow3.copy_average_values_from([flow, flow2], gm);
    writeln("flow3=", flow3);

    auto Q = new ConservedQuantities(gm);
    Q.mass = 99.0;
    Q.energies[0] = 301.0;
    writeln("Q=", Q);
    Q.clear_values();
    writeln("cleared Q=", Q);

    init_fvcore();
    writeln("faceName[faceIndex[\"south\"]]=", faceName[faceIndex["south"]]);
    writeln("updateSchemeName[updateSchemeIndex[\"classic_rk3\"]]=",
	    updateSchemeName[updateSchemeIndex["classic_rk3"]]);

    auto iface = new FVInterface(gm);
    writeln("iface=", iface);

    auto vtx = new FVVertex(gm);
    writeln("vtx=", vtx);

    writeln("done.");
}
