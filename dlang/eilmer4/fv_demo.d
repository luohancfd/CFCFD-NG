// fv_demo.d
// Exerside the finite-volume classes.
// PJ, 2014-07-17

import std.stdio;
import geom;
import gasmodelutil;
import flowstate;
import fvcore;

void main()
{
    writeln("test fv modules");
    auto gm = init_gas_model("ideal-air-gas-model.json");
    auto flow = new FlowState(gm, 100.0e3, [300.0,], Vector3(1.0,0.0,0.0));
    writeln("flow=", flow);

    init_fvcore();
    writeln("faceName[faceIndex[\"south\"]]=", faceName[faceIndex["south"]]);
    writeln("updateSchemeName[updateSchemeIndex[\"classic_rk3\"]]=",
	    updateSchemeName[updateSchemeIndex["classic_rk3"]]);
    writeln("done.");
}
