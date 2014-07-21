// block_demo.d
// Exercise the Structured-grid blocks as we build them.
// PJ 2014-07-20

import std.stdio;
import geom;
import gasmodelutil;
import globalconfig;
import flowstate;
import sblock;

void main()
{
    writeln("Block demo, structured blocks only for now...");
    GlobalConfig.gmodel = init_gas_model("sample-data/ideal-air-gas-model.json");
    GlobalConfig.verbosity_level = 2;
    auto flow = new FlowState(GlobalConfig.gmodel, 100.0e3, [300.0,], Vector3(1.0,0.0,0.0));
    writeln("flow=", flow);

    auto blk = new SBlock(1, "sample-data/sample-block.json");
    blk.assemble_arrays();
    blk.bind_faces_and_vertices_to_cells();
    writeln("blk=", blk);

    writeln("Done.");
}
