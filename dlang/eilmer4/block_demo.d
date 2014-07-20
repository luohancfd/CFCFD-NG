// block_demo.d
// Exercise the Structured-grid blocks as we build them.
// PJ 2014-07-20

import std.stdio;
import geom;
import gasmodelutil;
import flowstate;
import sblock;

void main()
{
    writeln("Block demo, structured blocks only for now...");
    auto gm = init_gas_model("sample-data/ideal-air-gas-model.json");
    auto flow = new FlowState(gm, 100.0e3, [300.0,], Vector3(1.0,0.0,0.0));
    writeln("flow=", flow);

    auto blk = new SBlock(1, "sample-data/sample-block.json");
    writeln("blk=", blk);

    writeln("Done.");
}
