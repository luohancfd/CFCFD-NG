/** e4_core.d
 * Eilmer4 compressible-flow simulation code, core coordination functions.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

import std.stdio;

import geom;
import gas;
import globalconfig;
import flowstate;
import sblock;

//-----------------------------------------------------------------------------------------

double init_simulation()
{
    writeln("Initialize simulation.");
    writeln("TODO fill in the REAL details.");
    string fileName = GlobalConfig.base_file_name ~ ".config";
    writeln("Read config from file: ", fileName);
    GlobalConfig.gmodel = init_gas_model("sample-data/ideal-air-gas-model.lua");
    GlobalConfig.flow_state ~= new FlowState(GlobalConfig.gmodel, 100.0e3, [300.0,],
					     Vector3(1.0,0.0,0.0));
    writeln("flow=", GlobalConfig.flow_state[0]);

    double sim_time;
    GlobalConfig.blk_data ~= new SBlock(1, "sample-data/sample-block.json");
    foreach (ref myblk; GlobalConfig.blk_data) {
	myblk.assemble_arrays();
	myblk.bind_faces_and_vertices_to_cells();
	writeln("myblk=", myblk);
	myblk.read_grid("sample-data/cone20.grid.b0000.t0000.gz", 0);
	myblk.write_grid("test-grid.txt.gz", 0.0);
	sim_time = myblk.read_solution("sample-data/cone20.flow.b0000.t0000.gz");
    }
    writeln("Done init_simulation.");
    return sim_time;
} // end init_simulation()

double integrate_in_time(double sim_time)
{
    writeln("Integrate in time.");
    writeln("TODO fill in the REAL details.");
    writeln("Done integrate_in_time.");
    return sim_time;
} // end integrate_in_time()

void finalize_simulation(double sim_time)
{
    writeln("Finalize the simulation.");
    writeln("TODO fill in the REAL details.");
    foreach (ref myblk; GlobalConfig.blk_data) {
	myblk.write_solution("test-flow.txt.gz", 1.0);
    }
    writeln("Done finalize_simulation.");
} // end finalize_simulation()
