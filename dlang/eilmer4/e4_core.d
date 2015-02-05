/** e4_core.d
 * Eilmer4 compressible-flow simulation code, core coordination functions.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

import std.stdio;
import std.json;
import std.file;
import std.conv;

import json_helper;
import geom;
import gas;
import globalconfig;
import flowstate;
import sblock;

//-----------------------------------------------------------------------------------------

// Flow condition array for use in boundary conditions.
static FlowState[] flow_state;

// Storage for the actual blocks of flow data.
static SBlock[] blk_data;

//-----------------------------------------------------------------------------------------

void read_control_file()
{
    if (GlobalConfig.verbosity_level > 1) writeln("read_control_file()");
    string fileName = GlobalConfig.base_file_name ~ ".control";
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
	writeln("Failed to read control file: ", fileName);
	exit(1);
    }
    JSONValue jsonData;
    try {
	jsonData = parseJSON!string(content);
    } catch (Exception e) {
	writeln("Failed to parse JSON from control file: ", fileName);
	exit(1);
    }
    GlobalConfig.Xorder = getJSONint(jsonData, "x_order", 2);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  Xorder: ", GlobalConfig.Xorder);
    }
} // end read_control_file()


double init_simulation()
{
    if (GlobalConfig.verbosity_level > 0) writeln("Begin init_simulation...");

    if (GlobalConfig.verbosity_level > 1) writeln("Read config file.");
    string fileName = GlobalConfig.base_file_name ~ ".config";
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
	writeln("Failed to read config file: ", fileName);
	exit(1);
    }
    JSONValue jsonData;
    try {
	jsonData = parseJSON!string(content);
    } catch (Exception e) {
	writeln("Failed to parse JSON from config file: ", fileName);
	exit(1);
    }
    GlobalConfig.title = jsonData["title"].str;
    string gasModelFile = jsonData["gas_model_file"].str;
    GlobalConfig.nBlocks = getJSONint(jsonData, "nblock", 0);
    GlobalConfig.dimensions = getJSONint(jsonData, "dimensions", 2);
    GlobalConfig.axisymmetric = getJSONbool(jsonData, "axisymmetric_flag", false);
    GlobalConfig.viscous = getJSONbool(jsonData, "viscous_flag", false);
    GlobalConfig.viscous = getJSONbool(jsonData, "viscous_flag", false);
    GlobalConfig.viscous_delay = getJSONdouble(jsonData, "viscous_delay", 0.0);
    GlobalConfig.viscous_factor_increment = getJSONdouble(jsonData, "viscous_factor_increment", 0.01);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  title: ", GlobalConfig.title);
	writeln("  gasModelFile: ", gasModelFile);
	writeln("  nBlocks: ", GlobalConfig.nBlocks);
	writeln("  dimensions: ", GlobalConfig.dimensions);
	writeln("  axisymmetric: ", GlobalConfig.axisymmetric);
	writeln("  viscous: ", GlobalConfig.viscous);
	writeln("  viscous_delay: ", GlobalConfig.viscous_delay);
	writeln("  viscous_factor_increment: ", GlobalConfig.viscous_factor_increment);
    }
    GlobalConfig.gmodel = init_gas_model(gasModelFile);

    writeln("TODO -------- fill in the REAL details ------------");
    flow_state ~= new FlowState(GlobalConfig.gmodel, 100.0e3, [300.0,],
				Vector3(1.0,0.0,0.0));
    writeln("flow=", flow_state[0]);
    double sim_time;
    blk_data ~= new SBlock(1, "sample-data/sample-block.json");
    foreach (ref myblk; blk_data) {
	myblk.assemble_arrays();
	myblk.bind_faces_and_vertices_to_cells();
	writeln("myblk=", myblk);
	myblk.read_grid("sample-data/cone20.grid.b0000.t0000.gz", 0);
	myblk.write_grid("test-grid.txt.gz", 0.0);
	sim_time = myblk.read_solution("sample-data/cone20.flow.b0000.t0000.gz");
    }

    if (GlobalConfig.verbosity_level > 0) writeln("Done init_simulation.");
    return sim_time;
} // end init_simulation()

double integrate_in_time(double sim_time)
{
    writeln("Integrate in time.");
    read_control_file(); // every step
    writeln("TODO fill in the REAL details.");
    writeln("Done integrate_in_time.");
    return sim_time;
} // end integrate_in_time()

void finalize_simulation(double sim_time)
{
    writeln("Finalize the simulation.");
    writeln("TODO fill in the REAL details.");
    foreach (ref myblk; blk_data) {
	myblk.write_solution("test-flow.txt.gz", 1.0);
    }
    writeln("Done finalize_simulation.");
} // end finalize_simulation()
