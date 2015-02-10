/** e4_core.d
 * Eilmer4 compressible-flow simulation code, core coordination functions.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;

import geom;
import gas;
import fvcore;
import globalconfig;
import readconfig;
import globaldata;
import flowstate;
import sblock;
import bc;

static int current_tindx = 0;

double init_simulation(int tindx)
{
    if (GlobalConfig.verbosity_level > 0) writeln("Begin init_simulation()...");
    read_config_file();
    current_tindx = tindx;
    double sim_time;
    foreach (ref myblk; myBlocks) {
	myblk.assemble_arrays();
	myblk.bind_faces_and_vertices_to_cells();
	writeln("myblk=", myblk);
	auto writer = appender!string();
	formattedWrite(writer, "grid/t%04d/%s.grid.b%04d.t%04d.gz",
		       tindx, GlobalConfig.base_file_name, myblk.id, tindx);
	auto fileName = writer.data();
	myblk.read_grid(fileName, 0);
	writer = appender!string();
	formattedWrite(writer, "flow/t%04d/%s.flow.b%04d.t%04d.gz",
		       tindx, GlobalConfig.base_file_name, myblk.id, tindx);
	fileName = writer.data();
	sim_time = myblk.read_solution(fileName);
    }
    foreach (ref myblk; myBlocks) {
	myblk.compute_primary_cell_geometric_data(0);
	myblk.compute_distance_to_nearest_wall_for_all_cells(0);
	myblk.compute_secondary_cell_geometric_data(0);
	myblk.identify_reaction_zones(0);
	myblk.identify_turbulent_zones(0);
	foreach (ref cell; myblk.active_cells) {
	    cell.encode_conserved(0, 0, myblk.omegaz);
	    // Even though the following call appears redundant at this point,
	    // fills in some gas properties such as Prandtl number that is
	    // needed for both the cfd_check and the BLomax turbulence model.
	    cell.decode_conserved(0, 0, myblk.omegaz);
	}
    }
    exchange_shared_boundary_data();
    if (GlobalConfig.verbosity_level > 0) writeln("Done init_simulation().");
    return sim_time;
} // end init_simulation()

void exchange_shared_boundary_data()
{
    foreach (ref myblk; myBlocks) {
	foreach (face; 0 .. (GlobalConfig.dimensions == 3 ? 6 : 4)) {
	    if (myblk.bc[face].type_code  == BCCode.full_face_exchange) {
		myblk.bc[face].do_copy_into_boundary();
	    }
	} // end foreach face
    } // end foreach myblk
} // end exchange_shared_boundary_data()

double integrate_in_time(double sim_time)
{
    writeln("Integrate in time.");
    read_control_file(); // every step
    writeln("TODO fill in the REAL details.");
    writeln("Done integrate_in_time().");
    return sim_time;
} // end integrate_in_time()

void finalize_simulation(double sim_time)
{
    writeln("Finalize the simulation.");
    writeln("TODO fill in the REAL details.");
    foreach (ref myblk; myBlocks) {
	myblk.write_solution("test-flow.txt.gz", 1.0);
    }
    writeln("Done finalize_simulation.");
} // end finalize_simulation()
