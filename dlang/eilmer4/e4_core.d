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
import std.algorithm;

import geom;
import gas;
import fvcore;
import globalconfig;
import readconfig;
import globaldata;
import flowstate;
import sblock;
import bc;

// State data for simulation.
// Needs to be seen by all of the coordination functions.
static int current_tindx = 0;
static double sim_time = 0;  // present simulation time, tracked by code
static int step = 0;
static double dt_global;     // simulation time step determined by code
static double dt_allow;      // allowable global time step determined by code
static double t_plot;        // time to write next soln
static bool output_just_written = true;
static double t_history;     // time to write next sample
static bool history_just_written = true;
static int wall_clock_at_start = 0; // seconds

string make_path_name(string file_type)(int blk_id, int tindx)
// Build a pathname for "grid" and "flow" files which are stored in
// subdirectories and labelled with the block id and time index counters.
{
    auto writer = appender!string();
    formattedWrite(writer, "%s/t%04d/%s.%s.b%04d.t%04d.gz",
		   file_type, tindx, GlobalConfig.base_file_name,
		   file_type, blk_id, tindx);
    return writer.data();
}

double init_simulation(int tindx)
{
    if (GlobalConfig.verbosity_level > 0) writeln("Begin init_simulation()...");
    wall_clock_at_start = 0; // TODO read system time in seconds
    read_config_file();  // most of the configuration is in here
    read_control_file(); // some of the configuration is in here
    current_tindx = tindx;
    double sim_time;
    foreach (ref myblk; myBlocks) {
	myblk.assemble_arrays();
	myblk.bind_faces_and_vertices_to_cells();
	writeln("myblk=", myblk);
	myblk.read_grid(make_path_name!"grid"(myblk.id, tindx), 0);
	sim_time = myblk.read_solution(make_path_name!"flow"(myblk.id, tindx));
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
	    if (myblk.bc[face].type_code == BCCode.full_face_exchange) {
		myblk.bc[face].do_copy_into_boundary();
	    }
	} // end foreach face
    } // end foreach myblk
} // end exchange_shared_boundary_data()

double integrate_in_time(double target_time, int maxWallClock)
{
    if (GlobalConfig.verbosity_level > 0) writeln("Integrate in time.");
    // The next time for output...
    t_plot = sim_time + GlobalConfig.dt_plot;
    t_history = sim_time + GlobalConfig.dt_history;
    // Flags to indicate that the saved output is fresh.
    // On startup or restart, it is assumed to be so.
    output_just_written = true;
    history_just_written = true;
    // Overall iteration count.
    step = 0;
    // Normally, we can terminate upon either reaching 
    // a maximum time or upon reaching a maximum iteration count.
    bool finished_time_stepping = 
	(sim_time >= min(target_time, GlobalConfig.max_time) || 
	 step >= GlobalConfig.max_step);
    //----------------------------------------------------------------
    //                 Top of main time-stepping loop
    //----------------------------------------------------------------
    while ( !finished_time_stepping ) {
        // 0. Alter configuration setting if necessary.
	if ( (step/GlobalConfig.control_count)*GlobalConfig.control_count == step ) {
	    read_control_file(); // Reparse the time-step control parameters occasionally.
	}

	// 1. Set the size of the time step.
	// TODO: check CFL constrain across blocks.
	dt_global = GlobalConfig.dt_init; 

        // 2. Attempt a time step.
	// 2a.
	// explicit or implicit update of the convective terms.
	gasdynamic_explicit_increment_with_fixed_grid();
	// 2b. Recalculate all geometry if moving grid.
	// 2c. Increment because of viscous effects may be done
	//     separately to the convective terms.
        // 2d. Chemistry step. 
	//     Allow finite-rate evolution of species due
        //     to chemical reactions

        // 3. Update the time record and (occasionally) print status.
        ++step;
        output_just_written = false;
        history_just_written = false;
        if ( (step / GlobalConfig.print_count) * GlobalConfig.print_count == step ) {
            // Print the current time-stepping status.
	    writeln("Step= ", step, " sim_time=", sim_time);
	}

        // 4. (Occasionally) Write out an intermediate solution
        if ( (sim_time >= t_plot) && !output_just_written ) {
	    ++current_tindx;
	    ensure_directory_is_present("something");
	    foreach (ref myblk; myBlocks) {
		auto file_name = make_path_name!"flow"(myblk.id, current_tindx);
		myblk.write_solution(file_name, sim_time);
	    }
	    output_just_written = true;
            t_plot += GlobalConfig.dt_plot;
        }

        // 5. For steady-state approach, check the residuals for mass and energy.

	// 6. Spatial filter may be applied occasionally.

        // 7. Loop termination criteria:
        //    (1) reaching a maximum simulation time or target time
        //    (2) reaching a maximum number of steps
        //    (3) finding that the "halt_now" parameter has been set 
	//        in the control-parameter file.
        //        This provides a semi-interactive way to terminate the 
        //        simulation and save the data.
	//    (4) Exceeding a maximum number of wall-clock seconds.
	//    (5) Having the temperature at one of the control points exceed 
	//        the preset tolerance.  
	//        This is mainly for the radiation-coupled simulations.
	//    (-) Exceeding an allowable delta(f_rad) / f_rad_org factor
	//
	//    Note that the max_time and max_step control parameters can also
	//    be found in the control-parameter file (which may be edited
	//    while the code is running).
        if ( sim_time >= target_time ) {
            finished_time_stepping = true;
            if( GlobalConfig.verbosity_level >= 1)
		writeln("Integration stopped: reached maximum simulation time.");
        }
        if (step >= GlobalConfig.max_step) {
            finished_time_stepping = true;
            if (GlobalConfig.verbosity_level >= 1)
		writeln("Integration stopped: reached maximum number of steps.");
        }
        if (GlobalConfig.halt_now == 1) {
            finished_time_stepping = true;
            if (GlobalConfig.verbosity_level >= 1)
		writeln("Integration stopped: Halt set in control file.");
        }
	int now = 0; // time(NULL);
	if (maxWallClock > 0 && (to!int(now - wall_clock_at_start) > maxWallClock)) {
            finished_time_stepping = true;
            if (GlobalConfig.verbosity_level >= 1)
		writeln("Integration stopped: reached maximum wall-clock time.");
	}
    } // end while !finished_time_stepping

    if (GlobalConfig.verbosity_level > 0) writeln("Done integrate_in_time().");
    return sim_time;
} // end integrate_in_time()

void finalize_simulation(double sim_time)
{
    if (GlobalConfig.verbosity_level > 0) writeln("Finalize the simulation.");
    if (!output_just_written) {
	++current_tindx;
	ensure_directory_is_present("something");
	foreach (ref myblk; myBlocks) {
	    auto file_name = make_path_name!"flow"(myblk.id, current_tindx);
	    myblk.write_solution(file_name, sim_time);
	}
    }
    if (GlobalConfig.verbosity_level > 0) writeln("Done finalize_simulation.");
} // end finalize_simulation()

//----------------------------------------------------------------------------

void ensure_directory_is_present(string path)
{
    writeln("TODO implelemt ensure_directory_is_present");
}

void gasdynamic_explicit_increment_with_fixed_grid()
{
    writeln("TODO implement explicit increment with fixed grid");
}
