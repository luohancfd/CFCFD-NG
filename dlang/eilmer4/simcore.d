/** simcore.d
 * Eilmer4 compressible-flow simulation code, core coordination functions.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

module simcore;

import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.algorithm;
import std.datetime;
import std.parallelism;

import fileutil;
import geom;
import gas;
import fvcore;
import globalconfig;
import readconfig;
import globaldata;
import flowstate;
import sblock;
import ssolidblock;
import solidfvcore;
import solidprops;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;

// State data for simulation.
// Needs to be seen by all of the coordination functions.
shared static int current_tindx = 0;
shared static double sim_time = 0;  // present simulation time, tracked by code
shared static int step = 0;
shared static double dt_global;     // simulation time step determined by code
shared static double dt_allow;      // allowable global time step determined by code
shared static double t_plot;        // time to write next soln
shared static bool output_just_written = true;
shared static double t_history;     // time to write next sample
shared static bool history_just_written = true;

 // For working how long the simulation has been running.
static SysTime wall_clock_start;

//----------------------------------------------------------------------------

double init_simulation(int tindx)
{
    if (GlobalConfig.verbosity_level > 0) writeln("Begin init_simulation()...");
    wall_clock_start = Clock.currTime();
    read_config_file();  // most of the configuration is in here
    read_control_file(); // some of the configuration is in here
    current_tindx = tindx;
    shared double sim_time;
    auto job_name = GlobalConfig.base_file_name;
    writeln("Parallel version running on ", totalCPUs, " CPUs");
    foreach (ref myblk; parallel(myBlocks,1)) {
	myblk.assemble_arrays();
	myblk.bind_interfaces_and_vertices_to_cells();
	myblk.bind_vertices_and_cells_to_interfaces();
	writeln("myblk=", myblk);
	myblk.read_grid(make_file_name!"grid"(job_name, myblk.id, tindx), 0);
	// I don't mind if blocks write over sim_time.  
	// They should all have the same value for it.
	sim_time = myblk.read_solution(make_file_name!"flow"(job_name, myblk.id, tindx));
    }
    foreach (ref myblk; parallel(myBlocks,1)) {
	myblk.compute_primary_cell_geometric_data(0);
	myblk.compute_distance_to_nearest_wall_for_all_cells(0);
	myblk.compute_secondary_cell_geometric_data(0);
	myblk.identify_reaction_zones(0);
	myblk.identify_turbulent_zones(0);
	myblk.set_grid_velocities(sim_time);
	foreach (ref cell; myblk.active_cells) {
	    cell.encode_conserved(0, 0, myblk.omegaz);
	    // Even though the following call appears redundant at this point,
	    // fills in some gas properties such as Prandtl number that is
	    // needed for both the cfd_check and the BLomax turbulence model.
	    cell.decode_conserved(0, 0, myblk.omegaz, myblk.gmodel);
	}
	myblk.set_cell_dt_chem(-1.0);
    }
    foreach (ref mySolidBlk; mySolidBlocks) {
	mySolidBlk.assembleArrays();
	mySolidBlk.bindFacesAndVerticesToCells();
	writeln("mySolidBlk= ", mySolidBlk);
	mySolidBlk.readGrid(make_file_name!"solid-grid"(job_name, mySolidBlk.id, tindx));
	mySolidBlk.readSolution(make_file_name!"solid"(job_name, mySolidBlk.id, tindx));
    }
    foreach (ref mySolidBlk; mySolidBlocks) {
	mySolidBlk.computePrimaryCellGeometricData();
	mySolidBlk.computeSecondaryCellGeometricData();
    }
    // [TODO] Is it appropriate to disable this?
    //exchange_shared_boundary_data();
    if (GlobalConfig.verbosity_level > 0) writeln("Done init_simulation().");
    return sim_time;
} // end init_simulation()

void update_times_file()
{
    auto writer = appender!string();
    formattedWrite(writer, "%04d %e %e\n", current_tindx, sim_time, dt_global);
    append(GlobalConfig.base_file_name ~ ".times", writer.data);
}

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
    // When starting a new calculation,
    // set the global time step to the initial value.
    dt_global = GlobalConfig.dt_init; 
    shared bool do_cfl_check_now = false;
    // Normally, we can terminate upon either reaching 
    // a maximum time or upon reaching a maximum iteration count.
    shared bool finished_time_stepping = 
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

	// 1. Set the size of the time step to be the minimum allowed for any active block.
	if (!GlobalConfig.fixed_time_step && 
	    (step/GlobalConfig.cfl_count)*GlobalConfig.cfl_count == step) {
	    // Check occasionally 
	    do_cfl_check_now = true;
	} // end if step == 0
	if (do_cfl_check_now) {
	    // Adjust the time step  
	    shared double dt_allow = 1.0e9; // Start with too large a guess to ensure it is replaced.
	    // [TODO] for parallel, need to reduce across myBlocks
	    foreach (ref myblk; myBlocks) {
		if (!myblk.active) continue;
		dt_allow = min(dt_allow, myblk.determine_time_step_size(dt_global)); 
	    }
	    // Change the actual time step, as needed.
	    if (dt_allow <= dt_global) {
		// If we need to reduce the time step, do it immediately.
		dt_global = dt_allow;
	    } else {
		// Make the transitions to larger time steps gentle.
		dt_global = min(dt_global*1.5, dt_allow);
		// The user may supply, explicitly, a maximum time-step size.
		dt_global = min(dt_global, GlobalConfig.dt_max);
	    }
	    do_cfl_check_now = false;  // we have done our check for now
	} // end if do_cfl_check_now 

        // 2. Attempt a time step.
	foreach (ref myblk; parallel(myBlocks,1)) {
	    myblk.set_grid_velocities(sim_time); 
	    // [TODO] We will need to attend to moving grid properly.
	}
	// 2a.
	// explicit or implicit update of the convective terms.
	gasdynamic_explicit_increment_with_fixed_grid();
	// 2b. Recalculate all geometry if moving grid.
	// 2c. Increment because of viscous effects may be done
	//     separately to the convective terms.
        // 2d. Chemistry step. 
	if ( GlobalConfig.reacting ) {
	    foreach (blk; parallel(myBlocks,1)) {
		if (!blk.active) continue;
		foreach ( i, cell; blk.active_cells) {
		    cell.chemical_increment(dt_global, 300.0, blk.gmodel, blk.reaction_update);
		}
	    }
	}
	// 2e.
	// Update the temperature values in the solid domain.
	solidDomainExplicitIncrement();
        // 3. Update the time record and (occasionally) print status.
        step = step + 1;
        output_just_written = false;
        history_just_written = false;
        if ( (step / GlobalConfig.print_count) * GlobalConfig.print_count == step ) {
            // Print the current time-stepping status.
	    auto writer = appender!string();
	    formattedWrite(writer, "Step=%7d t=%10.3e dt=%10.3e ", step, sim_time, dt_global);
	    long wall_clock_elapsed = (Clock.currTime() - wall_clock_start).total!"seconds"();
	    double wall_clock_per_step = to!double(wall_clock_elapsed) / step;
	    double WCtFT = (GlobalConfig.max_time - sim_time) / dt_global * wall_clock_per_step;
	    double WCtMS = (GlobalConfig.max_step - step) * wall_clock_per_step;
	    formattedWrite(writer, "WC=%d WCtFT=%.1f WCtMS=%.1f", 
			   wall_clock_elapsed, WCtFT, WCtMS);
	    writeln(writer.data);
	}

        // 4. (Occasionally) Write out an intermediate solution
        if ( (sim_time >= t_plot) && !output_just_written ) {
	    current_tindx = current_tindx + 1;
	    ensure_directory_is_present(make_path_name!"flow"(current_tindx));
	    auto job_name = GlobalConfig.base_file_name;
	    foreach (ref myblk; parallel(myBlocks,1)) {
		auto file_name = make_file_name!"flow"(job_name, myblk.id, current_tindx);
		myblk.write_solution(file_name, sim_time);
	    }
	    ensure_directory_is_present(make_path_name!"solid"(current_tindx));
	    foreach (ref mySolidBlk; mySolidBlocks) {
		auto fileName = make_file_name!"solid"(job_name, mySolidBlk.id, current_tindx);
		mySolidBlk.writeSolution(fileName, sim_time);
	    }
	    update_times_file();
	    output_just_written = true;
            t_plot = t_plot + GlobalConfig.dt_plot;
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
	auto wall_clock_elapsed = (Clock.currTime() - wall_clock_start).total!"seconds"();
	if (maxWallClock > 0 && (wall_clock_elapsed > maxWallClock)) {
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
	current_tindx = current_tindx + 1;
	ensure_directory_is_present(make_path_name!"flow"(current_tindx));
	auto job_name = GlobalConfig.base_file_name;
	foreach (ref myblk; myBlocks) {
	    auto file_name = make_file_name!"flow"(job_name, myblk.id, current_tindx);
	    myblk.write_solution(file_name, sim_time);
	}
	update_times_file();
    }
    writeln("Step= ", step, " final-t= ", sim_time);
    if (GlobalConfig.verbosity_level > 0) writeln("Done finalize_simulation.");
} // end finalize_simulation()

//----------------------------------------------------------------------------

void gasdynamic_explicit_increment_with_fixed_grid()
{
    shared double t0 = sim_time;
    shared bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega) &&
	!GlobalConfig.separate_update_for_k_omega_source;
    // Set the time-step coefficients for the stages of the update scheme.
    shared double c2 = 1.0; // default for predictor-corrector update
    shared double c3 = 1.0; // default for predictor-corrector update
    final switch ( GlobalConfig.gasdynamic_update_scheme ) {
    case GasdynamicUpdate.euler:
    case GasdynamicUpdate.pc: c2 = 1.0; c3 = 1.0; break;
    case GasdynamicUpdate.midpoint: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.classic_rk3: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.tvd_rk3: c2 = 1.0; c3 = 0.5; break;
    case GasdynamicUpdate.denman_rk3: c2 = 1.0; c3 = 0.5; break; 
    }
    // Preparation for the predictor-stage of inviscid gas-dynamic flow update.
    foreach (blk; parallel(myBlocks,1)) {
	if (!blk.active) continue;
	blk.clear_fluxes_of_conserved_quantities();
	foreach (cell; blk.active_cells) cell.clear_source_vector();
    }
    // First-stage of gas-dynamic update.
    shared int ftl = 0; // time-level within the overall convective-update
    shared int gtl = 0; // grid time-level remains at zero for the non-moving grid
    foreach (blk; myBlocks) { // # [TODO] [FIXME] parallel caused segmentation fault
	if (!blk.active) continue;
	blk.applyPreReconAction(sim_time, gtl, ftl);
    }
    // We've put this detector step here because it needs the ghost-cell data
    // to be current, as it should be just after a call to apply_convective_bc().
    if (GlobalConfig.flux_calculator == FluxCalculator.adaptive) {
	foreach (blk; parallel(myBlocks,1)) { blk.detect_shock_points(); }
    }
    foreach (blk; parallel(myBlocks,1)) {
	if (!blk.active) continue;
	blk.convective_flux();
	if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
	    blk.applyPreSpatialDerivAction(sim_time, gtl, ftl);
	    blk.flow_property_derivatives(gtl); 
	    blk.estimate_turbulence_viscosity();
	    blk.viscous_flux();
	} // end if viscous
	// [TODO] [FIXME] push the fine loop here into the Block class.
	// The cone20 run time has gone from 21 to 28 seconds, maybe because of
	// many memory barriers on the shared variables like gtl, ftl, with_k_omega.
	int local_ftl = ftl;
	int local_gtl = gtl;
	bool local_with_k_omega = with_k_omega;
	double local_dt_global = dt_global;
	double local_sim_time = sim_time;
	foreach (cell; blk.active_cells) {
	    cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
	    if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
		cell.add_viscous_source_vector(local_with_k_omega);
	    }
	    if (blk.udf_source_terms) {
		addUDFSourceTermsToCell(blk.udf_source_terms, cell, local_gtl, local_sim_time, blk.gmodel);
	    }
	    cell.time_derivatives(local_gtl, local_ftl, GlobalConfig.dimensions, local_with_k_omega);
	    bool force_euler = false;
	    cell.stage_1_update_for_flow_on_fixed_grid(local_dt_global, force_euler, local_with_k_omega);
	    cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz, blk.gmodel);
	} // end foreach cell
    } // end foreach blk
    //
    if ( number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) >= 2 ) {
	// Preparation for second-stage of gas-dynamic update.
	sim_time = t0 + c2 * dt_global;
	foreach (blk; parallel(myBlocks,1)) {
	    if (!blk.active ) continue;
	    blk.clear_fluxes_of_conserved_quantities();
	    foreach (cell; blk.active_cells) cell.clear_source_vector();
	}
	// Second stage of gas-dynamic update.
	ftl = 1;
	// We are relying on exchangine boundary data
	// as a pre-reconstruction activity.
	// NOTE FOR PJ: This replaces exchange. Good idea or no? 
	// Reply from PJ: YES
	foreach (blk; myBlocks) { // # [TODO] [FIXME] parallel caused segmentation fault
	    if (!blk.active) continue;
	    blk.applyPreReconAction(sim_time, gtl, ftl);
	}
	foreach (blk; parallel(myBlocks,1)) {
	    if (!blk.active) continue;
	    blk.convective_flux();
	    if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
		blk.applyPreSpatialDerivAction(sim_time, gtl, ftl);
		blk.flow_property_derivatives(gtl); 
		blk.estimate_turbulence_viscosity();
		blk.viscous_flux();
	    } // end if viscous
	    // [TODO] [FIXME] push the fine loop here into the Block class.
	    // The cone20 run time has gone from 21 to 28 seconds, maybe because of
	    // many memory barriers on the shared variables like gtl, ftl, with_k_omega.
	    int local_ftl = ftl;
	    int local_gtl = gtl;
	    bool local_with_k_omega = with_k_omega;
	    double local_dt_global = dt_global;
	    double local_sim_time = sim_time;
	    foreach (cell; blk.active_cells) {
		cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
		if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
		    cell.add_viscous_source_vector(local_with_k_omega);
		}
		if (blk.udf_source_terms) {
		    addUDFSourceTermsToCell(blk.udf_source_terms, cell, local_gtl, local_sim_time, blk.gmodel);
		}
		cell.time_derivatives(local_gtl, local_ftl, GlobalConfig.dimensions, local_with_k_omega);
		cell.stage_2_update_for_flow_on_fixed_grid(dt_global, local_with_k_omega);
		cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz, blk.gmodel);
	    } // end foreach cell
	} // end foreach blk
    } // end if number_of_stages_for_update_scheme >= 2 
    //
    if ( number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) >= 3 ) {
	// Preparation for third stage of gasdynamic update.
	sim_time = t0 + c3 * dt_global;
	foreach (blk; parallel(myBlocks,1)) {
	    if (!blk.active ) continue;
	    blk.clear_fluxes_of_conserved_quantities();
	    foreach (cell; blk.active_cells) cell.clear_source_vector();
	}
	ftl = 2;
	// We are relying on exchangine boundary data
	// as a pre-reconstruction activity.
	foreach (blk; myBlocks) { // [TODO] parallel
	    if (!blk.active) continue;
	    blk.applyPreReconAction(sim_time, gtl, ftl);
	}
	foreach (blk; parallel(myBlocks,1)) {
	    if (!blk.active) continue;
	    blk.convective_flux();
	    if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
		blk.applyPreSpatialDerivAction(sim_time, gtl, ftl);
		blk.flow_property_derivatives(gtl); 
		blk.estimate_turbulence_viscosity();
		blk.viscous_flux();
	    } // end if viscous
	    // [TODO] [FIXME] push the fine loop here into the Block class.
	    // The cone20 run time has gone from 21 to 28 seconds, maybe because of
	    // many memory barriers on the shared variables like gtl, ftl, with_k_omega.
	    int local_ftl = ftl;
	    int local_gtl = gtl;
	    bool local_with_k_omega = with_k_omega;
	    double local_dt_global = dt_global;
	    double local_sim_time = sim_time;
	    foreach (cell; blk.active_cells) {
		cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
		if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
		    cell.add_viscous_source_vector(local_with_k_omega);
		}
		if (blk.udf_source_terms) {
		    addUDFSourceTermsToCell(blk.udf_source_terms, cell, local_gtl, local_sim_time, blk.gmodel);
		}
		cell.time_derivatives(local_gtl, local_ftl, GlobalConfig.dimensions, local_with_k_omega);
		bool force_euler = false;
		cell.stage_2_update_for_flow_on_fixed_grid(local_dt_global, with_k_omega);
		cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz, blk.gmodel);
	    } // end foreach cell
	} // end foreach blk
    } // end if number_of_stages_for_update_scheme >= 3
    //
    // Get the end conserved data into U[0] for next step.
    size_t end_indx = 2;
    final switch (GlobalConfig.gasdynamic_update_scheme) {
    case GasdynamicUpdate.euler: end_indx = 1; break;
    case GasdynamicUpdate.pc: 
    case GasdynamicUpdate.midpoint: end_indx = 2; break;
    case GasdynamicUpdate.tvd_rk3:
    case GasdynamicUpdate.classic_rk3:
    case GasdynamicUpdate.denman_rk3: end_indx = 3; break;
    }
    foreach (blk; parallel(myBlocks,1)) {
	if (!blk.active) continue;
	foreach (cell; blk.active_cells) {
	    swap(cell.U[0], cell.U[end_indx]);
	}
    } // end foreach blk
    //
    // Finally, update the globally know simulation time for the whole step.
    sim_time = t0 + dt_global;
} // end gasdynamic_explicit_increment_with_fixed_grid()

void solidDomainExplicitIncrement()
{
    double t0 = sim_time;
    // Set the time-step coefficients for the stages of the update scheme.
    // Below are defaults for predictor-corrector update.
    double c2 = 1.0; 
    double c3 = 1.0;
    int tLevel = 0;
    // Apply boundary condition
    foreach (blk; mySolidBlocks) {
	if (!blk.active) continue;
	blk.applyPreSpatialDerivAction(sim_time, tLevel);
    }
    // First stage of update
    foreach (blk; mySolidBlocks) {
	if (!blk.active) continue;
	blk.computeSpatialDerivatives(tLevel);
	blk.applyPostFluxAction(sim_time, tLevel);
	blk.computeFluxes();
	foreach (cell; blk.activeCells) {
	    if (blk.udfSourceTerms) {
		addUDFSourceTermsToSolidCell(blk.udfSourceTerms, cell, sim_time);
	    }
	    cell.timeDerivatives(0, GlobalConfig.dimensions);
	    cell.stage1Update(dt_global);
	    cell.T[1] = updateTemperature(blk.sp, cell.e[1]);
	} // end foreach cell
    } // end foreach blk

    if ( numberOfStagesForUpdateScheme(GlobalConfig.solidDomainUpdateScheme) >= 2 ) {
	// Second stage of solid domain update
	tLevel = 1;
	// Apply boundary condition
	foreach (blk; mySolidBlocks) {
	    if (!blk.active) continue;
	    blk.applyPreSpatialDerivAction(sim_time, tLevel);
	}
	// Second stage of update
	foreach (blk; mySolidBlocks) {
	    if (!blk.active) continue;
	    blk.computeSpatialDerivatives(tLevel);
	    blk.applyPostFluxAction(sim_time, tLevel);
	    blk.computeFluxes();
	    foreach (cell; blk.activeCells) {
		if (blk.udfSourceTerms) {
		    addUDFSourceTermsToSolidCell(blk.udfSourceTerms, cell, sim_time);
		}
		cell.timeDerivatives(1, GlobalConfig.dimensions);
		cell.stage2Update(dt_global);
		cell.T[2] = updateTemperature(blk.sp, cell.e[2]);
	    } // end foreach cell
	} // end foreach blk
    } // end if numberOfStagesForUpdateScheme >= 2

    // Get the final bit of conserved data in U[0] for next step
    size_t endIdx = 2;
    final switch (GlobalConfig.solidDomainUpdateScheme) {
    case SolidDomainUpdate.euler: endIdx = 1; break;
    case SolidDomainUpdate.pc: endIdx = 2; break;
    }
    foreach (blk; mySolidBlocks) {
	if (!blk.active) continue;
	foreach (cell; blk.activeCells) {
	    cell.e[0] = cell.e[endIdx];
	    cell.T[0] = cell.T[endIdx];
	} 
    } // end foreach block
} // end solidDomainExplicitIncrement()
