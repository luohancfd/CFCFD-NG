/** readconfig.d
 * Eilmer4 compressible-flow simulation code, reading of JSON config and control files.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

import std.stdio;
import std.json;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;

import json_helper;
import geom;
import gas;
import kinetics;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import sblock;
import bc;

void read_config_file()
{
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

    // Now that we have parsed JSON data, dip into it to get config values.
    //
    GlobalConfig.title = jsonData["title"].str;
    GlobalConfig.gas_model_file = jsonData["gas_model_file"].str;
    GlobalConfig.gmodel = init_gas_model(GlobalConfig.gas_model_file);
    GlobalConfig.dimensions = getJSONint(jsonData, "dimensions", 2);
    GlobalConfig.axisymmetric = getJSONbool(jsonData, "axisymmetric", false);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  title: ", GlobalConfig.title);
	writeln("  gas_model_file: ", GlobalConfig.gas_model_file);
	writeln("  dimensions: ", GlobalConfig.dimensions);
	writeln("  axisymmetric: ", GlobalConfig.axisymmetric);
    }

    // Parameters controlling convective update
    //
    // GlobalConfig.gasdynamic_update_scheme set from .control file.
    GlobalConfig.adjust_invalid_cell_data =
	getJSONbool(jsonData, "adjust_invalid_cell_data", false);
    GlobalConfig.max_invalid_cells = getJSONint(jsonData, "max_invalid_cells", 0);
    try {
	string name = jsonData["thermo_interpolator"].str;
	GlobalConfig.thermo_interpolator = thermo_interpolator_from_name(name);
    } catch (Exception e) {
	GlobalConfig.thermo_interpolator = InterpolateOption.rhoe;
    }
    GlobalConfig.apply_limiter = getJSONbool(jsonData, "apply_limiter", true);
    GlobalConfig.extrema_clipping = getJSONbool(jsonData, "extreme_clipping", true);
    GlobalConfig.interpolate_in_local_frame = 
	getJSONbool(jsonData, "interpolate_in_local_frame", true);
    try {
	string name = jsonData["flux_calculator"].str;
	GlobalConfig.flux_calculator = flux_calculator_from_name(name);
    } catch (Exception e) {
	GlobalConfig.flux_calculator = FluxCalculator.adaptive;
    }
    GlobalConfig.shear_tolerance = getJSONdouble(jsonData, "shear_tolerance", 0.20);
    GlobalConfig.M_inf = getJSONdouble(jsonData, "M_inf", 0.01);
    GlobalConfig.compression_tolerance = 
	getJSONdouble(jsonData, "compression_tolerance", -0.30);
    GlobalConfig.moving_grid = getJSONbool(jsonData, "moving_grid", false);
    GlobalConfig.write_vertex_velocities = 
	getJSONbool(jsonData, "write_vertex_velocities", false);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  adjust_invalid_cell_data: ", GlobalConfig.adjust_invalid_cell_data);
	writeln("  max_invalid_cells: ", GlobalConfig.max_invalid_cells);
	writeln("  thermo_interpolator: ",
		thermo_interpolator_name(GlobalConfig.thermo_interpolator));
	writeln("  apply_limiter: ", GlobalConfig.apply_limiter);
	writeln("  extrema_clipping: ", GlobalConfig.extrema_clipping);
	writeln("  interpolate_in_local_frame: ", GlobalConfig.interpolate_in_local_frame);
	writeln("  flux_calculator: ", flux_calculator_name(GlobalConfig.flux_calculator));
	writeln("  shear_tolerance: ", GlobalConfig.shear_tolerance);
	writeln("  M_inf: ", GlobalConfig.M_inf);
	writeln("  compression_tolerance: ", GlobalConfig.compression_tolerance);
	writeln("  moving_grid: ", GlobalConfig.moving_grid);
	writeln("  write_vertex_velocities: ", GlobalConfig.write_vertex_velocities);
    }

    // Parameters controlling viscous/molecular transport
    //
    GlobalConfig.viscous = getJSONbool(jsonData, "viscous", false);
    GlobalConfig.viscous_delay = getJSONdouble(jsonData, "viscous_delay", 0.0);
    GlobalConfig.viscous_factor_increment = 
	getJSONdouble(jsonData, "viscous_factor_increment", 0.01);
    try {
	string name = jsonData["turbulence_model"].str;
	GlobalConfig.turbulence_model = turbulence_model_from_name(name);
    } catch (Exception e) {
	GlobalConfig.turbulence_model = TurbulenceModel.none;
    }
    GlobalConfig.turbulence_prandtl_number =
	getJSONdouble(jsonData, "turbulence_prandtl_number", 0.89);
    GlobalConfig.turbulence_schmidt_number =
	getJSONdouble(jsonData, "turbulence_schmidt_number", 0.75);
    GlobalConfig.max_mu_t_factor = getJSONdouble(jsonData, "max_mu_t_factor", 300.0);
    GlobalConfig.transient_mu_t_factor = getJSONdouble(jsonData, "transient_mu_t_factor", 1.0);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  viscous: ", GlobalConfig.viscous);
	writeln("  viscous_delay: ", GlobalConfig.viscous_delay);
	writeln("  viscous_factor_increment: ", GlobalConfig.viscous_factor_increment);
	writeln("  turbulence_model: ", turbulence_model_name(GlobalConfig.turbulence_model));
	writeln("  turbulence_prandtl_number: ", GlobalConfig.turbulence_prandtl_number);
	writeln("  turbulence_schmidt_number: ", GlobalConfig.turbulence_schmidt_number);
	writeln("  max_mu_t_factor: ", GlobalConfig.max_mu_t_factor);
	writeln("  transient_mu_t_factor: ", GlobalConfig.transient_mu_t_factor);
    }

    // Parameters controlling thermochemistry
    //
    GlobalConfig.reacting = getJSONbool(jsonData, "reacting", false);
    GlobalConfig.reactions_file = jsonData["reactions_file"].str;
    if ( GlobalConfig.reacting )
	GlobalConfig.reaction_update = new ReactionUpdateScheme(GlobalConfig.reactions_file, GlobalConfig.gmodel);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  reacting: ", GlobalConfig.reacting);
	writeln("  reactions_file: ", GlobalConfig.reactions_file);
    }

    // Parameters controlling other simulation options
    //
    GlobalConfig.control_count = getJSONint(jsonData, "control_count", 10);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  control_count: ", GlobalConfig.control_count);
    }
    // TODO -- still have other entries such as nheatzone, nreactionzone, ...

    // Now, configure blocks that make up the flow domain.
    //
    GlobalConfig.nBlocks = getJSONint(jsonData, "nblock", 0);
    if (GlobalConfig.verbosity_level > 1) { writeln("  nBlocks: ", GlobalConfig.nBlocks); }
    foreach (i; 0 .. GlobalConfig.nBlocks) {
	auto blk = new SBlock(i, jsonData["block_" ~ to!string(i)]);
	allBlocks ~= blk;
	myBlocks ~= blk; // Just make a copy, until we have to deal with MPI.
	if (GlobalConfig.verbosity_level > 1) {
	    writeln("  Block[", i, "]: ", myBlocks[i]);
	}
    }
} // end read_config_file()

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

    GlobalConfig.interpolation_order = getJSONint(jsonData, "interpolation_order", 2);
    try {
	string name = jsonData["gasdynamic_update_scheme"].str;
	GlobalConfig.gasdynamic_update_scheme = update_scheme_from_name(name);
    } catch (Exception e) {
	GlobalConfig.gasdynamic_update_scheme = GasdynamicUpdate.pc;
    }
    GlobalConfig.separate_update_for_viscous_terms =
	getJSONbool(jsonData, "separate_update_for_viscous_terms", false);
    GlobalConfig.max_step = getJSONint(jsonData, "max_step", 100);
    GlobalConfig.max_time = getJSONdouble(jsonData, "max_time", 1.0e-3);
    GlobalConfig.halt_now = getJSONint(jsonData, "halt_now", 0);
    GlobalConfig.print_count = getJSONint(jsonData, "print_count", 0);
    GlobalConfig.cfl_count = getJSONint(jsonData, "cfl_count", 0);
    GlobalConfig.dt_init = getJSONdouble(jsonData, "dt_init", 1.0e-6);
    GlobalConfig.dt_max = getJSONdouble(jsonData, "dt_max", 1.0-3);
    GlobalConfig.cfl_value = getJSONdouble(jsonData, "cfl_value", 0.5);
    GlobalConfig.stringent_cfl = getJSONbool(jsonData, "stringent_cfl", false);
    GlobalConfig.fixed_time_step = getJSONbool(jsonData, "fixed_time_step", false);
    GlobalConfig.dt_reduction_factor = getJSONdouble(jsonData, "dt_reduction_factor", 0.2);
    GlobalConfig.write_at_step = getJSONint(jsonData, "write_at_step", 0);
    GlobalConfig.dt_plot = getJSONdouble(jsonData, "dt_plot", 1.0e-3);
    GlobalConfig.dt_history = getJSONdouble(jsonData, "dt_history", 1.0e-3);
    if (GlobalConfig.verbosity_level > 1) {
	writeln("  interpolation_order: ", GlobalConfig.interpolation_order);
	writeln("  gasdynamic_update_scheme: ",
		gasdynamic_update_scheme_name(GlobalConfig.gasdynamic_update_scheme));
	writeln("  separate_update_for_viscous_terms: ", GlobalConfig.separate_update_for_viscous_terms);
	writeln("  max_step: ", GlobalConfig.max_step);
	writeln("  max_time: ", GlobalConfig.max_time);
	writeln("  halt_now: ", GlobalConfig.halt_now);
	writeln("  print_count: ", GlobalConfig.print_count);
	writeln("  cfl_count: ", GlobalConfig.cfl_count);
	writeln("  dt_init: ", GlobalConfig.dt_init);
	writeln("  dt_max: ", GlobalConfig.dt_max);
	writeln("  cfl_value: ", GlobalConfig.cfl_value);
	writeln("  stringent_cfl: ", GlobalConfig.stringent_cfl);
	writeln("  dt_reduction_factor: ", GlobalConfig.dt_reduction_factor);
	writeln("  fixed_time_step: ", GlobalConfig.fixed_time_step);
	writeln("  write_at_step: ", GlobalConfig.write_at_step);
	writeln("  dt_plot: ", GlobalConfig.dt_plot);
	writeln("  dt_history: ", GlobalConfig.dt_history);
    }
} // end read_control_file()
