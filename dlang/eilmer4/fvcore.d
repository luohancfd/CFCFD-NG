/**
 * fvcore.d
 * Core definitions for finite-volume cells, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvcore;

import std.conv;

// Symbolic names and indices for the cells' faces.
// The names of the faces of the structured-grid blocks will be the same.
enum
    north = 0,
    east = 1,
    south = 2,
    west = 3,
    top = 4,
    bottom = 5;

string[] face_name = [ "north", "east", "south", "west", "top", "bottom" ];
uint face_index(string name)
{
    switch ( name ) {
    case "north": return north;
    case "east": return east;
    case "south": return south;
    case "west": return west;
    case "top": return top;
    case "bottom": return bottom;
    default:
	throw new Error(text("Invalid face name:", name));
    }
} // end face_index

// Symbolic names for the time-stepping schemes used to update the gasdynamic eqn.
enum GasdynamicUpdate {
    euler, 
    pc,
    midpoint, 
    classic_rk3,
    tvd_rk3,
    denman_rk3
}

GasdynamicUpdate gasdynamic_update_scheme = GasdynamicUpdate.pc;

string gasdynamic_update_scheme_name(GasdynamicUpdate gdut)
{
    final switch ( gdut ) {
    case GasdynamicUpdate.euler: return "euler";
    case GasdynamicUpdate.pc: return "predictor-corrector";
    case GasdynamicUpdate.midpoint: return "midpoint";
    case GasdynamicUpdate.classic_rk3: return "classic-rk3";
    case GasdynamicUpdate.tvd_rk3: return "tvd-rk3";
    case GasdynamicUpdate.denman_rk3: return "denman-rk3";
    }
} // end gasdynamic_update_scheme_name()

int number_of_stages_for_update_scheme(GasdynamicUpdate gdut)
{
    final switch ( gdut ) {
    case GasdynamicUpdate.euler: return 1;
    case GasdynamicUpdate.pc: return 2;
    case GasdynamicUpdate.midpoint: return 2;
    case GasdynamicUpdate.classic_rk3: return 3;
    case GasdynamicUpdate.tvd_rk3: return 3;
    case GasdynamicUpdate.denman_rk3: return 3;
    }
} // end number_of_stages_for_update_scheme()

GasdynamicUpdate update_scheme_from_name(string name)
{
    switch ( name ) {
    case "euler": return GasdynamicUpdate.euler;
    case "pc": return GasdynamicUpdate.pc;
    case "predictor_corrector": return GasdynamicUpdate.pc;
    case "predictor-corrector": return GasdynamicUpdate.pc;
    case "midpoint": return GasdynamicUpdate.midpoint;
    case "classic_rk3": return GasdynamicUpdate.classic_rk3;
    case "classic-rk3": return GasdynamicUpdate.classic_rk3;
    case "tvd_rk3": return GasdynamicUpdate.tvd_rk3;
    case "tvd-rk3": return GasdynamicUpdate.tvd_rk3;
    case "denman_rk3": return GasdynamicUpdate.denman_rk3;
    case "denman-rk3": return GasdynamicUpdate.denman_rk3;
    default:
	throw new Error(text("Invalid gasdynamic update scheme name:", name));
    }
}  // end scheme_from_name()


// [TODO] think about the following...
enum CopyDataOption { all, flow, grid, cell_lengths_only }

// n_time_levels is used to size the time-derivative vectors.
enum
    n_time_levels = 4,
    n_interfaces_per_cell = 6,
    n_vertex_per_cell = 8;

// Minimum values for turbulent kinetic energy (m^2/s^2) and frequency (1/s)
// for applying limiters in the k-omega model.
enum
    small_tke = 0.1,
    small_omega = 1.0;

// Symbolic names for the types of flow-data reconstruction.
enum InterpolateOption { pt, rhoe, rhop, rhot }

string thermo_interpolator_name(InterpolateOption i)
{
    final switch ( i ) {
    case InterpolateOption.pt: return "pT";
    case InterpolateOption.rhoe: return "rhoe";
    case InterpolateOption.rhop: return "rhop";
    case InterpolateOption.rhot: return "rhoT";
    }
} // end thermo_interpolator_name()

InterpolateOption thermo_interpolator_from_name(string name)
{
    switch ( name ) {
    case "pT": return InterpolateOption.pt;
    case "pt": return InterpolateOption.pt;
    case "rhoe": return InterpolateOption.rhoe;
    case "rhop": return InterpolateOption.rhop;
    case "rhoT": return InterpolateOption.rhot;
    case "rhot": return InterpolateOption.rhot;
    default: return InterpolateOption.rhoe;
    }
} // end thermo_interpolator_from_name()

// Symbolic names for the flavours of our flux_calculators.
enum FluxCalculator {
    ausmdv, // Wada and Liou's flux calculator AIAA Paper 94-0083
    efm, // Mike Macrossan's EFM flux calculation
    ausm_plus_up, // Liou's 2006 all-speed flux calculator
    adaptive, // EFM near shocks, AUSMDV otherwise
    hlle // MHD HLLE approximate Riemann solver
}

string fluxcalc_name(FluxCalculator fcalc)
{
    final switch ( fcalc ) {
    case FluxCalculator.ausmdv: return "ausmdv";
    case FluxCalculator.efm: return "efm";
    case FluxCalculator.ausm_plus_up: return "ausm_plus_up";
    case FluxCalculator.adaptive: return "adaptive";
    case FluxCalculator.hlle: return "hlle";
    }
}

FluxCalculator fluxcalc_from_name(string name)
{
    switch ( name ) {
    case "ausmdv": return FluxCalculator.ausmdv;
    case "efm": return FluxCalculator.efm;
    case "ausm_plus_up": return FluxCalculator.ausm_plus_up;
    case "adaptive": return FluxCalculator.adaptive;
    case "hlle": return FluxCalculator.hlle;
    default:
	throw new Error(text("Invalid flux calculator name:", name));
    }
}
