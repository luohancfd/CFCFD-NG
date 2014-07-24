/**
 * fvcore.d
 * Core definitions for finite-volume cells, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvcore;

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
uint[string] face_index; // initialized with a call to init_fvcore()

// Symbolic names for the time-stepping schemes used to update the gasdynamic eqn.
enum 
    euler_update = 0, 
    pc_update = 1,
    midpoint_update = 2, 
    classic_rk3_update = 3,
    tvd_rk3_update = 4,
    denman_rk3_update = 5;

int gasdynamic_update_scheme = pc_update;

string[] gasdynamic_update_scheme_name = ["euler", "predictor-corrector", "midpoint", 
					  "classic-rk3", "tvd-rk3", "denman-rk3"];
int[] number_of_stages_for_update_scheme = [1, 2, 2, 3, 3, 3];
uint[string] gasdynamic_update_scheme_index; // initialized with a call to init_fvcore()


// [TODO] think about the following...
enum
    copy_all_data = 0,
    copy_flow_data = 1,
    copy_grid_data = 2,
    copy_cell_lengths_only = 3;

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

// Symbolic names for the flavours of our flux_calculators.
enum
    flux_ausmdv = 0, // Wada and Liou's flux calculator AIAA Paper 94-0083
    flux_efm = 1, // Mike Macrossan's EFM flux calculation
    flux_ausm_plus_up = 2, // Liou's 2006 all-speed flux calculator
    flux_adaptive = 3, // EFM near shocks, AUSMDV otherwise
    flux_hlle = 4; // MHD HLLE approximate Riemann solver

string[] fluxcalc_name = [ "ausmdv", "efm", "ausm_plus_up", "adaptive", "hlle" ];
uint[string] fluxcalc_index; // initialized with a call to init_fluxcalc()


void init_fvcore()
{
    face_index = ["north":north, "east":east, "south":south,
		  "west":west, "top":top, "bottom":bottom];
    gasdynamic_update_scheme_index = ["euler":euler_update,
				      "pc":pc_update,
				      "midpoint":midpoint_update, 
				      "classic_rk3":classic_rk3_update,
				      "tvd_rk3":tvd_rk3_update,
				      "denman_rk3":denman_rk3_update];
    fluxcalc_index = ["ausmdv":flux_ausmdv,
		      "efm":flux_efm, 
		      "ausm_plus_up":flux_ausm_plus_up,
		      "adaptive":flux_adaptive,
		      "hlle":flux_hlle];
}
