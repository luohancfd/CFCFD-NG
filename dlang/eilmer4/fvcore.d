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

string[] faceName = [ "north", "east", "south", "west", "top", "bottom" ];
uint[string] faceIndex;

enum 
    euler_update = 0, 
    pc_update = 1,
    midpoint_update = 2, 
    classic_rk3_update = 3,
    tvd_rk3_update = 4,
    denman_rk3_update = 5;

string[] updateSchemeName = [ "euler", "pc", "midpoint", 
			      "classic_rk3", "tvd_rk3", "denman_rk3" ];
uint[string] updateSchemeIndex;

enum
    copy_all_cell_data = 0,
    copy_flow_state_only = 1,
    copy_cell_lengths_only = 2;

// n_time_levels is used to size the time-derivative vectors.
enum
    n_time_levels = 4,
    n_interfaces_per_cell = 6,
    n_vertex_per_cell = 8;

// Minimum values for turbulent kinetic energy (m^2/s^2) and frequency (1/s)
// for applying limiters in the k-omega model.
enum
    minimum_tke = 0.1,
    minimum_omega = 1.0;


void init_fvcore()
{
    faceIndex = [ "north":north, "east":east, "south":south,
		  "west":west, "top":top, "bottom":bottom ];
    updateSchemeIndex = [ "euler":euler_update,
			  "pc":pc_update,
			  "midpoint":midpoint_update, 
			  "classic_rk3":classic_rk3_update,
			  "tvd_rk3":tvd_rk3_update,
			  "denman_rk3":denman_rk3_update ];
}
