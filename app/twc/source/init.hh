#ifndef TWC_INIT_HH
#define TWC_INIT_HH

#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"

using namespace std;

// Try writing a new structure to replicate the "inputs" Lua table that Zac used in Lua.
// The reason for this is so you don't have to use a heap of arguments to the grab_config and potentially further functions.

struct list_of_inputs {
	double stefanBoltzmann, D, L, M, N, y_h, k_11, k_12, k_21, k_22, dx, dy, rho_w, c_w, dt, dt_plot, t_max, Te_n, Te_s, Te_e, Te_w, he_n, he_s, he_e, he_w, emis_n, emis_s, emis_e, emis_w, dt_hist, T_init;
	int axi, flag_hist, nn_hist;
	string north, south, east, west;
	
};



// This is a new structure which contains the variables which involved in the update equation.  Separated from list_of_inputs
// structure for convenience.


struct list_of_vars {
	valarray<double> temps, fluxes, x_vals, y_vals, udf_temp_interp, udf_flux_interp;
	Valmatrix A, B, C;

};


// Functions
void grab_config(string fname, list_of_inputs &inputs);
int init_tr( list_of_inputs &inputs, list_of_vars &vars);
int axisym_checker(list_of_inputs &inputs);
int write_geom(list_of_inputs &inputs, list_of_vars &vars);
int init_coeff_matrix(list_of_inputs &inputs, list_of_vars &vars);
int init_tempprof(list_of_inputs &inputs, list_of_vars &vars);
int ar_checker(list_of_inputs &inputs);
int init_with_eilmer3(list_of_inputs &inputs, list_of_vars &vars);

#endif
