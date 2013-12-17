#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include "../util/config_parser.hh"
#include "../math/no_fuss_linear_algebra_jc.hh"



#include "init.hh"
#include "adv.hh"
#include "conjugateHeatUtil.hh"
#include <stdio.h>
#include "e3con.hh"

using namespace std;

int main(){

	list_of_inputs inputs;
	list_of_vars vars;
	// list_of_flow_fluxes flow_fluxes;

	init_tr(inputs, vars);

	double time_elapsed = 0.0;
	int print_number = 0;
	vector<double> q_wall; q_wall.resize(inputs.N*inputs.M, 0.0);

	if (inputs.north == "udf_lookup" || inputs.south == "udf_lookup" || inputs.east == "udf_lookup" || inputs.west == "udf_lookup") {
		read_lookup_flux(inputs, vars);
	}
	
	old_file_delete_catcher(inputs);
	cout << "Starting Solution..." << endl;

	cout << "Initial Temperature Profile... " << endl;
	wall_printer(vars.temps, inputs);


	while (time_elapsed < inputs.t_max) {
		update_coeff_matrix(vars, inputs);
		flux_calc(inputs, vars, q_wall);
		increment_in_time(vars, inputs);
		// wall_printer(vars.fluxes, inputs);
		// wall_printer(vars.temps, inputs);

		double check = fmod(time_elapsed, inputs.dt_plot);
		if (check < 1e-10 || check-inputs.dt < 1e-10 || time_elapsed == 0) { // This is a dodgy fix, fmod(1.1, 0.1) is returning 0.1???)
			// write_row(vars, inputs, 5);
			// write_column(vars, inputs, 1);
			write_soln(vars, inputs, time_elapsed, print_number);
			cout << "Solution successfully written after " << time_elapsed << " seconds." << endl;
			print_number++;
		}

		check = fmod(time_elapsed, inputs.dt_hist);

		if (inputs.flag_hist == 1 && ( check < 1e-10 || time_elapsed == 0 )) {
			cout << "writing history..." << endl;
			write_hist(vars, inputs, time_elapsed);
		}
		time_elapsed+=inputs.dt;
	}
	system("gnuplot \"../util/temp_hist_plot.p\"");





}

