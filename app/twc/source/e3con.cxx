#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/util/source/useful.h"
#include "init.hh"
#include "adv.hh"
#include "conjugateHeatUtil.hh"
#include "e3con.hh"

int grab_config_from_file(string fname, Wall_model &wm) {
	grab_config(fname, wm.inputs);
	return 0;
}


Wall_model* initialise_wall_model(string fname, double dt_plot) {

	Wall_model *wm = new Wall_model;
	
	//Read configuration file and store stuff in Wall_model.
	
	int flag = grab_config_from_file(fname, *wm);
	if ( flag != 0 ) {
		cerr << "There is a problem with reading the config file. ";
		cerr << "It is most likely that the location specified for the config file (fname) is incorrect. ";
		cerr << "Check file location.  Bailing out... " << endl;
		exit(MISMATCHED_DIMENSIONS);
    }	
	
	(*wm).inputs.dt_plot = dt_plot;

	return wm;
}

int sv(Wall_model &wm) {
	
	string thing = "good morning vietnam";
	
	wm.inputs.north = thing;
	cout << wm.inputs.north << endl;
	return 0;
	
}



int initialise_wall_node_positions(Wall_model &wm, const vector<double> &wall_xs, const vector<double> &wall_ys) {

	const double TOL = 1.0e-6;
	// Loop over x positions
	// Think about special values at edge.
	
	wm.inputs.M = wall_xs.size() + 2;
	wm.vars.x_vals.resize(wm.inputs.M*wm.inputs.N, 0.0);
	wm.vars.y_vals.resize(wm.inputs.M*wm.inputs.N, 0.0);
	wm.inputs.dx = wall_xs[1]-wall_xs[0];
	wm.inputs.dy = wm.inputs.D/(wm.inputs.N - 2.);
	
	
	
	double y_loc = wm.inputs.D + wm.inputs.y_h;
	double dx_test = 0;
	
	for (int i = 0; i < wm.inputs.M*wm.inputs.N; i+=wm.inputs.M) {
		int count = 0;
		for (int j = 1; j < wm.inputs.M-1; j++) {
			if ( count > 1 ) {
				dx_test = wall_xs[count] - wall_xs[count-1];
				if ( fabs(dx_test - wm.inputs.dx) >= TOL ) {
					cout << "For j = " << j << " Test = " << dx_test << ", actual = " << wm.inputs.dx << endl;
					cerr << "Clustering in x is not supported by the wall solver. ";
					cerr << "Please check that all cells in the fluid are equally spaced and sized. " ;
					cerr << endl;
					exit(MISMATCHED_DIMENSIONS);
				}
			}
			wm.vars.x_vals[i+j] = wall_xs[count];
			count++;			
		}
		for (int j = 0; j < wm.inputs.M; j++) {
			wm.vars.y_vals[i+j] = y_loc;
		}
		
		if ( i == 0 || i > wm.inputs.M*(wm.inputs.N-3)) {
			y_loc = y_loc - wm.inputs.dy/2.;			
		}
		else {
			y_loc = y_loc - wm.inputs.dy;
		}
	}
	
	for (int i = 0; i< wm.inputs.M*wm.inputs.N; i++) {
		int cols = (int) wm.inputs.M;
		if ( i%cols == 0 ) {
			wm.vars.x_vals[i] = wm.vars.x_vals[i+1] - wm.inputs.dx/2.;
		}
		else if ( i%cols == wm.inputs.M-1 ) {
			wm.vars.x_vals[i] = wm.vars.x_vals[i-1] + wm.inputs.dx/2.;
		}
	}
	
	double flow_dy = wall_ys[1] - wall_ys[0];
	double flow_dy_test = 0;
	
	for ( size_t i = 2; i < wall_ys.size() ; i++ ) {
		flow_dy_test = wall_ys[i] - wall_ys[i-1];
		if (fabs(flow_dy_test - flow_dy) >= TOL) {
			cerr << "Clustering in y is not supported by the wall solver. ";
			cerr << "Please check that all cells in the fluid are equally spaced and sized. " ;
			exit(MISMATCHED_DIMENSIONS);
		}
	}
	
	wall_printer(wm.vars.x_vals, wm.inputs);
	
	init_with_eilmer3(wm.inputs, wm.vars);

return 0;

}



int update_temperatures_from_fluxes(Wall_model &wm, double dt, const vector<double> q_wall, vector<double> &T_wall) {

	wm.inputs.dt = dt;
	
	update_coeff_matrix(wm.vars, wm.inputs);
	flux_calc(wm.inputs, wm.vars, q_wall);
	increment_in_time(wm.vars, wm.inputs);
	
	int count = wm.inputs.M*(wm.inputs.N-1)+1;
	for ( size_t i = 0; i < T_wall.size(); ++i ) {
		T_wall[i] = wm.vars.temps[count];
		count++;
	} 
	
return 0;


}


