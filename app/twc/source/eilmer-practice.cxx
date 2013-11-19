#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include "../util/config_parser.hh"
#include "../math/no_fuss_linear_algebra_jc.hh"
#include "../util/useful.h"



#include "init.hh"
#include "adv.hh"
#include "conjugateHeatUtil.hh"
#include <stdio.h>
#include "e3con.hh"

using namespace std;


int main(){
	vector<double> q_wall; q_wall.resize(8, 100000.);
	vector<double> T_wall; T_wall.resize(8, 300.);
	
	double dublist[] = {0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008};
	const vector<double> wall_xs(dublist, dublist + sizeof(dublist) / sizeof(double)); 
	const vector<double> wall_ys(dublist, dublist + sizeof(dublist) / sizeof(double)); 
	
	double dt = 0.001;
	double dt_plot = 0.1;
	double time = 1;
	double telapsed = 0;
	double time_to_print;
	int print_number = 1;
	string fname = "../user/config.cfg";
	
	Wall_model* wm = initialise_wall_model(fname, dt_plot);
	initialise_wall_node_positions(*wm, wall_xs, wall_ys);
	
	old_file_delete_catcher((*wm).inputs);
	
	time_to_print = (*wm).inputs.dt_plot;
	
	
	while ( telapsed < time ) {
		update_temperatures_from_fluxes(*wm, dt, q_wall, T_wall);
		telapsed+=dt;
		
		cout << "New step... " << endl;
		// wall_printer((*wm).vars.temps, (*wm).inputs);
		for (int i = 0; i< T_wall.size(); i++) {
			cout << setw(7) << setprecision(5)<< T_wall[i];
		}
		cout << endl;
		
		time_to_print = time_to_print - (*wm).inputs.dt;
		cout << "dt_plot = " << (*wm).inputs.dt_plot << " dt = " << (*wm).inputs.dt << " Time to plot = " << time_to_print << endl;
		double tol = 1.0e-9;
		if (time_to_print < tol) {
			write_soln((*wm).vars, (*wm).inputs, telapsed, print_number);
			time_to_print = (*wm).inputs.dt_plot;
			print_number++;
		}
		
	}
	

}
