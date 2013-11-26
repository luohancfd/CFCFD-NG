#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/util/source/useful.h"

#include "init.hh"
#include "adv.hh"
#include "e3con.hh"
#include "conjugateHeatUtil.hh"



int read_lookup_flux(list_of_inputs &inputs, list_of_vars &vars) {

	/*
	--------------------------------------------------------------------------
	--------------------------------------------------------------------------
	This function is written to be used in "flux_calc".  It reads data for 
	flux/temp used as part of the udf_lookup boundary condition.  

	It is used later on to determine the flux out of a particular node given
	its temperature.  This is a radiative flux.

	Data is stored in vector<double>s called "udf_temp_interp" and 
	"udf_flux_interp", which are both part of the vars structure.

	---------------------------------------------------------------------------
	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------------
	---------------------------------------------------------------------------
	*/


	const char* filename = "../user/flux-temp-lookup.txt";

	int count = 0;
	unsigned int numlines = 0;
	string buffer;
	ifstream tf(filename);

	while(getline(tf, buffer)) {	
		numlines++;
	}

	vars.udf_temp_interp.resize(numlines, 0.0);
	vars.udf_flux_interp.resize(numlines, 0.0);

	ifstream tft;
	tft.open(filename);

	bool filename_exists = file_exists(filename);
	if (filename_exists) {

		while(!tft.eof()) {
			
			tft >> vars.udf_temp_interp[count] >> vars.udf_flux_interp[count];	
			count++;
		}

	tft.close();
	}	

	else {
		cout << "No lookup table exists.  Please place the required file in location " << filename << endl;
		exit(MISMATCHED_DIMENSIONS);
	}	
	return SUCCESS;
}



int flux_calc( list_of_inputs &inputs, list_of_vars &vars, vector<double> q_wall) {

	/*
	--------------------------------------------------------------------------
	--------------------------------------------------------------------------
	This is a function to calculate the flux for any boundary of the wall.

	It uses values from the "inputs" structure (dims, timestep etc) and writes 
	fluxes to the "fluxes" vector which is part of the "vars" structure (both 
	of these are first used in "init.cxx" and are declared in "init.hh").

	The way it is set up at the moment (08/10/2013), the flow fluxes are brought
	in using the flow_fluxes structure, which is declared as part of "e3con.hh".

	
	---------------------------------------------------------------------------

	What does this function have to do:

		- Look at each of the north, south, east and west boundaries and
		  determine which boundary condition is applied.

		- Calculate fluxes for that boundary

		- Store the fluxes for each boundary in one overall vector vars.fluxes

	---------------------------------------------------------------------------
	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------------
	---------------------------------------------------------------------------
	*/

	vars.fluxes.resize(inputs.M*inputs.N, 0.0);

	// Iterator variable initialize
	int counter;

	// Create arrays for the fluxes through each of the four boundaries

	vector<double> nflux, sflux, eflux, wflux;
	nflux.resize(inputs.M, -1.0);
	sflux.resize(inputs.M, -1.0);
	eflux.resize(inputs.N, -1.0);
	wflux.resize(inputs.N, -1.0);


	// North wall

	if (inputs.north == "adiabatic") {
		for (int i = 0; i<inputs.M; i++){
			nflux[i] = 0.0;
			}
	}

	else if (inputs.north == "conv_rad") {
		counter = 0;
		for (int i = 0; i<inputs.M; i++) {
			nflux[counter] = inputs.he_n * (inputs.Te_n - vars.temps[i]) + inputs.emis_n*inputs.stefanBoltzmann*(pow(inputs.Te_n, 4) - pow(vars.temps[i], 4));
			counter++;
		}
	}

	else if (inputs.north == "flow") {
		for (int i = 1; i<inputs.M-1; i++) {
			nflux[i] = q_wall[i-1];
		}
		nflux[0] = q_wall[0];
		nflux[nflux.size()] = q_wall[q_wall.size()];
	}

	else if (inputs.north == "udf_flux") {
		for (int i = 0; i< inputs.M; i++) {
			nflux[i] = 300000.0;
		}
	}

	else if (inputs.north == "wall_fixedTBC") {

		double yf1, yf2;

		if (inputs.axi == 1) {
			yf1 = (0.5*(vars.y_vals[0] + vars.y_vals[inputs.M]))/vars.y_vals[0];
			yf2 = vars.y_vals[inputs.M]/vars.y_vals[0];
		}
		else {
			yf1 = 1;
			yf2 = 1;
		}

		for (int i = 0; i< inputs.M; i++) {
			nflux[i] = yf1*2.*inputs.k_22 * (vars.temps[i] - vars.temps[i+inputs.M])/inputs.dy; 			
		}
		for (int i = 1; i< inputs.M-1; i++) {
			nflux[i]+= yf2*(0.25*inputs.dy*inputs.k_11/(inputs.dx*inputs.dx))*(2.0*vars.temps[i] - vars.temps[i+1] - vars.temps[i-1]);
		}
	}

	else if (inputs.north == "udf_lookup") {
		for (int i = 0; i< inputs.M; i++) {
			nflux[i] = -1*fabs(linear_interp(vars.udf_temp_interp, vars.udf_flux_interp, vars.temps[i]));
		}
	}
		




	// South wall

	if (inputs.south == "adiabatic") {
		for (int i = 0; i<inputs.M; i++){
			sflux[i] = 0.0;
			}
	}
	
	else if (inputs.south == "conv_rad") {
		counter = 0;
		for (int i = inputs.M*(inputs.N-1); i<inputs.M*inputs.N; i++) {
			sflux[counter] = inputs.he_s * (inputs.Te_s - vars.temps[i]) + inputs.emis_s*inputs.stefanBoltzmann*(pow(inputs.Te_s, 4) - pow(vars.temps[i], 4));
			counter++;
		}
	}

	else if (inputs.south == "flow") {
		for (int i = 1; i<inputs.M-1; i++) {
			sflux[i] = q_wall[i-1];
		}
		sflux[0] = q_wall[0];
		sflux[sflux.size()-1] = q_wall[q_wall.size()-1];
		//		print_vector(sflux);
	}

	else if (inputs.south == "udf_flux") {
		for (int i = 0; i< inputs.M; i++) {
			sflux[i] = 100000.0;
		}
	}

	else if (inputs.south == "wall_fixedTBC") {
		int new_counter = 0;
		double yf1, yf2;
		if (inputs.axi == 1) {
			yf1 = (0.5*(vars.y_vals[inputs.M*inputs.N-1] + vars.y_vals[inputs.M*inputs.N-inputs.M-1]))/vars.y_vals[inputs.M*inputs.N-1];
			yf2 = vars.y_vals[inputs.M*inputs.N-inputs.M-1]/vars.y_vals[inputs.M*inputs.N-1];
		}
		else {
			yf1 = 1;
			yf2 = 1;
		}
		for (int i = inputs.M*(inputs.N-1); i< inputs.M*inputs.N; i++) {
			sflux[new_counter] = yf2*2.*inputs.k_22 * (vars.temps[i] - vars.temps[i-inputs.M])/inputs.dy;
			if (i != inputs.M*(inputs.N-1) && i != inputs.M*inputs.N-1) {
				sflux[new_counter] += 0.25*yf1*inputs.k_11*(inputs.dy/(inputs.dx*inputs.dx))*(2*vars.temps[i]-vars.temps[i-1]-vars.temps[i+1]);
			}
			new_counter++;
		}
	}

	else if (inputs.south == "udf_lookup") {
		int new_counter = 0;

		for (int i = inputs.M*(inputs.N-1); i< inputs.M*inputs.N; i++) {
			sflux[new_counter] = -1*fabs(linear_interp(vars.udf_temp_interp, vars.udf_flux_interp, vars.temps[i]));
			new_counter++;
		}
	}

	// East wall

	if (inputs.east == "adiabatic") {
		for (int i = 0; i<inputs.N; i++){
			eflux[i] = 0.0;
			}
	}

	else if (inputs.east == "conv_rad") {
		counter = 0;
		for (int i = inputs.M-1.; i<inputs.M*inputs.N; i+=inputs.M) {
			eflux[counter] = inputs.he_e * (inputs.Te_e - vars.temps[i]) + inputs.emis_e*inputs.stefanBoltzmann*(pow(inputs.Te_e, 4) - pow(vars.temps[i], 4));
			counter++;
		}
	}

	else if (inputs.east == "flow") {
		for (int i = 1; i<inputs.N-1; i++) {
			eflux[i] = q_wall[i-1];
		}
		eflux[0] = q_wall[0];
		eflux[eflux.size()] = q_wall[q_wall.size()];
	}

	else if (inputs.east == "udf_flux") {
		for (int i = 0; i< inputs.N; i++) {
			eflux[i] = 300000.0;
		}
	}

	else if (inputs.east == "wall_fixedTBC") {
		int new_counter = 0;
		double yf1, yf2;

		for (int i = inputs.M-1.; i<inputs.M*inputs.N; i+=inputs.M) {

			if (inputs.axi == 1) {
				yf1 = (0.5*(vars.y_vals[i-inputs.M]+vars.y_vals[i]))/vars.y_vals[i];
				yf2 = (0.5*(vars.y_vals[i+inputs.M]+vars.y_vals[i]))/vars.y_vals[i];
			}
			else {
				yf1 = 1;
				yf2 = 1;
			}
			eflux[new_counter] = 2.*inputs.k_11 * (vars.temps[i] - vars.temps[i-1])/inputs.dx;
			if (i != inputs.M-1 && i != inputs.M*inputs.N-1) {
				eflux[new_counter] += 0.25*inputs.k_22*(inputs.dx/(inputs.dy*inputs.dy))*(2*vars.temps[i]-yf1*vars.temps[i-inputs.M]-yf2*vars.temps[i+inputs.M]);
			}
			new_counter++;
		}
	}


	else if (inputs.east == "udf_lookup") {
		int new_counter = 0;
		for (int i = inputs.M-1.; i<inputs.M*inputs.N; i+=inputs.M) {
			eflux[new_counter] = -1*fabs(linear_interp(vars.udf_temp_interp, vars.udf_flux_interp, vars.temps[i]));			
			new_counter++;
		}
	}


	// West wall

	if (inputs.west == "adiabatic") {
		for (int i = 0; i<inputs.N; i++){
			wflux[i] = 0.0;
			}
	}

	else if (inputs.west == "conv_rad") {
		counter = 0;
		for (int i = 0.; i<inputs.M*inputs.N; i+=inputs.M) {
			wflux[counter] = inputs.he_w * (inputs.Te_w - vars.temps[i]) + inputs.emis_w*inputs.stefanBoltzmann*(pow(inputs.Te_w, 4) - pow(vars.temps[i], 4));
			counter++;
		}
	}

	else if (inputs.west == "flow") {
		for (int i = 1; i<inputs.N-1; i++) {
			wflux[i] = q_wall[i-1];
		}
		wflux[0] = q_wall[0];
		wflux[wflux.size()] = q_wall[q_wall.size()];
	}

	else if (inputs.west == "udf_flux") {
		for (int i = 0; i< inputs.N; i++) {
			wflux[i] = 300000.0;
		}
	}

	else if (inputs.west == "wall_fixedTBC") {
		int new_counter = 0;
		double yf1, yf2;

		for (int i = 0.; i<inputs.M*inputs.N; i+=inputs.M) {

			if (inputs.axi == 1) {
				yf1 = (0.5*(vars.y_vals[i-inputs.M]+vars.y_vals[i]))/vars.y_vals[i];
				yf2 = (0.5*(vars.y_vals[i+inputs.M]+vars.y_vals[i]))/vars.y_vals[i];
			}
			else {
				yf1 = 1;
				yf2 = 1;
			}
			wflux[new_counter] = 2.*inputs.k_11 * (vars.temps[i] - vars.temps[i+1])/inputs.dx;
			if (i != 0 && i != inputs.M*(inputs.N-1)) {
				wflux[new_counter] += 0.25*inputs.k_22*(inputs.dx/(inputs.dy*inputs.dy))*(2*vars.temps[i]-yf1*vars.temps[i-inputs.M]-yf2*vars.temps[i+inputs.M]);
			}
			new_counter++;
		}
	}

	else if (inputs.west == "udf_lookup") {
		int new_counter = 0;
		for (int i = 0.; i<inputs.M*inputs.N; i+=inputs.M) {
			wflux[new_counter] = -1*fabs(linear_interp(vars.udf_temp_interp, vars.udf_flux_interp, vars.temps[i]));
			new_counter++;
		}
	}
	
	

	
	
	
	
	// Take the individual vectors and put them into vars.fluxes
	// Note that the corner nodes have no flux through them intentionally
	// so that the flux at any node is only one dimensional.

	// They should be sufficiently far from the area of interest that this
	// won't affect the solution.

	double yfp;
	double yfm;


	// North face
	counter = 0;

	for (int i = 0; i<inputs.M; i++) {
		if (inputs.axi == 0) {
			yfp = 1.0;
			yfm = 1.0;
		}
		else {
			yfp = 0.0;
			yfm = (vars.y_vals[i]-inputs.dy/4.)/(vars.y_vals[i]-inputs.dy/8.);
		}
		vars.fluxes[i]+=nflux[counter]*4.*inputs.dt*yfm/(inputs.rho_w*inputs.c_w*inputs.dy);
		counter++;
	}

	// South face
	counter = 0;
	for (int i = inputs.M*(inputs.N-1); i<inputs.M*inputs.N; i++) {
		if (inputs.axi == 0) {
			yfp = 1.0;
			yfm = 1.0;
		}
		else {
			yfp = (vars.y_vals[i]+inputs.dy/4.)/(vars.y_vals[i]+inputs.dy/8.);
			yfm = 0.0;
		}
		vars.fluxes[i]+=sflux[counter]*4.*inputs.dt*yfp/(inputs.rho_w*inputs.c_w*inputs.dy);
		counter++;
	}

	// East face
	counter = 0;
	for (int i = inputs.M-1; i<inputs.M*inputs.N; i+=inputs.M) {
		vars.fluxes[i]+=eflux[counter]*4.*inputs.dt/(inputs.rho_w * inputs.c_w * inputs.dx);
		counter++;
	}

	// West face
	counter = 0;
	for (int i = 0; i<inputs.M*inputs.N; i+=inputs.M) {
		vars.fluxes[i]+=wflux[counter]*4.*inputs.dt/(inputs.rho_w*inputs.c_w*inputs.dx);
		counter++;
	}

	return SUCCESS;
}


int update_coeff_matrix(list_of_vars &vars, list_of_inputs &inputs) {


	/*
	--------------------------------------------------------------------------
	--------------------------------------------------------------------------
	This function updates the coefficient matrix with the timestep.

	
	---------------------------------------------------------------------------

	17/10/2013- The vector of variables passed between this function and Eilmer3
				hasn't been written yet so I haven't passed the proper reference
				to this function for the location of dt.  For the moment, the 
				constant value of inputs.dt is used.

	
	20/10/2013- I think I've decided to add dt to inputs.dt because otherwise
				there are probably too many to change.  Don't know about this
				though.
	---------------------------------------------------------------------------
	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------------
	---------------------------------------------------------------------------
	*/

	for (int i = 0; i< inputs.M*inputs.N; i++) {
		for (int j = 0; j< inputs.M*inputs.N; j++ ) {
			vars.A.set(i, j, vars.B.get(i,j)*inputs.dt + vars.C.get(i,j));			
		}	
	}

	return SUCCESS;
}

int increment_in_time(list_of_vars &vars, list_of_inputs &inputs) {

/*
--------------------------------------------------------------------------
--------------------------------------------------------------------------
Runs the solution T[p+1] = [A]*T[p] + [Q].  This function can be called by
external solvers to solve the wall.

---------------------------------------------------------------------------
\author Jared Clifford
\version October 2013

---------------------------------------------------------------------------
---------------------------------------------------------------------------
*/

	vector<double> intermediate_temps; // Intermediate storage location

	intermediate_temps.resize(vars.temps.size(), 0.0);
        vector_mul(vars.A, vars.temps, intermediate_temps); // [A]*T[p]

	add_vectors(vars.temps, intermediate_temps, vars.fluxes); // [A]*T[p] + [Q]

	return 0;
}
	











