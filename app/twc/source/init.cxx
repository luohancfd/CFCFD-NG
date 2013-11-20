#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

#include "../../../lib/util/source/config_parser.hh"
#include "../../../lib/nm/source/no_fuss_linear_algebra.hh"
#include "../../../lib/util/source/useful.h"

#include "init.hh"
#include "e3con.hh"
#include "conjugateHeatUtil.hh"

using namespace std;

void grab_config(string fname, list_of_inputs &inputs) {

    /* 
       
       ---------------------------------------------------------------------
       ---------------------------------------------------------------------
       This function grabs data from the user-defined configuration 
    file called config.cfg.  It uses the configuration file parser
    which comes with Eilmer3.
    
    The original files config_parser.cxx and config_parser.hh are
    located in "cfcfd3/lib/util/source".  They have been copied and
    moved to one directory above this one ie "../util".

	There are now a bunch of error catching if loops in here (added 15/10/2013)

	---------------------------------------------------------------------

	Updated to include a filename (string fname) as an argument
	which tells the solver the location of the configuration file.
	(18/11/2013)	

	---------------------------------------------------------------------

    The original parser: Author Rowan Gollan, 2006.

    This function: Author Jared Clifford, 2013. 

	---------------------------------------------------------------------
	---------------------------------------------------------------------
	
*/

    bool truth_val;
    
    //cout << endl << endl <<"--------------------------------------------------------------------" << endl;
    //cout << "\t \t \t \t  ERRORS \t \t \t \t " << endl;
    //cout << "--------------------------------------------------------------------" << endl;
    
    // Open the configuration file
    ConfigParser cfg( fname );
    
    // Bring in all the data from the config.cfg file.
    truth_val = cfg.parse_double( string("Constants"), string("stefanBoltzmann"),
				  inputs.stefanBoltzmann, -1.0 );
    if (inputs.stefanBoltzmann == -1.0) {
	// // // cerr << "No input provided for Stefan Boltzmann Constant. Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }

    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Geometry"), string("D"),
				  inputs.D, -1.0 );
    if (inputs.D < 0.0) {
	cerr << "No input provided for plate width in NS direction D. Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }

    // -------------------------------------------------------------------------

    truth_val = cfg.parse_double( string("Geometry"), string("L"),
				  inputs.L, -1.0 );
    if (inputs.L < 0.0) {
	cerr << "No input provided for plate width in EW direction flow L. Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }

    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Geometry"), string("M"),
				  inputs.M, -1.0 );
    if (inputs.M < 5.0) {
	cerr << "Number of nodes in EW direction is too low/not provided. Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Geometry"), string("N"),
				  inputs.N, -1.0 );
    if (inputs.N < 5.0) {
	cerr << "Number of nodes in NS direction is too low/not provided. Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }

    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Geometry"), string("y_h"),
				  inputs.y_h, -1.0 ); // This isn't important if not provided since
    // axisym_checker fills it in anyway.
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_int( string("Geometry"), string("axi"),
			       inputs.axi, -1.0 );
    if (inputs.axi == -1.0) {
	cerr << "Axisymmetric Flag is not specified. Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Wall-Material-Properties"), string("k_11"),
				  inputs.k_11, -1.0 );
    if (inputs.k_11 < 0.0) {
	cerr << "No value/negative value for thermal conductivity k_11 provided.  Default value of 0.0 is applied." << endl;
	inputs.k_11 = 0.0;
	// exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Wall-Material-Properties"), string("k_12"),
				  inputs.k_12, -1.0 );
    if (inputs.k_12 < 0.0) {
	cerr << "No value/negative value for thermal conductivity k_12 provided.  Default value of 0.0 is applied." << endl;
	inputs.k_12 = 0.0;
	// exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Wall-Material-Properties"), string("k_21"),
				  inputs.k_21, -1.0 );
    if (inputs.k_21 < 0.0) {
	cerr << "No value/negative value for thermal conductivity k_21 provided.  Default value of 0.0 is applied." << endl;
	inputs.k_21 = 0.0;
	// exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Wall-Material-Properties"), string("k_22"),
				  inputs.k_22, -1.0 );
    if (inputs.k_22 < 0.0) {
	cerr << "No value/negative value for thermal conductivity k_22 provided.  Default value of 0.0 is applied." << endl;
	inputs.k_22 = 0.0;
	// exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Wall-Material-Properties"), string("rho_w"),
				  inputs.rho_w, -1.0 );
    if (inputs.rho_w < 0.0) {
	cerr << "Non physical value/No value for wall density provided.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------

    truth_val = cfg.parse_double( string("Wall-Material-Properties"), string("c_w"),
				  inputs.c_w, -1.0 );
    if (inputs.c_w < 0.0) {
	cerr << "Non physical value/No value for wall thermal conductivity provided.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }

    // -------------------------------------------------------------------------
    
    
    truth_val = cfg.parse_double( string("Timestepping"), string("dt"),
				  inputs.dt, -1.0 );
    if (inputs.dt < 0.0) {
	cerr << "Non physical value/No value for timestep provided.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Timestepping"), string("dt_plot"),
					inputs.dt_plot, -1.0 );
    if (inputs.dt < 0.0) {
	cerr << "Non physical value/No value for plotting timestep provided.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Timestepping"), string("t_max"),
				  inputs.t_max, -1.0 );
    if (inputs.dt < 0.0) {
	cerr << "Non physical value/No value for maximum time provided.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_string( string("Boundary-Conditions"), string("north"),
				  inputs.north, "Brian" );
    
    if (inputs.north == "Brian") {
		cerr << "No input for North BC type.  Exiting..." << endl;
		exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_string( string("Boundary-Conditions"), string("east"),
				  inputs.east, "Brian" );
    if (inputs.east == "Brian") {
	cerr << "No input for East BC type.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }

    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_string( string("Boundary-Conditions"), string("south"),
				  inputs.south, "Brian" );
    if (inputs.south == "Brian") {
	cerr << "No input for North BC type.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_string( string("Boundary-Conditions"), string("west"),
				  inputs.west, "Brian" );
    if (inputs.west == "Brian") {
	cerr << "No input for West BC type.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("Te_n"),
				  inputs.Te_n, -1.0 );
    if (inputs.north == "conv_rad" && inputs.Te_n < 0.0) {
	cerr << "North boundary specified as conv_rad but no positive environment temp given.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.north != "conv_rad" && inputs.Te_n > 0.0) {
	//cerr << "Unnecessary value for north boundary environment temp specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("Te_s"),
				  inputs.Te_s, -1.0 );
    if (inputs.south == "conv_rad" && inputs.Te_s < 0.0) {
	cerr << "South boundary specified as conv_rad but no positive environment temp given.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.south != "conv_rad" && inputs.Te_s > 0.0) {
	//cerr << "Unnecessary value for south boundary environment temp specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("Te_e"),
				  inputs.Te_e, -1.0 );
    if (inputs.east == "conv_rad" && inputs.Te_e < 0.0) {
	cerr << "East boundary specified as conv_rad but no positive environment temp given.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.east != "conv_rad" && inputs.Te_e > 0.0) {
	//cerr << "Unnecessary value for east boundary environment temp specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("Te_w"),
				  inputs.Te_w, -1.0 );
    if (inputs.west == "conv_rad" && inputs.Te_w < 0.0) {
	cerr << "West boundary specified as conv_rad but no positive environment temp given.  Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.west != "conv_rad" && inputs.Te_w > 0.0) {
	//cerr << "Unnecessary value for west boundary environment temp specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("he_n"),
				  inputs.he_n, -1.0 );
    if (inputs.north == "conv_rad" && inputs.he_n < 0.0) {
	cerr << "North boundary specified as conv_rad but no physical environment convection coefficient given.  Convection ignored from this BC." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.north != "conv_rad" && inputs.he_n > 0.0) {
	//cerr << "Unnecessary value for north boundary environment convection coefficient specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("he_s"),
				  inputs.he_s, -1.0 );
    if (inputs.south == "conv_rad" && inputs.he_s < 0.0) {
	cerr << "South boundary specified as conv_rad but no physical environment convection coefficient given.  Convection ignored from this BC." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.south != "conv_rad" && inputs.he_s > 0.0) {
	//cerr << "Unnecessary value for south boundary environment convection coefficient specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------

    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("he_e"),
				  inputs.he_e, -1.0 );
    if (inputs.east == "conv_rad" && inputs.he_e < 0.0) {
	cerr << "East boundary specified as conv_rad but no physical environment convection coefficient given.  Convection ignored from this BC." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.east != "conv_rad" && inputs.he_e > 0.0) {
	//cerr << "Unnecessary value for east boundary environment convection coefficient specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("he_w"),
				  inputs.he_w, -1.0 );
    if (inputs.west == "conv_rad" && inputs.he_w < 0.0) {
	cerr << "West boundary specified as conv_rad but no physical environment convection coefficient given.  Convection ignored from this BC." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    else if (inputs.west != "conv_rad" && inputs.he_w > 0.0) {
	//cerr << "Unnecessary value for west boundary environment convection coefficient specified.  Value is ignored in calcultions." << endl;
    }
    
    // -------------------------------------------------------------------------
    
    
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("emis_n"),
				  inputs.emis_n, -1.0 );
    if (inputs.north == "conv_rad" && inputs.emis_n <= 0.0) {
	cerr << "North boundary specified as conv_rad but no physical surface emissivity coefficient given.  Radiation ignored from this BC." << endl;
    }
    else if (inputs.north != "conv_rad" && inputs.emis_n > 0.0) {
	//cerr << "Unnecessary value for north boundary environment surface emissivity specified.  Value is ignored in calcultions." << endl;
    }
    
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("emis_s"),
				  inputs.emis_s, -1.0 );
    if (inputs.south == "conv_rad" && inputs.emis_s <= 0.0) {
	cerr << "South boundary specified as conv_rad but no physical surface emissivity coefficient given.  Radiation ignored from this BC." << endl;
    }
    else if (inputs.south != "conv_rad" && inputs.emis_s > 0.0) {
	//cerr << "Unnecessary value for south boundary environment surface emissivity specified.  Value is ignored in calcultions." << endl;
    }
    
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("emis_e"),
				  inputs.emis_e, -1.0 );
    if (inputs.east == "conv_rad" && inputs.emis_e <= 0.0) {
	cerr << "East boundary specified as conv_rad but no physical surface emissivity coefficient given.  Radiation ignored from this BC." << endl;
    }
    else if (inputs.east != "conv_rad" && inputs.emis_e > 0.0) {
	//cerr << "Unnecessary value for east boundary environment surface emissivity specified.  Value is ignored in calcultions." << endl;
    }
    
    
    // -------------------------------------------------------------------------

    truth_val = cfg.parse_double( string("Boundary-Conditions"), string("emis_w"),
					inputs.emis_w, -1.0 );
    if (inputs.west == "conv_rad" && inputs.emis_w <= 0.0) {
	cerr << "West boundary specified as conv_rad but no physical surface emissivity coefficient given.  Radiation ignored from this BC." << endl;
    }
    else if (inputs.west != "conv_rad" && inputs.emis_w > 0.0) {
	//cerr << "Unnecessary value for west boundary environment surface emissivity specified.  Value is ignored in calcultions." << endl;
    }
    
    
    // -------------------------------------------------------------------------
    
    if (inputs.north == "conv_rad" && inputs.he_n <= 0.0 && inputs.emis_n <= 0.0) {
	cerr << "Insufficient information provided to compute north conv_rad boundary" << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    if (inputs.south == "conv_rad" && inputs.he_s <= 0.0 && inputs.emis_s <= 0.0) {
	cerr << "Insufficient information provided to compute south conv_rad boundary" << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    if (inputs.east == "conv_rad" && inputs.he_e <= 0.0 && inputs.emis_e <= 0.0) {
	cerr << "Insufficient information provided to compute east conv_rad boundary" << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    if (inputs.west == "conv_rad" && inputs.he_w <= 0.0 && inputs.emis_w <= 0.0) {
	cerr << "Insufficient information provided to compute west conv_rad boundary" << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_int( string("History"), string("flag_hist"),
			       inputs.flag_hist, 0 );
    // A value of 1 means that history is recorded.  Default value of 0 means
    // that no history is recorded.
    
    
    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_int( string("History"), string("nn_hist"),
			       inputs.nn_hist, -1.0 );
    if (inputs.flag_hist == 0) {
	if (inputs.nn_hist < 0 || inputs.nn_hist > inputs.M*inputs.N-1) {
	    cerr << "No node specified for recording history. Exiting..." << endl;
	    exit(MISMATCHED_DIMENSIONS);
	    
	}	
    }

    // -------------------------------------------------------------------------
    
    truth_val = cfg.parse_double( string("History"), string("dt_hist"),
				  inputs.dt_hist, -1.0 );
    if (inputs.dt_hist < 0 || inputs.dt_hist > inputs.t_max ) {
	cerr << "Inappropriate timestep specified for recording history. Exiting..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
        
    truth_val = cfg.parse_double( string("Initial-Conditions"), string("T_init"),
				  inputs.T_init, -1.0 );
    if (inputs.T_init < 0) {
	cerr << "No initial temperature set for the wall.  Bailing out..." << endl;
	exit(MISMATCHED_DIMENSIONS);
    }
    
    // -------------------------------------------------------------------------
    
    
    //cout << "--------------------------------------------------------------------" << endl;
    //cout << "--------------------------------------------------------------------" << endl << endl << endl;
    // Empty test to stop compiler warning
    if (truth_val)
	;
 
}

int init_with_eilmer3(list_of_inputs &inputs, list_of_vars &vars){
    ar_checker(inputs);
    axisym_checker(inputs);
    init_tempprof(inputs, vars);
    vars.fluxes.resize(inputs.M*inputs.N, 0.0);
    init_coeff_matrix(inputs, vars);
    return SUCCESS;
}

int init_tr( list_of_inputs &inputs, list_of_vars &vars) {
	
/* 

	---------------------------------------------------------------------
	---------------------------------------------------------------------
	This is the highest level function written into init.cxx.  It runs
	all the relevant initializing functions (also contained in this file)
	in the correct order.

	---------------------------------------------------------------------

    \author Jared Clifford
	\version October 2013
	---------------------------------------------------------------------
	---------------------------------------------------------------------

*/

	string fname = "../user/config.cfg";

	grab_config(fname, inputs); // Reads the config file input

	axisym_checker(inputs); // Check axisymmetry

	init_tempprof(inputs, vars); // Initialize the temps vector

	vars.fluxes.resize(inputs.M*inputs.N, 0.0); // Initialize the fluxes vector
	vars.x_vals.resize(inputs.M*inputs.N, 0.0); // Initiaize x_vals storage vector
	vars.y_vals.resize(inputs.M*inputs.N, 0.0); // Initiaize y_vals storage vector

	write_geom(inputs, vars); // Populates the x_vals and y_vals vectors
	
	ar_checker(inputs); // Checks aspect ratio

	init_coeff_matrix(inputs, vars); // Initializes the [B] and [C] parts of the coefficient matrix
	return SUCCESS;
}


int axisym_checker(list_of_inputs &inputs) {

/*

	---------------------------------------------------------------------
	---------------------------------------------------------------------
	This is a utility function that checks for planar/axisymmetric and 
	updates y_h accordingly.  It also exits the program if the axisymmetric
	flag is on but no y_h value is specified.

	---------------------------------------------------------------------

	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------
	---------------------------------------------------------------------

*/

	if (inputs.axi == 0) {
		inputs.y_h = 0.;
		cout << "Planar solution will be computed." << endl;
	}
	else {
		cout << endl << "Axisymmetric solution will be computed." << endl << endl;
		if (inputs.y_h <= 0.) {
			cerr << "^But, no axis of symmetry specified.  Program exiting..." << endl << endl;
			exit(MISMATCHED_DIMENSIONS);
		}
	}
	return SUCCESS;
}

int write_geom(list_of_inputs &inputs, list_of_vars &vars) {

/*

	---------------------------------------------------------------------
	---------------------------------------------------------------------
	This function write x and y values for each node.  It updates
	the vectors x_vals and y_vals which are part of the "vars" 
	structure.

	---------------------------------------------------------------------

	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------
	---------------------------------------------------------------------

*/
	
	inputs.dx = inputs.L/(inputs.M - 2.0); // uniform x spacing
	inputs.dy = inputs.D/(inputs.N - 2.0); // uniform y spacing

	// Note that the -2.0 accounts for the half cells on the borders.
	

	double y_loc = inputs.D+inputs.y_h;
	double x_loc;
	for (int i = 0; i < inputs.M*inputs.N; i+=inputs.M) {
		x_loc = 0;
		
		for (int j = 0; j < inputs.M; j++) {
			vars.x_vals[i+j] = x_loc;
			if (i>inputs.M*inputs.N-(inputs.M+1.)) {
				vars.y_vals[i+j] = inputs.y_h;
			}
			else {
				vars.y_vals[i+j] = y_loc + inputs.y_h;
			}
			if (j == 0  || j == inputs.M-2.0) { 
				x_loc = x_loc + inputs.dx/2.0;
			}
			else {
				x_loc = x_loc + inputs.dx;
			}
		}
		if ( i == 0 || i > (inputs.M*inputs.N - (inputs.M * 3.0))) {
			y_loc = y_loc - inputs.dy/2.0;
		}
		else {
			y_loc = y_loc - inputs.dy;
		}
			
		
	}

	return SUCCESS;

}




int init_tempprof(list_of_inputs &inputs, list_of_vars &vars) {

/*

  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  Initializes the temperature field.  This can be read from a text
  file or changed using the for and if statements as part of this 
  function.
  ---------------------------------------------------------------------
  
  \author Jared Clifford
  \version October 2013
  
  ---------------------------------------------------------------------
  ---------------------------------------------------------------------
  
*/
    
    vars.temps.resize(inputs.M*inputs.N, inputs.T_init); // Initializes the temp field.
    
    const char* filename = "../user/init_temp.txt"; //Filename of initial temp field
    bool temp_text_exists = file_exists(filename); // Evaluates whether file exists
    
    if (temp_text_exists == 1) {
	cout << "Temp profile needs to be read-- function needs to be written " << endl;
    }
/*
	else {
		for (int i = 0; i<inputs.M*inputs.N; i++) {
			int num = (int) inputs.M;
			if ( i%num == 0 ) {
				vars.temps[i] = 368.0;
			}
			else {
				vars.temps[i] = 293.0;
			}
			//if (i<inputs.M) {
			//	vars.temps[i] = 368;
			//}
			//else {
			//	vars.temps[i] = 293;
			//}
		}
	}
*/
	
    return SUCCESS;

}

int ar_checker(list_of_inputs &inputs) {

/*

	---------------------------------------------------------------------
	---------------------------------------------------------------------
	This function checks the aspect ratio of the cells and bails out if 
	AR > 5.

	---------------------------------------------------------------------

	\author Jared Clifford
	\version November 2013

	---------------------------------------------------------------------
	---------------------------------------------------------------------

*/

	double ar_limit = 5.0;

	if (inputs.dx/inputs.dy < 1/ar_limit || inputs.dx/inputs.dy >ar_limit) {
		cout << "Aspect ratio = " << inputs.dx/inputs.dy << endl;
		cerr << "Aspect Ratio is greater than 5.  Revise cell spacing.  Bailing out...";
		exit(MISMATCHED_DIMENSIONS);
	}

	return SUCCESS;
}

int init_coeff_matrix(list_of_inputs &inputs, list_of_vars &vars) {

/* 	
	---------------------------------------------------------------------
	---------------------------------------------------------------------

	This function is designed to form the coefficient matrix from the input
	parameters contained in the "inputs" structure.

	Because of the discretization required to allow each boundary to be an 
	interface with the flow, there are a lot of special case nodes which
	have to specified in this formulation.

	Special formulations are all hard coded.

  
	Note that the matrix does not ever need to be inverted since an explicit 
	solution has been developed. 

	
	Note that the matrix is formulated using [A] = dt*[B] + [C].  This has
	been done for the purpose of quickly updating the coefficient matrix
	due to changing the timestep.  This prevents the full program having to 
	form each element at each timestep, since the matrices [B] and [C] are
	invariant.

	The accompanying program adv_coeff_matrix performs the [A] = dt*[B] + [C]
	for the local timestep dt.  This ONLY CHANGES THE MATRIX [A]- [B] AND [C] 
	stay constant once they are set using this function.

	---------------------------------------------------------------------

	\author Jared Clifford
	\version October 2013

	---------------------------------------------------------------------
	---------------------------------------------------------------------
*/

// Size matrices
vars.A.resize(inputs.M*inputs.N, inputs.M*inputs.N);
vars.B.resize(inputs.M*inputs.N, inputs.M*inputs.N);
vars.C.resize(inputs.M*inputs.N, inputs.M*inputs.N);

// These are the gamma parameters as defined for each node.
double gam11, gam12, gam21, gam22; 

// These stand for "y fraction plus" and "y fraction minus". They are used to 
// account for the axisymmetric case.  For the planar case they are always equal to 1.0.
double yfp, yfm; 


/*
	-------------------------------------------------------------------------
	-------------------------------------------------------------------------

	The initial part of this code writes the [B] matrix.

	-------------------------------------------------------------------------
	
	\author Jared Clifford
	\version October 2013

	-------------------------------------------------------------------------
	-------------------------------------------------------------------------

*/


// These are alpha parameters (thermal diffusivity, generally defined for the entire 
// formulation.
double alph11 = inputs.k_11/(inputs.rho_w*inputs.c_w);
double alph12 = inputs.k_12/(inputs.rho_w*inputs.c_w);
double alph21 = inputs.k_21/(inputs.rho_w*inputs.c_w);
double alph22 = inputs.k_22/(inputs.rho_w*inputs.c_w);



// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
// 					NORTH AND INNER NORTH BOUNDARIES				   //
// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //


// North-West Corner (Node 0)------------------------------------------

int node = 0;
if (inputs.axi == 0) {
	yfp = 1.0;
	yfm = 1.0;
}
else {
	yfp = 0.0;
	yfm = (vars.y_vals[0]-inputs.dy/4.)/(vars.y_vals[0]-inputs.dy/8.);
}

gam11 = 8*alph11/pow(inputs.dx,2);
gam22 = 8*alph22/pow(inputs.dy,2);
gam12 = 4*alph12/(inputs.dx*inputs.dy);
gam21 = 4*alph21/(inputs.dx*inputs.dy);

vars.B.set(0,0,(-gam22*yfm+gam21*yfm-gam11+gam12));
vars.B.set(0,1,(gam12+gam11-gam21*yfm));
vars.B.set(0,inputs.M,(gam22*yfm+gam21*yfm-gam12));
vars.B.set(0,inputs.M+1.,-1.*(gam12+gam21*yfm));



// North-East Corner --------------------------------------------------
node = inputs.M-1;

// yfp and yfm stay the same, so do the gamma values

vars.B.set(node,node,(-gam22*yfm-gam21*yfm-gam11-gam12));
vars.B.set(node,node-1.,(gam21*yfm+gam11-gam12));
vars.B.set(node,node+inputs.M,(gam22*yfm-gam21*yfm+gam12));
vars.B.set(node,node+inputs.M-1.,(gam12+gam21*yfm));



// North-North-West Corner (Node 1)------------------------------------------

node = 1.;
if (inputs.axi == 0) {
	yfp = 1.0;
	yfm = 1.0;
}
else {
	yfp = 0.0;
	yfm = (vars.y_vals[node]-inputs.dy/4.)/(vars.y_vals[node]-inputs.dy/8.);
}

gam11 = (4/3.)*alph11/pow(inputs.dx,2);
gam22 = 8.*alph22/pow(inputs.dy,2);
gam12 = (4/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (4/3.)*alph21/(inputs.dx*inputs.dy);

vars.B.set(node, node,(-gam22*yfm-3*gam11));
vars.B.set(node, node-1.,(gam21*yfm+2*gam11-gam12));
vars.B.set(node, node+1.,(gam11-gam21*yfm+gam12));
vars.B.set(node, node+inputs.M,(gam22*yfm));
vars.B.set(node, node+inputs.M+1.,-1.*(gam12+gam21*yfm));
vars.B.set(node, node+inputs.M-1.,(gam12+gam21*yfm));


// North-North-East Corner ---------------------------------------------------

node = inputs.M-2.;

vars.B.set(node, node,(-gam22*yfm-3*gam11));
vars.B.set(node, node-1.,(gam21*yfm+gam11-gam12));
vars.B.set(node, node+1.,(2*gam11-gam21*yfm+gam12));
vars.B.set(node, node+inputs.M,(gam22*yfm));
vars.B.set(node, node+inputs.M+1.,-1.*(gam12+gam21*yfm));
vars.B.set(node, node+inputs.M-1.,(gam12+gam21*yfm));


// The rest of the North Border -----------------------------------------------



gam11 = alph11/pow(inputs.dx,2);
gam22 = 8.*alph22/pow(inputs.dy,2);
gam12 = alph12/(inputs.dx*inputs.dy);
gam21 = alph21/(inputs.dx*inputs.dy);


for (int i = 2; i<inputs.M-2.; i++) {
	
	node = i;
	vars.B.set(node, node,(-gam22*yfm-2*gam11));
	vars.B.set(node, node-1.,(gam21*yfm-gam12+gam11));
	vars.B.set(node, node+1.,(gam11-gam21*yfm+gam12));
	vars.B.set(node, node+inputs.M,(gam22*yfm));
	vars.B.set(node, node+inputs.M+1.,-1.*(gam12+gam21*yfm));
	vars.B.set(node, node+inputs.M-1.,(gam12+gam21*yfm));

}


// North-West-West Corner -----------------------------------------------------

node = inputs.M;

if (inputs.axi == 0) {
	yfp = 1.0;
	yfm = 1.0;
}
else {
	yfp = (vars.y_vals[node]+inputs.dy/4.)/(vars.y_vals[node]-inputs.dy/8.);
	yfm = (vars.y_vals[node]-inputs.dy/2.)/(vars.y_vals[node]-inputs.dy/8.);
}

gam11 = 8.*alph11/pow(inputs.dx,2);
gam22 = (4/3.)*alph22/pow(inputs.dy,2);
gam12 = (4/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (4/3.)*alph21/(inputs.dx*inputs.dy);


vars.B.set(node, node,(-gam22*(2*yfp+yfm)-gam21*(yfp-yfm)-gam11));
vars.B.set(node, node+1.,(gam21*(yfp-yfm)+gam11));
vars.B.set(node, node+inputs.M,(gam22*yfm+gam21*yfm-gam12));
vars.B.set(node, node+inputs.M+1.,-1.*(gam12+gam21*yfm));
vars.B.set(node, node-inputs.M,(2*gam22*yfp-gam21*yfp+gam12));
vars.B.set(node, node-inputs.M+1.,(gam12+gam21*yfp));


// North-East-East Corner -----------------------------------------------------

node = 2*inputs.M -1.;

vars.B.set(node, node,(-gam22*(2*yfp+yfm)+gam21*(yfp-yfm)-gam11));
vars.B.set(node, node-1.,(gam21*(yfm-yfp)+gam11));
vars.B.set(node, node+inputs.M,(gam22*yfm-gam21*yfm+gam12));
vars.B.set(node, node+inputs.M-1.,(gam12+gam21*yfm));
vars.B.set(node, node-inputs.M,(2*gam22*yfp+gam21*yfp-gam12));
vars.B.set(node, node-inputs.M-1.,-1.*(gam12+gam21*yfp));


// Inner North-West Corner -----------------------------------------------------

node = inputs.M+1.;

gam11 = (4/3.)*alph11/pow(inputs.dx,2);
gam22 = (4/3.)*alph22/pow(inputs.dy,2);
gam12 = (4/9.)*alph12/(inputs.dx*inputs.dy);
gam21 = (4/9.)*alph21/(inputs.dx*inputs.dy);

vars.B.set(node, node, (-gam22*(2*yfp+yfm)-3*gam11));
vars.B.set(node, node+1, (gam21*(yfp-yfm)+gam11));
vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-2*gam11));
vars.B.set(node, node+inputs.M, (gam22*yfm));
vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
vars.B.set(node, node-inputs.M, (2*gam22*yfp));
vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));


// Inner North-East Corner -----------------------------------------------------

node = 2*inputs.M-2.;

vars.B.set(node, node, (-gam22*(2*yfp+yfm)-3*gam11));
vars.B.set(node, node+1, (gam21*(yfp-yfm)+2*gam11));
vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-gam11));
vars.B.set(node, node+inputs.M, (gam22*yfm));
vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
vars.B.set(node, node-inputs.M, (2*gam22*yfp));
vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));



// Inner North Boundary ---------------------------------------------------------

gam11 = alph11/pow(inputs.dx,2);
gam22 = (4/3.)*alph22/pow(inputs.dy,2);
gam12 = (1/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (1/3.)*alph21/(inputs.dx*inputs.dy);

for (int i = inputs.M+2; i<2*inputs.M-2.;i++) {

	node = i;
	vars.B.set(node, node, (-gam22*(2*yfp+yfm)-2*gam11));
	vars.B.set(node, node+1, (gam21*(yfp-yfm)+gam11));
	vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-gam11));
	vars.B.set(node, node+inputs.M, (gam22*yfm));
	vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
	vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
	vars.B.set(node, node-inputs.M, (2*gam22*yfp));
	vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
	vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));	

}

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
// 					SOUTH AND INNER SOUTH BOUNDARIES				   //
// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //


// South-West Corner --------------------------------------------------

node = inputs.M*(inputs.N-1.);

if (inputs.axi == 0) {
	yfp = 1.0;
	yfm = 1.0;
}
else {
	yfp = (vars.y_vals[node]+inputs.dy/4.)/(vars.y_vals[node]+inputs.dy/8.);
	yfm = 0.0;
} 


gam11 = 8.*alph11/pow(inputs.dx,2);
gam22 = 8.*alph22/pow(inputs.dy,2);
gam12 = 4.*alph12/(inputs.dx*inputs.dy);
gam21 = 4.*alph21/(inputs.dx*inputs.dy);


vars.B.set(node,node,(-gam22*yfp-gam21*yfp-gam11-gam12));
vars.B.set(node,node+1.,(gam21*yfp+gam11-gam12));
vars.B.set(node,node-inputs.M,(gam22*yfp-gam21*yfp+gam12));
vars.B.set(node,node-inputs.M+1.,(gam12+gam21*yfp));


// South-East Corner --------------------------------------------------

node = inputs.M*inputs.N-1.;

vars.B.set(node,node,(-gam22*yfp+gam21*yfp-gam11+gam12));
vars.B.set(node,node-1.,(gam12+gam11-gam21*yfp));
vars.B.set(node,node-inputs.M,(gam22*yfp+gam21*yfp-gam12));
vars.B.set(node,node-inputs.M-1.,-1.*(gam12+gam21*yfp));



// South-South-West Corner --------------------------------------------------

node = inputs.M*(inputs.N-1.)+1.;

gam11 = (4/3.)*alph11/pow(inputs.dx,2);
gam22 = (8.)*alph22/pow(inputs.dy,2);
gam12 = (4/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (4/3.)*alph21/(inputs.dx*inputs.dy);


vars.B.set(node, node, (-gam22*yfp-3*gam11));
vars.B.set(node, node+1., (gam21*yfp+gam11-gam12));
vars.B.set(node, node-1., (2*gam11+gam12-gam21*yfp));
vars.B.set(node, node-inputs.M, (gam22*yfp));
vars.B.set(node, node-inputs.M+1., (gam12+gam21*yfp));
vars.B.set(node, node-inputs.M-1., -1.*(gam12+gam21*yfp));

// South-South-East Corner --------------------------------------------------

node = inputs.M*inputs.N-2.;


vars.B.set(node, node, (-gam22*yfp-3*gam11));
vars.B.set(node, node+1., (gam21*yfp+2*gam11+gam12));
vars.B.set(node, node-1., (gam11-gam12-gam21*yfp));
vars.B.set(node, node-inputs.M, (gam22*yfp));
vars.B.set(node, node-inputs.M+1., (gam12+gam21*yfp));
vars.B.set(node, node-inputs.M-1., -1.*(gam12+gam21*yfp));


// The rest of the South Boundary -------------------------------------------


gam11 = alph11/pow(inputs.dx,2);
gam22 = 8.*alph22/pow(inputs.dy,2);
gam12 = alph12/(inputs.dx*inputs.dy);
gam21 = alph21/(inputs.dx*inputs.dy);


for (int i = inputs.M*(inputs.N-1.)+2.; i < inputs.M*inputs.N-2.;i++) {
	
	node = i;

	vars.B.set(node, node, (-gam22*yfp-2*gam11));
	vars.B.set(node, node+1., (gam21*yfp+gam11-gam12));
	vars.B.set(node, node-1., (gam11+gam12-gam21*yfp));
	vars.B.set(node, node-inputs.M, (gam22*yfp));
	vars.B.set(node, node-inputs.M+1., (gam12+gam21*yfp));
	vars.B.set(node, node-inputs.M-1., -1.*(gam12+gam21*yfp));

}


// South-West-West Corner --------------------------------------------------

node = inputs.M*(inputs.N-2.);

if (inputs.axi == 0) {
	yfp = 1.0;
	yfm = 1.0;
}
else {
	yfp = (vars.y_vals[node]+inputs.dy/2.)/(vars.y_vals[node]+-inputs.dy/8.);
	yfm = (vars.y_vals[node]-inputs.dy/4.)/(vars.y_vals[node]+inputs.dy/8.);
} 

gam11 = (8.)*alph11/pow(inputs.dx,2);
gam22 = (4/3.)*alph22/pow(inputs.dy,2);
gam12 = (4/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (4/3.)*alph21/(inputs.dx*inputs.dy);


vars.B.set(node, node, (-gam21*(yfp-yfm)-gam22*(yfp+2*yfm)-gam11));
vars.B.set(node, node+1., (gam21*(yfp-yfm)+gam11));
vars.B.set(node, node-inputs.M, (gam22*yfp-gam21*yfp+gam12));
vars.B.set(node, node-inputs.M+1., (gam12+gam21*yfp));
vars.B.set(node, node+inputs.M, (2*gam22*yfm+gam21*yfm-gam12));
vars.B.set(node, node+inputs.M+1., -1.*(gam12+gam21*yfm));


// South-East-East Corner --------------------------------------------------

node = inputs.M*(inputs.N-1.)-1.;

vars.B.set(node, node, (-gam22*(yfp+2.*yfm)+gam21*(yfp-yfm)-gam11));
vars.B.set(node, node-1., (gam21*(yfm-yfp)+gam11));
vars.B.set(node, node-inputs.M, (gam22*yfp+gam21*yfp-gam12));
vars.B.set(node, node-inputs.M-1., -1.*(gam12+gam21*yfp));
vars.B.set(node, node+inputs.M, (2.*gam22*yfm-gam21*yfm+gam12));
vars.B.set(node, node+inputs.M-1., (gam12+gam21*yfm));


// Inside South-West Corner --------------------------------------------------

node = inputs.M*(inputs.N-2.)+1.;

gam11 = (4/3.)*alph11/pow(inputs.dx,2);
gam22 = (4/3.)*alph22/pow(inputs.dy,2);
gam12 = (4/9.)*alph12/(inputs.dx*inputs.dy);
gam21 = (4/9.)*alph21/(inputs.dx*inputs.dy);

vars.B.set(node, node, (-gam22*(yfp+2*yfm)-3*gam11));
vars.B.set(node, node+1, (gam21*(yfp-yfm)+gam11));
vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-2*gam11));
vars.B.set(node, node+inputs.M, (2*gam22*yfm));
vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
vars.B.set(node, node-inputs.M, (gam22*yfp));
vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));



// Inside South-East Corner --------------------------------------------------

node = inputs.M*(inputs.N-1.)-2.;

vars.B.set(node, node, (-gam22*(yfp+2*yfm)-3*gam11));
vars.B.set(node, node+1, (gam21*(yfp-yfm)+2.*gam11));
vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-gam11));
vars.B.set(node, node+inputs.M, (2.*gam22*yfm));
vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
vars.B.set(node, node-inputs.M, (gam22*yfp));
vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));



// Inner South Boundary ---------------------------------------------------------

gam11 = alph11/pow(inputs.dx,2);
gam22 = (4/3.)*alph22/pow(inputs.dy,2);
gam12 = (1/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (1/3.)*alph21/(inputs.dx*inputs.dy);



for (int i = inputs.M*(inputs.N-2)+2; i<inputs.M*(inputs.N-1)-2;i++) {

	node = i;

	vars.B.set(node, node, (-2*gam11-gam22*(yfp+2.*yfm)));
	vars.B.set(node, node+1, (gam21*(yfp-yfm)+gam11));
	vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-gam11));
	vars.B.set(node, node+inputs.M, (2.*gam22*yfm));
	vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
	vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
	vars.B.set(node, node-inputs.M, (gam22*yfp));
	vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
	vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));	

}




// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
// 						WEST AND INNER WEST BOUNDARIES				   //
// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //


// West Boundary --------------------------------------------------------

gam11 = 8.*alph11/pow(inputs.dx,2);
gam22 = alph22/pow(inputs.dy,2);
gam12 = alph12/(inputs.dx*inputs.dy);
gam21 = alph21/(inputs.dx*inputs.dy);

for (int i = inputs.M*2; i < inputs.M*(inputs.N-2.); i+=inputs.M) {
	
	node = i;
	
	if (inputs.axi == 0) {
		yfp = 1.0;
		yfm = 1.0;
	}
	else {
		yfp = (vars.y_vals[node]+inputs.dy/2.)/(vars.y_vals[node]);
		yfm = (vars.y_vals[node]-inputs.dy/2.)/(vars.y_vals[node]);
	} 


	vars.B.set(node, node, (-gam22*(yfp+yfm)-gam21*(yfp-yfm)-gam11));
	vars.B.set(node, node+1., (gam21*(yfp-yfm)+gam11));
	vars.B.set(node, node-inputs.M, (gam22*yfp-gam21*yfp+gam12));
	vars.B.set(node, node-inputs.M+1., (gam12+gam21*yfp));
	vars.B.set(node, node+inputs.M, (gam22*yfm+gam21*yfm-gam12));
	vars.B.set(node, node+inputs.M+1., -1.*(gam12+gam21*yfm));

}

// Inner West Boundary --------------------------------------------------

gam11 = (4/3.)*alph11/pow(inputs.dx,2);
gam22 = alph22/pow(inputs.dy,2);
gam12 = (1/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (1/3.)*alph21/(inputs.dx*inputs.dy);

for (int i = inputs.M*2+1.; i < inputs.M*(inputs.N-2.)+1.; i+=inputs.M) {
	
	node = i;
	
	if (inputs.axi == 0) {
		yfp = 1.0;
		yfm = 1.0;
	}
	else {
		yfp = (vars.y_vals[node]+inputs.dy/2.)/(vars.y_vals[node]);
		yfm = (vars.y_vals[node]-inputs.dy/2.)/(vars.y_vals[node]);
	} 


	vars.B.set(node, node, (-gam22*(yfp+yfm)-3*gam11));
	vars.B.set(node, node+1, (gam21*(yfp-yfm)+gam11));
	vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-2*gam11));
	vars.B.set(node, node+inputs.M, (gam22*yfm));
	vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
	vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
	vars.B.set(node, node-inputs.M, (gam22*yfp));
	vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
	vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));	

}


// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
// 						EAST AND INNER EAST BOUNDARIES				   //
// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //


// East Boundary --------------------------------------------------------

gam11 = 8.*alph11/pow(inputs.dx,2);
gam22 = alph22/pow(inputs.dy,2);
gam12 = alph12/(inputs.dx*inputs.dy);
gam21 = alph21/(inputs.dx*inputs.dy);

for (int i = inputs.M*3-1; i < inputs.M*(inputs.N-2.); i+=inputs.M) {
	
	node = i;
	
	if (inputs.axi == 0) {
		yfp = 1.0;
		yfm = 1.0;
	}
	else {
		yfp = (vars.y_vals[node]+inputs.dy/2.)/(vars.y_vals[node]);
		yfm = (vars.y_vals[node]-inputs.dy/2.)/(vars.y_vals[node]);
	} 


	vars.B.set(node, node, (-gam22*(yfp+yfm)+gam21*(yfp-yfm)-gam11));
	vars.B.set(node, node-1., -1.*(gam21*(yfp-yfm)-gam11));
	vars.B.set(node, node-inputs.M, (gam22*yfp+gam21*yfp-gam12));
	vars.B.set(node, node-inputs.M-1., -1.*(gam12+gam21*yfp));
	vars.B.set(node, node+inputs.M, (gam22*yfm-gam21*yfm+gam12));
	vars.B.set(node, node+inputs.M-1., (gam12+gam21*yfm));

}



// Inner East Boundary --------------------------------------------------

gam11 = (4/3.)*alph11/pow(inputs.dx,2);
gam22 = alph22/pow(inputs.dy,2);
gam12 = (1/3.)*alph12/(inputs.dx*inputs.dy);
gam21 = (1/3.)*alph21/(inputs.dx*inputs.dy);

for (int i = inputs.M*3-2; i < inputs.M*(inputs.N-2.)-1.; i+=inputs.M) {
	
	node = i;
	
	if (inputs.axi == 0) {
		yfp = 1.0;
		yfm = 1.0;
	}
	else {
		yfp = (vars.y_vals[node]+inputs.dy/2.)/(vars.y_vals[node]);
		yfm = (vars.y_vals[node]-inputs.dy/2.)/(vars.y_vals[node]);
	} 


	vars.B.set(node, node, (-gam22*(yfp+yfm)-3*gam11));
	vars.B.set(node, node+1, (gam21*(yfp-yfm)+2*gam11));
	vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-gam11));
	vars.B.set(node, node+inputs.M, (gam22*yfm));
	vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
	vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
	vars.B.set(node, node-inputs.M, (gam22*yfp));
	vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
	vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));	

}



// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
// 								EVERYTHING ELSE						   //
// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //


for (int i = inputs.M*2+2.; i<inputs.M*(inputs.N-2.)-2.; i++) {

	int num = (int) inputs.M;

	if ( i%num != 0 && i%num != 1 && i%num != num-1 && i%num != num-2 ) {

		node = i;
		
		gam11 = alph11/pow(inputs.dx,2);
		gam22 = alph22/pow(inputs.dy,2);
		gam12 = (1/4.)*alph12/(inputs.dx*inputs.dy);
		gam21 = (1/4.)*alph21/(inputs.dx*inputs.dy);	

		if (inputs.axi == 0) {
			yfp = 1.0;
			yfm = 1.0;
		}
		else {
			yfp = (vars.y_vals[node]+inputs.dy/2.)/(vars.y_vals[node]);
			yfm = (vars.y_vals[node]-inputs.dy/2.)/(vars.y_vals[node]);
		}

 
		vars.B.set(node, node, (-gam22*(yfp+yfm)-2*gam11));
		vars.B.set(node, node+1, (gam21*(yfp-yfm)+gam11));
		vars.B.set(node, node-1, -1.*(gam21*(yfp-yfm)-gam11));
		vars.B.set(node, node+inputs.M, (gam22*yfm));
		vars.B.set(node, node+inputs.M+1, -1.*(gam12+gam21*yfm));
		vars.B.set(node, node+inputs.M-1, (gam12+yfm*gam21));
		vars.B.set(node, node-inputs.M, (gam22*yfp));
		vars.B.set(node, node-inputs.M+1, (gam12+yfp*gam21));
		vars.B.set(node, node-inputs.M-1, -1.*(gam12+yfp*gam21));	
	}
}

/*
	-------------------------------------------------------------------------
	-------------------------------------------------------------------------

	The second part of this code writes the [C] matrix. It is just the identity
	matrix of dims M*N x M*N.

	-------------------------------------------------------------------------
	
	\author Jared Clifford
	\version October 2013

	-------------------------------------------------------------------------
	-------------------------------------------------------------------------

*/

	for (int i = 0; i< inputs.M*inputs.N; i++) {
		for (int j = 0; j< inputs.M*inputs.N; j++ ) {
			if (i == j) {
				vars.C.set(i, j, 1.0);
			}
		}	
	}

/*
	-------------------------------------------------------------------------
	-------------------------------------------------------------------------

	The last step is to find [A] = dt*[B] + [C].

	15/10/2013 for the moment, inputs.dt is used as a fixed timestep.  This
	will be updated in future tests to account for the varying timestep.

	-------------------------------------------------------------------------
	
	\author Jared Clifford
	\version October 2013

	-------------------------------------------------------------------------
	-------------------------------------------------------------------------

*/

	for (int i = 0; i< inputs.M*inputs.N; i++) {
		for (int j = 0; j< inputs.M*inputs.N; j++ ) {
			vars.A.set(i, j, vars.B.get(i,j)*inputs.dt + vars.C.get(i,j));			
		}	
	}
	return SUCCESS;
}


