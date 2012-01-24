/*  \file reaction_system_test.cxx
 *  \brief Testing program for the finite-rate chemistry module.
 *
 *  \author Rowan J Gollan
 *  \version 25-Feb-2006
 *
 */

/** \file tpgm_test.cxx
 * \brief Testing program of the thermally perfect gas mix.
 * \author RJG
 * \version 14-Feb-06
 *
 **/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "../../util/source/config_parser.hh"
#include "../../util/source/useful.h"
#include "physical_constants.hh"
#include "gas.hh"
//#include "rr_coeffs.hh"
//#include "reaction.hh"
//#include "../../nm/source/ode_solver.hh"
#include "reaction_scheme.hh"


using namespace std;

void printUsage()
{
    cout << "Usage: chresspi.x input_file\n";
    return;	
}


void perform_single_step( const string infile, Gas_data &Q, double step_interval );
void integrate_to_final_time( const string infile, Gas_data &Q, double dt_final, double dt_plot,
			      const string outfile, const string output_type, int plot_every );

int main( int argc, char *argv[] )
{
    if( argc != 2 ) {
	printUsage();
	return 0;
    }

    string infile( argv[1] );

    ConfigParser cfg( infile );

    //--- Parse options from gas_state options ---//

    string gas_model;
    string nf("notfound");
    if( ! cfg.parse_string( "gas_state", "gas_model", gas_model, nf ) ) {
	cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
	     << "Error reading gas_model in section [gas_state] of file:"
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string input_file;
    if( ! cfg.parse_string( "gas_state", "input_file", input_file, nf ) ) {
	cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
	     << "Error reading input_file in section [gas_state] of file:"
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    // Now we can set the gas model (also set mf array sizes)
    set_type_of_gas( gas_model, input_file );
    Gas_model *g = get_gas_model_pointer();
    Gas_data Q(g);
    set_array_sizes_in_gas_data(Q, g->nsp, 0);

    // Now grab the gas state.
    if( ! cfg.parse_double( "gas_state", "p", Q.p, 0.0) ) {
	cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
	     << "Error reading p in section [gas_state] of file: "
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( ! cfg.parse_double( "gas_state", "T", Q.T, 0.0) ) {
	cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
	     << "Error reading T in section [gas_state] of file: "
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    ostringstream mf;
    for( int isp = 0; isp < g->nsp; ++isp ) {
	mf.str("");
	mf << "f-" << isp;
	if( ! cfg.parse_double( "gas_state", mf.str(), Q.f[isp], 0.0) ) {
	    cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error reading " << mf << " in section [gas_state] of file: "
		 << infile << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
    }
    Q.dt_chem = 0.0; Q.dt_therm = 0.0;

    g->EOS_pT( Q, true );
   //  cout << "Initial conditions...\n";
//     print_gas_data( Q );

    // --- Parse options from control section --- //
    if( ! cfg.parse_double( "control", "dt_initial", Q.dt_chem, -1.0 ) ) {
	cerr << "Error reading dt_initial in section [control] of file: "
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    string reac_file;
    if( ! cfg.parse_string( "control", "reac_file", reac_file, "reactions.rsi") ) {
	cerr << "Error reading reac_file in section [control] of file: "
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    set_reaction_scheme( "chresspi_system", reac_file); 

    string action;
    if( ! cfg.parse_string( "control", "action", action, nf) ) {
	cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
	     << "Error reading action in [control] section of file: "
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    if( action == "single_step" ) {
	double step_interval;
	if( ! cfg.parse_double( "control", "step_interval", step_interval, -1.0 ) ) {
	     cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error reading step_interval in section [control] of file: "
		 << infile << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	if( step_interval <= 0.0 ) {
	    cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error step_interval= " << step_interval << endl
		 << "This is less than or equal to zero, and therefore an invalid double\n"
		 << "for the step interval (should be positive and finite).\n"
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}

	perform_single_step( infile, Q, step_interval );
    }
    else if( action == "integrate_to_final_time" ) {
	string outfile;
	string odefault( "chresspi.out.data" );
	if( ! cfg.parse_string( "control", "output_file", outfile, odefault ) ) {
	    cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error reading out_file in section [control] of file: "
		 << infile << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	string output_type;
	if( ! cfg.parse_string( "control", "output_type", output_type, "mass_f" ) ) {
	    cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error reading output_type in section [control] of file: "
		 << infile << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}

	double dt_final;
	if( ! cfg.parse_double( "control", "dt_final", dt_final, -1.0 ) ) {
	    cerr << "Error reading dt_final in section [control] of file: "
		 << infile << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	if( dt_final <= 0.0 ) {
	    cerr << "Error dt_final= " << dt_final << endl
		 << "This is less than or equal to zero, and therefore an invalid double\n"
		 << "for integrating to final time.\n"
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	double dt_thermo;
	if( ! cfg.parse_double( "control", "dt_thermo", dt_thermo, -1.0 ) ) {
	    cerr << "Error reading dt_thermo in section [control] of file: "
		 << infile << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	if( dt_thermo <= 0.0 ) {
	    cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error dt_thermo= " << dt_thermo << endl
		 << "This is less than or equal to zero, and therefore an invalid double\n"
		 << "for plotting intervals.\n"
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}

	if( dt_thermo > dt_final ) {
	    cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error: dt_plot= " << dt_thermo << " and dt_final= " << dt_final << endl
		 << "dt_thermo is greater than dt_final and so no data would be plotted.\n"
		 << "This is probably not what you want.\n"
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}

	int plot_every;
	if( ! cfg.parse_int( "control", "plot_every", plot_every, 10 ) ) {
	    cerr << "Error reading plot_every in section [control] of file: "
		 << infile << endl
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}
	
	if( plot_every <= 0 ) {
	    cerr << __FILE__ << ":main()(" << (__LINE__ - 1) << endl
		 << "Error plot_every= " << plot_every << endl
		 << "This is less than or equal to zero, and therefore an invalid int for then\n"
		 << "for plotting intervals.\n"
		 << "Bailing out!\n";
	    exit(BAD_INPUT_ERROR);
	}

	
	integrate_to_final_time( infile, Q, dt_final, dt_thermo, outfile, output_type, plot_every );
    }
    else {
	cerr << "Error: unknown action= " << action << " listed in [control] of file: "
	     << infile << endl
	     << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    clear_gas_model_pointer();
    clear_reaction_scheme_pointer();
    return 0;

}

void perform_single_step( const string infile, Gas_data &Q, double step_interval )
{
    cout << "Input gas state...\n";
    print_gas_data( Q );

    cout << "Perform integration over single interval: dt= " <<  step_interval << endl << endl;
    
    perform_chemical_increment( Q, step_interval );

    cout << "Output gas state...\n";
    print_gas_data( Q );

}

void integrate_to_final_time( const string infile, Gas_data &Q, double dt_final, double dt_thermo,
			      const string outfile, const string output_type, int plot_every )
{
    Gas_model *g = get_gas_model_pointer();

    double t = 0.0;
    ofstream of( outfile.c_str() );
    
    of << setprecision(10) << showpoint;
    of << "# " << outfile << endl
       << "# Columns:\n"
       << "# 1: t \n"
       << "# 2: T \n"
       << "# 3: p \n";
    if( output_type == "mass_f" ) {
	for( int isp = 0; isp < g->nsp; ++isp ) of << "# " << isp+4 << ": f[" << isp << "]\n";
    }
    else if( output_type == "mole_f" ) {
	for( int isp = 0; isp < g->nsp; ++isp ) of << "# " << isp+4 << ": mole_f[" << isp << "]\n";
    }
    else if( output_type == "moles" ) {
	for( int isp = 0; isp < g->nsp; ++isp ) of << "# " << isp+4 << ": moles[" << isp << "]\n";
    }
    else if( output_type == "number_density" ) {
	for( int isp = 0; isp < g->nsp; ++isp ) of << "# " << isp+4 << ": nd[" << isp << "]\n";
    }
    else {
	cout << "Output type: " << output_type << " is unknown - defaulting to mass fractions.\n";
	for( int isp = 0; isp < g->nsp; ++isp ) of << "# " << isp+4 << ": f[" << isp << "]\n";
    }
    int count = 1;
    while( t < dt_final ) {
	perform_chemical_increment( Q, dt_thermo);
	t += dt_thermo;
	if( count % plot_every == 0 ) {
	    of << t << "  " << Q.T << "  " << Q.p;
	    if( output_type == "mass_f" ) {
		for( int isp = 0; isp < g->nsp; ++isp ) of << " " << Q.f[isp];
	    }
	    else if( output_type == "mole_f" ) {
		for( int isp = 0; isp < g->nsp; ++isp ) of << " " << Q.c[isp]/Q.c_tot;
	    }
	    else if( output_type == "moles" ) {
		for( int isp = 0; isp < g->nsp; ++isp ) of << " " << Q.c[isp];
	    }
	    else if( output_type == "number_density" ) {
		for( int isp = 0; isp < g->nsp; ++isp ) of << " " << Q.c[isp]*PC_Na;
	    }
	    else { // Just print mass_f
		for( int isp = 0; isp < g->nsp; ++isp ) of << " " << Q.f[isp];
	    }
	    of << endl;
	}
	if( count % 1000 == 0 ) cout << "count: " << count << endl;
	count = count + 1;
    }

   //  cout << "Gas state at end of integration...\n";
//     print_gas_data( Q );

    of.close();
  //   cout << "Single-point integration finished at t= " << t << endl;
//     cout << "Thermodynamic properties were re-evaluated every dt= " << dt_thermo << endl;
//     cout << "File: " << outfile << " created.\n";
//     cout << "Done.\n";

    return;
}



