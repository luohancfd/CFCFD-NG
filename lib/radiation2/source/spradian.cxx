/** \file spradian.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 8-Mar-12 : initial framework
 *
 **/
 
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <algorithm>

#include "spradian.hh"
#include "radiation_constants.hh"
#include "../../util/source/lua_service.hh"

using namespace std;

Spradian::
Spradian( lua_State * L )
 : RadiationSpectralModel( L )
{
    this->initialise(L);
}

Spradian::
Spradian( const string input_file )
 : RadiationSpectralModel( input_file )
{
    // 1. Get spectral_model string from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, input_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Spradian::Spradian()\n";
	ost << "Error in input file: " << input_file << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "spectral_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Spradian::Spradian():\n";
	ost << "Error in the declaration of spectral_data: a table is expected.\n";
	input_error(ost);
    }
    
    this->initialise(L);

    lua_pop(L,1);	// pop spectral_data
    lua_close(L);
}

void Spradian::initialise( lua_State * L )
{
    string spradian_template_filename = get_string(L, -1, "spradian_template");

    if ( ECHO_RAD_INPUT > 0 )
	cout << "spradian_template_filename: " << spradian_template_filename << endl;

    read_spradian_template_file( spradian_template_filename );

    iT  = get_int(L,-1,"iT");
    iTe = get_int(L,-1,"iTe");

    lua_getfield(L, -1, "radiators" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Spradian::Spradian():\n";
	ost << "Error in the declaration of radiators: a table is expected.\n";
	input_error(ost);
    }

    nrad = lua_objlen(L, -1);

    if ( nrad==0 ) {
	cout << "Spradian::Spradian()" << endl
             << "No radiators have been requested, exiting program" << endl;
        exit( BAD_INPUT_ERROR );
    }

    for ( int irad = 0; irad < nrad; ++irad ) {
	lua_rawgeti(L, -1, irad+1); // A Lua list is offset one from the C++ vector index
	const char* rad = luaL_checkstring(L, -1);
	rad_names.push_back(string(rad));
	lua_pop(L, 1);
    }

    lua_pop(L,1);	// pop radiators

    // Now construct the radiators
    for ( int irad = 0; irad < nrad; ++irad ) {
	radiators.push_back( create_new_spradian_radiator( L, rad_names[irad] ) );
	if ( radiators.back()->type == "electron_radiator" ) {
	    e_index = radiators.back()->isp;
	}
    }
}

Spradian::~Spradian()
{
    // 1. Delete radiators
    for ( size_t irad=0; irad<radiators.size(); ++irad )
	delete radiators[irad];
}

string Spradian::str() const
{
    return "Spradian";
}

double
Spradian::
integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved )
{
    double j_total = 0.0;
    
    return j_total;
}

double
Spradian::
variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved )
{
    cout << "Spradian::variably_integrated_emission_for_gas_state()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X )
{
    // 0. Make sure the vectors in CoeffSpectra are sized to zero.
    //    This is required for adaptive spectral distributions.
    X.nu.clear();
    X.j_nu.clear();
    X.kappa_nu.clear();

    // 1. Create the spradian control file
    create_spradian_control_file( Q );
    
    // 2. Run the spradian executable
    system("spradian > spradian.out");

    // 3. Pick up the solution and insert it into CoeffSpectra
    ifstream specfile( "par_res.txt" );
    if ( !specfile.is_open() ) {
	cout << "Spradian::spectra_for_gas_state()" << endl
             << "Could not open spradian spectra file 'par_res.dat'." << endl
             << "Exiting program." << endl;
	exit( FAILURE );
    }

    // Remaining lines should be spectral data
    // Note that the spradian data starts from the lower wavelength whereas the
    // CoeffSpectra class starts from the highest wavelength, and spradian outputs
    // j_lambda whereas we need j_nu (hence the conversion)
    double lambda_ang, nu, j_lambda, kappa;
    while ( specfile.good() ) {
        specfile >> lambda_ang >> j_lambda >> kappa;
        // cout << "lambda_ang = " << lambda_ang << ", j_lambda = " << j_lambda << ", kappa = " << kappa << endl;
        nu = lambda2nu( lambda_ang/10.0 );
        X.nu.push_back( nu );
        X.j_nu.push_back( j_lambda * RC_c_SI / nu / nu );
        X.kappa_nu.push_back( kappa );
    }

    specfile.close();

    // We want ascending frequencies for consistency with photaura model
    if ( X.nu.front() > X.nu.back() ) {
        reverse( X.nu.begin(), X.nu.end() );
        reverse( X.j_nu.begin(), X.j_nu.end() );
        reverse( X.kappa_nu.begin(), X.kappa_nu.end() );
    }

    return;
}

void
Spradian::
spectral_distribution_for_gas_state(Gas_data &Q, vector<double> &nus)
{
    cout << "Spradian::spectral_distribution_for_gas_state()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
write_line_widths( Gas_data &Q )
{
    cout << "Spradian::write_line_widths()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
prep_rad_pop_files()
{
    cout << "Spradian::prep_rad_pop_files()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
append_current_rad_pops( double x )
{
    cout << "Spradian::append_current_rad_pops()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
write_QSS_analysis_files( Gas_data &Q, int index )
{
    cout << "Spradian::write_QSS_analysis_files()" << endl
         << "This function is not available for the spradian radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Spradian::
read_spradian_template_file( string spradian_template_filename )
{
    ifstream ptfile( spradian_template_filename.c_str() );

    spradian_template_file_buffer << ptfile.rdbuf();

    ptfile.close();

    if ( ECHO_RAD_INPUT>1 ) {
        cout << "template file:" << endl
	     << spradian_template_file_buffer.str()
	     << endl;
    }
}

void
Spradian::
create_spradian_control_file( Gas_data &Q )
{
    ofstream pcfile( "spradian.con" );

    pcfile << spradian_template_file_buffer.str();

    // Now put down Tt, Te
    pcfile << " " << Q.T[iT]  << endl
	   << " " << Q.T[iTe] << endl;

    // Now put down the partial densities and Tr and Tv for molecules
    for ( int irad=0; irad<nrad; ++irad ) {
	SpradianRadiator * R = radiators[irad];
	double number_density = Q.rho * Q.massf[R->isp] / R->m_w * RC_Na;
	if ( number_density < 1.0 ) number_density = 1.0;
	pcfile << number_density;
	if ( R->type=="diatomic_radiator" || R->type=="triatomic_radiator" )
	    pcfile << "   " << Q.T[R->iTr] << "   " << Q.T[R->iTv];
	pcfile << endl;
    }

    pcfile.close();
}
