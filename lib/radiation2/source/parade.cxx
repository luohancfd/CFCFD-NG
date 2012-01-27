/** \file photaura.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 28-Jul-09 : initial framework
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

#include "parade.hh"
#include "radiation_constants.hh"
#include "../../util/source/lua_service.hh"

using namespace std;

Parade::
Parade( lua_State * L )
 : RadiationSpectralModel( L )
{
    this->initialise(L);
}

Parade::
Parade( const string input_file )
 : RadiationSpectralModel( input_file )
{
    // 1. Get spectral_model string from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, input_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Parade::Parade()\n";
	ost << "Error in input file: " << input_file << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "spectral_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Parade::Parade():\n";
	ost << "Error in the declaration of spectral_data: a table is expected.\n";
	input_error(ost);
    }
    
    this->initialise(L);

    lua_pop(L,1);	// pop spectral_data
    lua_close(L);
}

void Parade::initialise( lua_State * L )
{
    string parade_template_filename = get_string(L, -1, "parade_template");

    if ( ECHO_RAD_INPUT > 0 )
	cout << "parade_template_filename: " << parade_template_filename << endl;

    read_parade_template_file( parade_template_filename );

    iT  = get_int(L,-1,"iT");
    iTe = get_int(L,-1,"iTe");

    lua_getfield(L, -1, "radiators" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Parade::Parade():\n";
	ost << "Error in the declaration of radiators: a table is expected.\n";
	input_error(ost);
    }

    nrad = lua_objlen(L, -1);

    if ( nrad==0 ) {
	cout << "Parade::Parade()" << endl
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
	radiators.push_back( create_new_parade_radiator( L, rad_names[irad] ) );
	if ( radiators.back()->type == "electron_radiator" ) {
	    e_index = radiators.back()->isp;
	}
    }
}

Parade::~Parade()
{
    // 1. Delete radiators
    for ( size_t irad=0; irad<radiators.size(); ++irad )
	delete radiators[irad];
}

string Parade::str() const
{
    return "Parade";
}

double
Parade::
integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved )
{
    double j_total = 0.0;
    
    return j_total;
}

double
Parade::
variably_integrated_emission_for_gas_state( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u, bool spectrally_resolved )
{
    cout << "Parade::variably_integrated_emission_for_gas_state()" << endl
         << "This function is not available for the parade radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Parade::
spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X )
{
    // 0. Make sure the vectors in CoeffSpectra are sized to zero.
    //    This is required for adaptive spectral distributions.
    X.nu.clear();
    X.j_nu.clear();
    X.kappa_nu.clear();

    // 1. Create the parade control file
    create_parade_control_file( Q );
    
    // 2. Run the parade executable
    system("parade > parade.out");

    // 3. Pick up the solution and insert it into CoeffSpectra
    ifstream specfile( "par_res.dat" );
    if ( !specfile.is_open() ) {
	cout << "Parade::spectra_for_gas_state()" << endl
             << "Could not open parade spectra file 'par_res.dat'." << endl
             << "Exiting program." << endl;
	exit( FAILURE );
    }

    string line;
    // First five lines are headers
    for ( int i=0; i<5; ++i ) {
	getline( specfile, line );
	// cout << line << endl;
    }

    // Remaining lines should be spectral data
    // Note that the parade data starts from the lower wavelength whereas the
    // CoeffSpectra class starts from the highest wavelength, and parade outputs
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
Parade::
spectral_distribution_for_gas_state(Gas_data &Q, vector<double> &nus)
{
    cout << "Parade::spectral_distribution_for_gas_state()" << endl
         << "This function is not available for the parade radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Parade::
write_line_widths( Gas_data &Q )
{
    cout << "Parade::write_line_widths()" << endl
         << "This function is not available for the parade radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Parade::
prep_rad_pop_files()
{
    cout << "Parade::prep_rad_pop_files()" << endl
         << "This function is not available for the parade radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Parade::
append_current_rad_pops( double x )
{
    cout << "Parade::append_current_rad_pops()" << endl
         << "This function is not available for the parade radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Parade::
write_QSS_analysis_files( Gas_data &Q, int index )
{
    cout << "Parade::write_QSS_analysis_files()" << endl
         << "This function is not available for the parade radiation model." << endl
         << "Exiting program!" << endl;
    exit( BAD_INPUT_ERROR );
}

void
Parade::
read_parade_template_file( string parade_template_filename )
{
    ifstream ptfile( parade_template_filename.c_str() );

    parade_template_file_buffer << ptfile.rdbuf();

    ptfile.close();

    if ( ECHO_RAD_INPUT>0 ) {
        cout << "template file:" << endl
	     << parade_template_file_buffer.str()
	     << endl;
    }
}

void
Parade::
create_parade_control_file( Gas_data &Q )
{
    ofstream pcfile( "parade.con" );

    pcfile << parade_template_file_buffer.str();

    // Now put down Tt, Te
    pcfile << " " << Q.T[iT]  << endl
	   << " " << Q.T[iTe] << endl;

    // Now put down the partial densities and Tr and Tv for molecules
    for ( int irad=0; irad<nrad; ++irad ) {
	ParadeRadiator * R = radiators[irad];
	double number_density = Q.rho * Q.massf[R->isp] / R->m_w * RC_Na;
	if ( number_density < 1.0 ) number_density = 1.0;
	pcfile << number_density;
	if ( R->type=="diatomic_radiator" || R->type=="triatomic_radiator" )
	    pcfile << "   " << Q.T[R->iTr] << "   " << Q.T[R->iTv];
	pcfile << endl;
    }

    pcfile.close();
}
