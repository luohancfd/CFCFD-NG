/** \file parade.cxx
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
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

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
    // system return value
    int srv;

    string control_template_filename = get_string(L, -1, "control_template");

    if ( ECHO_RAD_INPUT > 0 )
	cout << "control_template_filename: " << control_template_filename << endl;

    read_control_template_file( control_template_filename );

    string snbopt_template_filename = get_string(L, -1, "snbopt_template");

    if ( ECHO_RAD_INPUT > 0 )
        cout << "snbopt_template_filename: " << snbopt_template_filename << endl;

    if ( snbopt_template_filename.compare("none")!=0 )
        read_snbopt_template_file( snbopt_template_filename );

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

    // create working directories for each thread
    if ( omp_get_thread_num()==0 ) {
        for ( int i=0; i<omp_get_max_threads(); ++i ) {
            ostringstream pathname;
            pathname << "parade_working_dir_" << i;
            if ( access(pathname.str().c_str(), F_OK) != 0 ) {
		ostringstream cmd;
		cmd << "mkdir " << pathname.str();
		srv = system(cmd.str().c_str());
            }
        }
    }

    UNUSED_VARIABLE(srv);
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

ParadeRadiator *
Parade::
get_radiator_pointer_from_name( string name )
{
    ParadeRadiator * R = 0;

    for ( int irad=0; irad<nrad; ++irad ) {
        if ( radiators[irad]->name == name ) {
            R = radiators[irad];
            break;
        }
    }

    if ( !R ) {
        cout << "Parade::get_radiator_pointer_from_name()" << endl
             << "Radiator with name: " << name << " not found." << endl
             << "Exiting program." << endl;
        exit( BAD_INPUT_ERROR );
    }

    return R;
}

double
Parade::
integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved )
{
    CoeffSpectra X;
    this->spectra_for_gas_state(Q,X);
    
    return X.integrate_emission_spectra();
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
    // system return value
    int srv;
    
    // Move to the working directory
    ostringstream oss;
    oss << "parade_working_dir_" << omp_get_thread_num();
    srv = chdir(oss.str().c_str());

    // Set the adaptive flag
    X.adaptive = adaptive_spectral_grid;

    // 0. Make sure the vectors in CoeffSpectra are sized to zero.
    //    This is required for adaptive spectral distributions.
    X.nu.clear();
    X.j_nu.clear();
    X.kappa_nu.clear();

    // 1. Create the parade control file
    create_parade_control_files( Q );
    
    // 2. Run the parade executable
    srv = system("parade > parade.out 2>&1");

    if ( srv!=0 ) {
        cout << "Parade::spectra_for_gas_state()" << endl
             << "Something went wrong with the system call." << endl;
        exit( FAILURE );
    }

    // 3. Pick up the solution and insert it into CoeffSpectra
    ifstream specfile( "par_res.txt" );
    if ( !specfile.is_open() ) {
	cout << "Parade::spectra_for_gas_state()" << endl
             << "Could not open parade spectra file 'par_res.txt'." << endl
             << "Exiting program." << endl;
	exit( FAILURE );
    }

    // Remaining lines should be spectral data
    // Note that the parade data starts from the lower wavelength whereas the
    // CoeffSpectra class starts from the highest wavelength, and parade outputs
    // j_lambda whereas we need j_nu (hence the conversion)
    double lambda_ang, nu, j_lambda, kappa;
    while ( specfile >> lambda_ang >> j_lambda >> kappa ) {
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

    // print the size of the spectral vectors
#   if DEBUG_RAD > 0
    cout << "X.nu.size() = " << X.nu.size() << endl;
#   endif

    // Move out of the working directory
    srv = chdir("..");

    UNUSED_VARIABLE( srv );

    return;
}

void
Parade::
write_line_widths( Gas_data &Q )
{}

void
Parade::
prep_rad_pop_files()
{}

void
Parade::
append_current_rad_pops( double x )
{}

void
Parade::
write_QSS_analysis_files( Gas_data &Q, int index )
{}

void
Parade::
read_control_template_file( string control_template_filename )
{
    ifstream ptfile( control_template_filename.c_str() );

    control_template_file_buffer << ptfile.rdbuf();

    ptfile.close();

    if ( ECHO_RAD_INPUT>1 ) {
        cout << "template file:" << endl
	     << control_template_file_buffer.str()
	     << endl;
    }
}

void
Parade::
read_snbopt_template_file( string snbopt_template_filename )
{
    ifstream sotfile( snbopt_template_filename.c_str() );

    snbopt_template_file_buffer << sotfile.rdbuf();

    sotfile.close();

    if ( ECHO_RAD_INPUT>1 ) {
        cout << "template file:" << endl
             << snbopt_template_file_buffer.str()
             << endl;
    }
}

void
Parade::
create_parade_control_files( Gas_data &Q )
{
#   if USE_FLO_INPUT_FILES
    // Flowfield style: gas data in the .flow files
    // firstly the con file
    // ofstream pcfile( "parade.con" );  // commented out these 3 lines Sept 2014.
    // pcfile << control_template_file_buffer.str();
    // pcfile.close();
    // Grid file
    // FIXME: this should be done as a initialisation step
    ofstream gfile( "grid.flo" );
    gfile << "TINA" << endl
          << "           1          1" << endl
          << "      0.00000000000000000E+00      0.00000000000000000E+00" << endl;
    gfile.close();
    // Number density file
    ofstream dfile( "dens.flo" );
    dfile << "           1          1          " << nrad << endl;
    for ( int irad=0; irad<nrad; ++irad ) {
	ParadeRadiator * R = radiators[irad];
	double number_density = Q.rho * Q.massf[R->isp] / R->m_w * RC_Na;
	if ( number_density < 1.0 ) number_density = 1.0;
	dfile << " " << number_density;
    }
    dfile << endl;
    dfile.close();
    // Temperature file
    ofstream tfile( "temp.flo" );
    //determine number of temperatures
    // FIXME: this should be done as a initialisation step
    int ntemps = 2;
    for ( int irad=0; irad<nrad; ++irad ) {
	ParadeRadiator * R = radiators[irad];
	if ( R->type=="diatomic_radiator" || R->type=="triatomic_radiator" )
            ntemps += 2;
    }
    tfile << "           1          1          " << ntemps << endl;
    // Now put down Tt, Te
    tfile << " " << Q.T[iT]  << " " << Q.T[iTe];
    // Now put down Tv and Tr for each molecule
    for ( int irad=0; irad<nrad; ++irad ) {
	ParadeRadiator * R = radiators[irad];
	if ( R->type=="diatomic_radiator" || R->type=="triatomic_radiator" )
	    tfile << "   " << Q.T[R->iTr] << "   " << Q.T[R->iTv];
    }
    tfile << endl;
    tfile.close();
#   else
    // Single cell style: all gas data in the .con file
    ofstream pcfile( "parade.con" );
    pcfile << control_template_file_buffer.str();
    // Firstly put down Tt, Te
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
#   endif

    if ( snbopt_template_file_buffer.str().length()>0 ) {
        ofstream sofile( "SNBOPT" );
        sofile << snbopt_template_file_buffer.str();
        sofile.close();
    }
}
