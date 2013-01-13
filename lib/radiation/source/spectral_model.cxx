/** \file spectral_model.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 10-July-09
 *
 **/

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <math.h>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "spectral_model.hh"
#include "equilibrium_air.hh"
#include "photaura.hh"
#if WITH_SPRADIAN == 1
#include "spradian.hh"
#endif
#include "parade.hh"
#include "radiation_constants.hh"

using namespace std;

RadiationSpectralModel::
RadiationSpectralModel() {}

RadiationSpectralModel::
RadiationSpectralModel( lua_State * L )
{
    lua_getglobal(L, "spectral_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "RadiationSpectralModel::RadiationSpectralModel():\n";
	ost << "Error in the declaration of spectral_data: a table is expected.\n";
	input_error(ost);
    }
    
    lambda_min = get_positive_number(L, -1, "lambda_min");
    lambda_max = get_positive_number(L, -1, "lambda_max");
    
    // Sanity check
    if ( lambda_min >= lambda_max ) {
    	ostringstream ost;
	ost << "RadiationSpectralModel::RadiationSpectralModel():\n";
	ost << "Error in the declaration of spectral limits\n";
	ost << "lambda_min = " << lambda_min << " > lambda_max = " << lambda_max
	    << "\n";
	input_error(ost);
    }
    
    spectral_points = get_positive_int( L, -1, "spectral_points" );
    spectral_blocks = get_positive_int( L, -1, "spectral_blocks" );
    
    if ( ECHO_RAD_INPUT > 0 ) {
    	cout << "lambda_min = " << lambda_min << " nm" << endl
    	     << "lambda_max = " << lambda_max << " nm" << endl
    	     << "spectral_points = " << spectral_points << endl
    	     << "spectral_blocks = " << spectral_blocks << endl;
    }
    
    reset_spectral_params(0);
    
    lua_pop(L,1);	// pop spectral_data
}

RadiationSpectralModel::
RadiationSpectralModel( const string input_file )
{
    // 1. Get spectral_model string from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, input_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "RadiationSpectralModel::RadiationSpectralModel()\n";
	ost << "Error in input file: " << input_file << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "spectral_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "RadiationSpectralModel::RadiationSpectralModel():\n";
	ost << "Error in the declaration of spectral_data: a table is expected.\n";
	input_error(ost);
    }
    
    lambda_min = get_positive_number(L, -1, "lambda_min");
    lambda_max = get_positive_number(L, -1, "lambda_max");
    
    // Sanity check
    if ( lambda_min >= lambda_max ) {
    	ostringstream ost;
	ost << "RadiationSpectralModel::RadiationSpectralModel():\n";
	ost << "Error in the declaration of spectral limits\n";
	ost << "lambda_min = " << lambda_min << " > lambda_max = " << lambda_max
	    << "\n";
	input_error(ost);
    }
    
    spectral_points = get_positive_int( L, -1, "spectral_points" );
    spectral_blocks = get_positive_int( L, -1, "spectral_blocks" );
    
    if ( ECHO_RAD_INPUT > 0 ) {
    	cout << "lambda_min = " << lambda_min << " nm" << endl
    	     << "lambda_max = " << lambda_max << " nm" << endl
    	     << "spectral_points = " << spectral_points << endl
    	     << "spectral_blocks = " << spectral_blocks << endl;
    }
    
    reset_spectral_params(0);
    
    lua_pop(L,1);	// pop spectral_data
    
    lua_pop(L,1);	// pop spectral_data
    lua_close(L);
}

RadiationSpectralModel::
~RadiationSpectralModel() {}



void
RadiationSpectralModel::
reset_spectral_params( int isb )
{
    // If there is only one spectral block just copy the global range
    if ( spectral_blocks==1 ) {
	lambda_min_star = lambda_min;
	lambda_max_star = lambda_max;
	spectral_points_star = spectral_points;
	spectral_block = 0;
	return;
    }
    
    double sps_av = double(spectral_points) / double(spectral_blocks);
    // If this is the first block make it the largest with the remainder points
    // NOTE: this needs to be done so as to make sure buffers are large enough
    if ( isb==0 ) {
	spectral_points_star = spectral_points - ( spectral_blocks - 1 )*int(sps_av);
    }
    else {
	spectral_points_star = int(sps_av);
    }
    
    // Calculate upper and lower wavelengths
    // NOTE: need to account for odd size of first block
    int inu_lower = 0;
    if ( isb>0 ) inu_lower = spectral_points - ( spectral_blocks - 1 )*int(sps_av) + isb*int(sps_av);
    int inu_upper = inu_lower + spectral_points_star - 1;
    
    double dnu = ( lambda2nu(lambda_min) - lambda2nu(lambda_max) ) 
		/ double ( spectral_points - 1 );
    delta_nu = dnu;
    double nu_min_star = lambda2nu(lambda_max) + inu_lower * dnu;
    double nu_max_star = lambda2nu(lambda_max) + inu_upper * dnu;
    
    // Set lambda_min_star and lambda_max_star
    // NOTE: lambda_min uses nu_max and vice-versa
    lambda_max_star = nu2lambda(nu_min_star);
    lambda_min_star = nu2lambda(nu_max_star);
    
    // Set the spectral block index
    spectral_block = isb;
    
    // cout << "lambda_min_star = " << lambda_min_star << ", lambda_max_star = " << lambda_max_star << endl;
    // cout << "spectral_points_star = " << spectral_points_star << endl;
    
    return;
}

void
RadiationSpectralModel::
prep_rad_pop_files()
{
    cout << "RadiationSpectralModel::prep_rad_pop_files()" << endl
         << "This function has No functionality in base class" << endl
         << "Bailing out!" << endl;
    exit( FAILURE );
}

void
RadiationSpectralModel::
append_current_rad_pops( double x )
{
    cout << "RadiationSpectralModel::append_current_rad_pops()" << endl
         << "This function has No functionality in base class" << endl
         << "Bailing out!" << endl;
    exit( FAILURE );
}

void
RadiationSpectralModel::
write_QSS_analysis_files( Gas_data &Q, int index )
{
    cout << "RadiationSpectralModel::write_QSS_analysis_files()" << endl
         << "This function has No functionality in base class" << endl
         << "Bailing out!" << endl;
    exit( FAILURE );
}

RadiationSpectralModel * create_radiation_spectral_model( const string input_file )
{
    // 0. Create a Radiation_spectral_model pointer
    RadiationSpectralModel * rsm = 0;
    
    // 1. Get spectral_model string from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, input_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "set_radiation_spectral_model():\n";
	ost << "Error in input file: " << input_file << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "spectral_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "set_radiation_spectral_model():\n";
	ost << "Error in the declaration of spectral_data: a table is expected.\n";
	input_error(ost);
    }
    
    string spectral_model = get_string(L, -1, "spectral_model");
    
    if ( ECHO_RAD_INPUT > 0 ) 
	cout << "--- Creating a new " << spectral_model << " spectral model"
   	     << " from file: " << input_file << " ---" << endl;
    
    // 2. Create the spectral model
    if( spectral_model == "equilibrium_air" ) {
	rsm = new EquilibriumAir();
    }
    else if( spectral_model == "photaura" ) {
	rsm = new Photaura(L);
    }
    else if( spectral_model == "spradian" ) {
#	if WITH_SPRADIAN == 1
	rsm = new Spradian(L);
#       else
	cout << "Code not built with Spradian spectral radiation model available." << endl
	     << "Recompile using WITH_SPRADIAN==1." << endl
	     << "Exiting program." << endl;
	exit(BAD_INPUT_ERROR);
#	endif
    }
    else if( spectral_model == "parade" ) {
	rsm = new Parade(L);
    }
    else {
	cout << "The specified spectral radiation model: " << spectral_model << endl;
	cout << "is not available of no yet implemented.\n";
	cout << "Bailing Ou1!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    lua_pop(L,1);	// pop spectral_data
    lua_close(L);
    
    return rsm;
}

lua_State* initialise_radiation_lua_State()
{
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    return L;
}
