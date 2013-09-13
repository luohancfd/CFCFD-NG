/** \file photoionisation.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 2-March-10: Initial implementation
 *
 *  \brief Definitions for photo-ionisation cross-section model classes
 *
 **/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "photoionisation.hh"
#include "radiation_constants.hh"

using namespace std;

PhotoIonisationCrossSectionModel::PhotoIonisationCrossSectionModel( string name )
 : name( name ) {}
 
PhotoIonisationCrossSectionModel::~PhotoIonisationCrossSectionModel() {}

NoPICSModel::NoPICSModel()
 : PhotoIonisationCrossSectionModel( "NoPICSModel" ) {}
 
NoPICSModel::~NoPICSModel() {}

double NoPICSModel::eval( double nu )
{
    return 0.0;
}

HydrogenicModel::HydrogenicModel( double n_eff, int Z, double I )
 : PhotoIonisationCrossSectionModel( "HydrogenicModel" ), 
   n_eff( n_eff ), Z( double(Z) ), I( I )
{
    constB = I * Z * Z;
    constC = n_eff*n_eff;
    constD = pow( Z, 4 );
    constE = pow( n_eff, 5 );
}

HydrogenicModel::~HydrogenicModel() {}

static const double constA = 2.8141909681753553e+29;

double HydrogenicModel::eval( double nu )
{
    // Ref 1: Johnston (2006 ) PhD Thesis p. 70
    // Ref 2: Zeldovich and Razier (1966) p. 261-272
    
    // See value defined above
    // double constA = 64.0 * pow( M_PI, 4 ) * RC_m * pow( RC_e, 10 ) / ( 3.0 * sqrt(3.0) * RC_c * RC_h**6 );
    
    // Gaunt factor
    // Ref: Zeldovich and Razier (1966) p. 266
    double E = RC_h_SI * nu;
    double G = 1.0 - 0.173 * pow( E / constB , 1.0/3.0) *
		   ( 2.0 / constC * constB / E - 1.0 );
    
    // Cross-section in cm**2
    // CHECKME: - are pow() functions inefficient here?
    double sigma_bf = constA * constD / ( ( nu*nu*nu ) * constE ) * G;
    
#if DEBUG_RAD > 0
    if ( isnan(sigma_bf) || isinf(sigma_bf) ) {
    	cout << "sigma_bf = " << sigma_bf << ", constA = " << constA << ", G = " << G << endl;
    	cout << "pow( RC_h_SI * nu / ( I * Z * Z ), 1.0/3.0) = " << pow( RC_h_SI * nu / ( I * Z * Z ), 1.0/3.0) << endl;
    	cout << "( 2.0 / ( n_eff*n_eff ) * I * Z * Z / (RC_h_SI * nu) - 1.0 ) = " << ( 2.0 / ( n_eff*n_eff ) * I * Z * Z / (RC_h_SI * nu) - 1.0 ) << endl;
    	cout << "n_eff = " << n_eff << endl;
    	exit(FAILURE);
    }
#endif
    
    return sigma_bf;
}

PICS_step::PICS_step( double nu_a, double nu_b, double sigma_bf )
 : nu_a( nu_a ), nu_b( nu_b ), sigma_bf( sigma_bf )
{}

PICS_step::~PICS_step() {}

JohnstonStepModel::JohnstonStepModel( lua_State * L, int ilev )
 : PhotoIonisationCrossSectionModel( "JohnstonStepModel" )
{
    int nsteps = get_int( L, -1, "nsteps" );
    
    for ( int istep=0; istep<nsteps; ++istep ) {
    	ostringstream step_label;
    	step_label << "step_" << istep;
    	lua_getfield(L,-1,step_label.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "JohnstonStepModel::JohnstonStepModel()\n";
	    ost << "Error locating " << step_label.str() << " table" << endl;
	    input_error(ost);
	}
	vector<double> step_data;
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    step_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop istep
	// Check the size of the step data vector
	if ( step_data.size()!=4 ) {
	    ostringstream oss;
	    oss << "JohnstonStepModel::JohnstonStepModel()" << endl
	        << "Step data expected to have 4 elements." << endl;
	    input_error( oss );
	}	    
 	// Create the step if it corresponds to this ilev
 	if ( ilev==int(step_data[0]) ) {
 	    double nu_a = step_data[1]*RC_e_SI/RC_h_SI;     // convert eV -> Hz
 	    double nu_b = step_data[2]*RC_e_SI/RC_h_SI;	    // convert eV -> Hz
 	    double sigma_bf = step_data[3] * 1.0e-18;	    // convert cm**2 x 1.0e18 -> cm**2 
 	    steps.push_back( new PICS_step( nu_a, nu_b, sigma_bf ) );
 	}
    }
    
    cout << "JohnstonStepModel::JohnstonStepModel()" << endl
         << "- Created " << steps.size() << " steps for level " << ilev << endl;
}

JohnstonStepModel::~JohnstonStepModel()
{
    for ( size_t istep=0; istep<steps.size(); ++istep )
    	delete steps[istep];
}

double JohnstonStepModel::eval( double nu )
{
    for ( size_t istep=0; istep<steps.size(); ++istep ) {
    	if ( nu>=steps[istep]->nu_a && nu<=steps[istep]->nu_b )
    	    return steps[istep]->sigma_bf;
    	// else continue;
    }
    
    // if we get here the cross-section is zero
    return 0.0;
}

JohnstonThresholdModel::JohnstonThresholdModel( lua_State * L, int ilev )
 : PhotoIonisationCrossSectionModel( "JohnstonThresholdModel" )
{
    int nthresholds = get_int( L, -1, "nthresholds" );
    
    for ( int ithreshold=0; ithreshold<nthresholds; ++ithreshold ) {
    	ostringstream threshold_label;
    	threshold_label << "threshold_" << ithreshold;
    	lua_getfield(L,-1,threshold_label.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "JohnstonThresholdModel::JohnstonThresholdModel()\n";
	    ost << "Error locating " << threshold_label.str() << " table" << endl;
	    input_error(ost);
	}
	vector<double> threshold_data;
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    threshold_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop ithreshold
	
	// Check the size of the threshold data vector
	if ( threshold_data.size()!=4 ) {
	    ostringstream oss;
	    oss << "JohnstonThresholdModel::JohnstonThresholdModel()" << endl
	        << "Threshold data expected to have 4 elements." << endl;
	    input_error( oss );
	}
	
 	// Initialise the threshold params if they corresponds to this ilev
 	if ( ilev==int(threshold_data[0]) ) {
 	    nu_t = threshold_data[1]*RC_e_SI/RC_h_SI;     // convert eV -> Hz
 	    sigma_bf_t = threshold_data[2] * 1.0e-18;	  // convert cm**2 x 1.0e18 -> cm**2
 	    theta = threshold_data[3];
 	    // break the threshold loop as there is only one per level
 	    break;
 	}
    }
}

JohnstonThresholdModel::~JohnstonThresholdModel() {}

double JohnstonThresholdModel::eval( double nu )
{
    return sigma_bf_t * pow( nu_t / nu, theta );
}

TOPBaseModel::TOPBaseModel( lua_State * L, int ilev )
 : PhotoIonisationCrossSectionModel( "TOPBaseModel" )
{
    int npoints = get_int( L, -1, "npoints" );
    
    for ( int i=0; i<npoints; ++i ) {
    	ostringstream datapoint_label;
    	datapoint_label << "point_" << i;
    	lua_getfield(L,-1,datapoint_label.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "TOPBaseModel::TOPBaseModel()\n";
	    ost << "Error locating " << datapoint_label.str() << " table" << endl;
	    input_error(ost);
	}
	vector<double> point_data;
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    point_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop ithreshold
	
	// Check the size of the threshold data vector
	if ( point_data.size()!=2 ) {
	    ostringstream oss;
	    oss << "TOPBaseModel::TOPBaseModel()" << endl
	        << "Cross-section data expected to have 2 elements." << endl;
	    input_error( oss );
	}
	
 	// Initialise the threshold params if they corresponds to this ilev
 	if ( ilev==int(point_data[0]) ) {
 	    double nu = point_data[1]*RC_Ry_SI*RC_c_SI;	  // convert Ry -> Hz
 	    double sigma_bf = point_data[2] * 1.0e-18;	  // convert cm**2 x 1.0e18 -> cm**2
 	    nu_list.push_back( nu );
 	    sigma_list.push_back( sigma_bf );
 	}
    }
    
    // initialise i_prev to 0
    i_prev = 0;
    
    // initialise N_points
    N_points = (int) nu_list.size();
}

TOPBaseModel::~TOPBaseModel() {}

double TOPBaseModel::eval( double nu )
{
    cout << "TOPBaseModel::eval( " << nu << " ) " << endl;
    cout << "nu.front() = " << nu_list.front() << ", nu.back() = " << nu_list.back() << endl;
    cout << "i_prev = " << i_prev << endl;
    
    // need to find bounding data points
    if ( nu < nu_list.front() ) return 0.0;
    else if ( nu > nu_list.back() ) return 0.0;
    else {
    	// start from i_prev
        for ( int i=i_prev; i<(N_points-1); ++i ) {
            cout << "i = " << i << endl;
    	    if ( nu >= nu_list[i] && nu < nu_list[i+1] )
    	        return 0.5 * ( sigma_list[i] + sigma_list[i+1] );
    	}
	// start from 0
    	for ( int i=0; i<(i_prev-1); ++i ) {
    	    cout << "i = " << i << endl;
    	    if ( nu >= nu_list[i] && nu < nu_list[i] ) {
    	        return 0.5 * ( sigma_list[i] + sigma_list[i+1] );
    	    }
        }
    }
    
    // if we get to here then the search failed
    cout << "TOPBaseModel::eval()" << endl
         << "Failed to find cross-section data for nu = " << nu << endl;
    exit( FAILURE );
}


PhotoIonisationCrossSectionModel*
create_new_PICS_model( lua_State * L, int ilev, double n_eff, int Z, double I )
{
    PhotoIonisationCrossSectionModel * PICS_model;
    
    string PICS_model_type = get_string( L, -1, "model" );
    
    if ( PICS_model_type=="none" ) {
    	PICS_model = new NoPICSModel();
    }
    else if ( PICS_model_type=="hydrogenic" ) {
    	PICS_model = new HydrogenicModel( n_eff, Z,  I );
    }
    else if ( PICS_model_type=="JohnstonModel" ) {
    	int nstep_levs = get_int(L,-1,"nstep_levs");
    	if ( ilev < nstep_levs ) {
    	    PICS_model = new JohnstonStepModel(L, ilev);
    	}
    	else {
    	    PICS_model = new JohnstonThresholdModel(L, ilev);
    	}
    }
    else if ( PICS_model_type=="TOPBaseModel" ) {
        ostringstream level_label;
        level_label << "ilev_" << ilev;
        lua_getfield(L,-1,level_label.str().c_str());
        if ( !lua_istable(L, -1) )
            PICS_model = new NoPICSModel();
        else
            PICS_model = new TOPBaseModel( L, ilev );
        lua_pop(L,1);
    }
    else {
    	cout << "create_new_PICS_model()" << endl
    	     << "Model with name '" << PICS_model_type << "' not understood." << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    return PICS_model;
}
