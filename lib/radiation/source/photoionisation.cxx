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

PhotoIonisationCrossSectionModel::PhotoIonisationCrossSectionModel( string name, double E_min )
 : name( name ), E_min( E_min ) {}
 
PhotoIonisationCrossSectionModel::~PhotoIonisationCrossSectionModel() {}

NoPICSModel::NoPICSModel()
 : PhotoIonisationCrossSectionModel( "NoPICSModel" ) {}
 
NoPICSModel::~NoPICSModel() {}

void NoPICSModel::spectral_distribution( vector<double> &nus )
{
    return;
}

double NoPICSModel::eval( double nu )
{
    return 0.0;
}

HydrogenicModel::HydrogenicModel( lua_State * L, double n_eff, int Z, double I, double E_min )
 : PhotoIonisationCrossSectionModel( "HydrogenicModel", E_min ),
   n_eff_( n_eff ), Z_( double(Z) ), I_( I )
{
    constB = I_ * Z_ * Z_;
    constC = n_eff_ * n_eff_;
    constD = pow( Z_, 4 );
    constE = pow( n_eff_, 5 );

    // FIXME: put these in the Lua file
    E_max = RC_k_SI * 1.0e5;
    nnus = 100;
}

HydrogenicModel::~HydrogenicModel() {}

void HydrogenicModel::spectral_distribution( vector<double> &nus )
{
    // 1. Determine some limits
    double nu_min = E_min / RC_h_SI;
    double nu_max = E_max / RC_h_SI;
    double dnu = ( nu_max - nu_min ) / double(nnus-1);

    // 2. Equidistant points in frequency space
    double nu = nu_min;
    while( nu <= nu_max ) {
        nus.push_back(nu);
        nu += dnu;
    }

    return;
}

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
    	cout << "pow( RC_h_SI * nu / ( I * Z * Z ), 1.0/3.0) = " 
	     << pow( RC_h_SI * nu / ( I_ * Z_ * Z_ ), 1.0/3.0) << endl;
    	cout << "( 2.0 / ( n_eff*n_eff ) * I * Z * Z / (RC_h_SI * nu) - 1.0 ) = " 
	     << ( 2.0 / ( n_eff_ * n_eff_ ) * I_ * Z_ * Z_ / (RC_h_SI * nu) - 1.0 ) << endl;
    	cout << "n_eff = " << n_eff_ << endl;
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
    // Initialise E_min to a large value
    E_min = 9.9e9;

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
 	    // Determine if this step corresponds to the minimum energy
 	    if ( nu_a*RC_h_SI < E_min ) E_min = nu_a*RC_h_SI;
 	    if ( nu_b*RC_h_SI < E_min ) E_min = nu_b*RC_h_SI;
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

void JohnstonStepModel::spectral_distribution( vector<double> &nus )
{
    for ( size_t istep=0; istep<steps.size(); ++istep ) {
        nus.push_back( steps[istep]->nu_a );
        nus.push_back( steps[istep]->nu_b );
    }

    return;
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

JohnstonThresholdModel::JohnstonThresholdModel( lua_State * L, int ilev, double E_min )
 : PhotoIonisationCrossSectionModel( "JohnstonThresholdModel", E_min )
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

    // FIXME: put these in the Lua file
    E_max = RC_k_SI * 1.0e5;
    nnus = 100;
}

JohnstonThresholdModel::~JohnstonThresholdModel() {}

void JohnstonThresholdModel::spectral_distribution( vector<double> &nus )
{
    // 1. Determine some limits
    double nu_min = E_min / RC_h_SI;
    double nu_max = E_max / RC_h_SI;
    double dnu = ( nu_max - nu_min ) / double(nnus-1);

    // 2. Equidistant points in frequency space
    double nu = nu_min;
    while( nu <= nu_max ) {
	nus.push_back(nu);
	nu += dnu;
    }

    return;
}

double JohnstonThresholdModel::eval( double nu )
{
    return sigma_bf_t * pow( nu_t / nu, theta );
}

TOPBaseModel::TOPBaseModel( lua_State * L, int ilev )
 : PhotoIonisationCrossSectionModel( "TOPBaseModel" )
{
    // Initialise E_min to a large value - it will be found below
    E_min = 9.9e9;

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
	
 	// Initialise the PICS params
        double nu = point_data[0]*RC_Ry_SI*RC_c_SI;	  // convert Ry -> Hz
        double sigma_bf = point_data[1] * 1.0e-18;	  // convert cm**2 x 1.0e18 -> cm**2
        nu_list.push_back( nu );
        sigma_list.push_back( sigma_bf );
        // Check if this point corresponds to E_min
        if ( RC_h_SI*nu < E_min ) E_min = RC_h_SI*nu;
    }
    
    // initialise i_prev to 0
    i_prev = 0;
    
    // initialise N_points
    N_points = (int) nu_list.size();
}

TOPBaseModel::~TOPBaseModel() {}

void TOPBaseModel::spectral_distribution( vector<double> &nus )
{
    for ( int i=0; i<N_points; ++i)
        nus.push_back( nu_list[i] );

    return;
}

double TOPBaseModel::eval( double nu )
{
#   if 0
    cout << "TOPBaseModel::eval( " << nu << " ) " << endl;
    cout << "N_points = " << N_points << endl;
    cout << "nu.front() = " << nu_list.front() << ", nu.back() = " << nu_list.back() << endl;
    cout << "i_prev = " << i_prev << endl;
#   endif
    
    // need to find bounding data points
    if ( nu < nu_list.front() ) return 0.0;
    else if ( nu > nu_list.back() ) return sigma_list.back();
    else if ( nu < nu_list[i_prev] || nu > nu_list[i_prev+1] ) {
        // Use bisection method
	int i_left = 0;
	int i_right = N_points-1;

	int i_mid=(i_left+i_right)/2;

	for(i_mid=(i_left+i_right)/2; abs(i_left-i_right) > 1; i_mid=(i_left+i_right)/2) {
	    double f_mid = nu_list[i_mid] - nu;
	    if (f_mid <= 0.0) {
	        i_left = i_mid; // use right interval
	    } else {
	        i_right = i_mid; // use left interval
	    }
	}
        i_prev = i_mid;
    }

    // sanity check
    if ( nu < nu_list[i_prev] || nu > nu_list[i_prev+1] ) {
	cout << "nu is out of range " << endl;
	exit(0);
    }

    // Linear interpolation between the found points
    double sigma = ( sigma_list[i_prev+1] - sigma_list[i_prev] ) /
		           ( nu_list[i_prev+1] - nu_list[i_prev] ) * ( nu - nu_list[i_prev]) +
		           sigma_list[i_prev];

    // cout << "nu = " << nu << ", nu_list[i_prev] = " << nu_list[i_prev] << ", nu_list[i_prev+1] = " << nu_list[i_prev+1] << ", sigma = " << sigma << endl;

    return sigma;
}

PhotoIonisationCrossSectionModel*
create_Hydrogenic_PICS_model( lua_State * L, int Z, double I, double E )
{
    PhotoIonisationCrossSectionModel * PICS_model;

    // hydrogenic effective principal quantum number
    double E_min = I - E;
    double n_eff = sqrt( RC_H_ionise_J / E_min );
    if ( !isfinite(n_eff) )
	PICS_model = new NoPICSModel();
    else
        PICS_model = new HydrogenicModel( L, n_eff, Z,  I, E_min );

    return PICS_model;
}

PhotoIonisationCrossSectionModel*
create_Johnston_PICS_model( lua_State * L, int ilev, double I, double E )
{
    PhotoIonisationCrossSectionModel * PICS_model;

    int nstep_levs = get_int(L,-1,"nstep_levs");
    if ( ilev < nstep_levs ) {
        PICS_model = new JohnstonStepModel( L, ilev );
    }
    else {
        double E_min = I - E;
        PICS_model = new JohnstonThresholdModel(L, ilev, E_min );
    }

    return PICS_model;
}

PhotoIonisationCrossSectionModel*
create_TOPBase_PICS_model( lua_State * L, int ilev, int Z, double I, double E )
{
    PhotoIonisationCrossSectionModel * PICS_model;

    bool hydrogenic_fill = get_boolean( L, -1, "hydrogenic_fill" );
    ostringstream level_label;
    level_label << "ilev_" << ilev;
    lua_getfield(L,-1,level_label.str().c_str());
    if ( !lua_istable(L, -1) ) {
        if ( hydrogenic_fill )
            PICS_model = create_Hydrogenic_PICS_model( L, Z, I, E );
        else
            PICS_model = new NoPICSModel();
    }
    else
        PICS_model = new TOPBaseModel( L, ilev );
    lua_pop(L,1);

    return PICS_model;
}

PhotoIonisationCrossSectionModel*
create_new_PICS_model( lua_State * L, int ilev, int Z, double I, double E )
{
    PhotoIonisationCrossSectionModel * PICS_model;
    
    string PICS_model_type = get_string( L, -1, "model" );
    
    if ( PICS_model_type=="none" ) {
    	PICS_model = new NoPICSModel();
    }
    else if ( PICS_model_type=="hydrogenic" ) {
	PICS_model = create_Hydrogenic_PICS_model( L, Z, I, E );
    }
    else if ( PICS_model_type=="JohnstonModel" ) {
    	PICS_model = create_Johnston_PICS_model( L, ilev, I, E );
    }
    else if ( PICS_model_type=="TOPBaseModel" ) {
        PICS_model = create_TOPBase_PICS_model( L, ilev, Z, I, E );
    }
    else {
    	cout << "create_new_PICS_model()" << endl
    	     << "Model with name '" << PICS_model_type << "' not understood." << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    return PICS_model;
}
