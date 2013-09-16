/** \file poshax_radiation_transport.cxx
 *  \ingroup poshax2
 *  \brief Definitions for the poshax radiation transport model class and functions.
 **/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdlib>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../../lib/radiation/source/spectral_model.hh"
#include "../../../lib/util/source/useful.h"

#include "poshax_radiation_transport.hh"

#define VERBOSE_RADIATION_TRANSPORT 0

using namespace std;

/* --------- Model: "PoshaxRadiationTransportModel: ----------- */
PoshaxRadiationTransportModel::PoshaxRadiationTransportModel( lua_State *L )
{
    spectrally_resolved_ = (int) get_boolean(L,-1,"spectrally_resolved");
    electronic_mode_factor_ = get_number(L,-1,"electronic_mode_factor");
}

PoshaxRadiationTransportModel::~PoshaxRadiationTransportModel()
{
    delete rsm_;
}

int PoshaxRadiationTransportModel::set_radiation_spectral_model( string file_name )
{
    // Create one spectral model
    rsm_ = create_radiation_spectral_model( file_name );
    
    return SUCCESS;
}

/* --------- Model: "OpticallyThin" --------- */

OpticallyThin::OpticallyThin( lua_State *L )
: PoshaxRadiationTransportModel(L) {}

OpticallyThin::~OpticallyThin() {}

string
OpticallyThin::str() const
{
    return "OpticallyThin";
}

double
OpticallyThin::eval_Q_rad( Gas_data &Q )
{
    // No flowfield reabsorption, 100% emission
    return ( - 4.0 * M_PI ) * \
    rsm_->radiative_integrated_emission_for_gas_state(Q, spectrally_resolved_);
}

/* --------- Model: "OpticallyThick" --------- */

OpticallyThick::OpticallyThick( lua_State *L )
: PoshaxRadiationTransportModel(L) {}

OpticallyThick::~OpticallyThick() {}

string
OpticallyThick::str() const
{
    return "OpticallyThick";
}

double
OpticallyThick::eval_Q_rad( Gas_data &Q )
{
    // 100% flowfield reabsorption
    return 0.0;
}

/* --------- Model: "OpticallyVariable" --------- */

OpticallyVariable::OpticallyVariable( lua_State *L )
: PoshaxRadiationTransportModel(L)
{
    wavel_switch_ = get_number( L, -1, "optical_switch" );
    Lambda_lower_ = get_number( L, -1, "lower_escape_factor" );
    Lambda_upper_ = get_number( L, -1, "upper_escape_factor" );
    
    cout << "wavel_switch = " << wavel_switch_ << endl;
    cout << "Lambda_lower = " << Lambda_lower_ << endl;
    cout << "Lambda_upper = " << Lambda_upper_ << endl;
}

OpticallyVariable::~OpticallyVariable() {}

string
OpticallyVariable::str() const
{
    return "OpticallyVariable";
}

double
OpticallyVariable::eval_Q_rad( Gas_data &Q )
{
    // No flowfield reabsorption, 100% emission
    return ( - 4.0 * M_PI ) * \
    rsm_->radiative_variably_integrated_emission_for_gas_state(Q, wavel_switch_, 
    	                                                       Lambda_lower_,
    	                                                       Lambda_upper_,
    	                                                       spectrally_resolved_);
}

PoshaxRadiationTransportModel
* create_poshax_radiation_transport_model( const string file_name )
{
    // 0. Initialise a Radiation_transport_model pointer
    PoshaxRadiationTransportModel * rtm;
    
    // 1. Get transport model name from lua file
    lua_State *L = initialise_radiation_lua_State();

    if( luaL_dofile(L, file_name.c_str()) != 0 ) {
	ostringstream ost;
	ost << "set_poshax_radiation_transport_model():\n";
	ost << "Error in input file: " << file_name << endl;
	input_error(ost);
    }
    
    lua_getglobal(L,"transport_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "set_radiation_transport_model():\n";
	ost << "Error in the declaration of transport_data - a table is expected.\n";
	input_error(ost);
    }
    
    string transport_model = get_string(L, -1, "transport_model");
    
    // 2. Create the transport model
    if ( ECHO_RAD_INPUT > 0 ) 
	cout << "--- Creating a new " << transport_model << " transport model"
   	     << " from file: " << file_name << " ---" << endl;
        
    if( transport_model == "optically thin" ) {
	rtm = new OpticallyThin(L);
    }
    else if( transport_model == "optically thick" ) {
	rtm = new OpticallyThick(L);
    }
    else if( transport_model == "optically variable" ) {
	rtm = new OpticallyVariable(L);
    }
    else {
	cout << "The specified radiation transport model: " << transport_model
	     << endl
	     << "is not available in poshax2 or no yet implemented.\n"
	     << "Bailing Out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    lua_pop(L,1);	// pop transport_data
    lua_close(L);
    
    // 3. Create the spectral model(s)
    rtm->set_radiation_spectral_model( file_name );
    
    return rtm;
}

