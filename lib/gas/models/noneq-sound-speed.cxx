// Author: Daniel F Potter
// Version: 
//   21-Sep-2009
//      Initial coding.
//   23-Sep-2009
//       Some refactoring to keep style consistent
//       with rest of module. (RJG)

#include <sstream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "noneq-sound-speed.hh"

using namespace std;

Noneq_sound_speed_model::
Noneq_sound_speed_model(Gas_model &gm, lua_State *L)
 : gm_( gm )
{
    // Search thermal modes for heavy-particle and free-electron translational
    // temperature indices
    lua_getglobal(L, "thermal_modes");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Noneq_sound_speed_model::Noneq_sound_speed_model():\n";
	ost << "Error in the declaration of thermal_modes: a table is expected.\n";
	input_error(ost);
    }
    
    iT_ = -1; iTe_ = -1;
    
    int ntm = lua_objlen(L, -1);
    
    for ( int itm = 0; itm < ntm; ++itm ) {
	lua_rawgeti(L, -1, itm+1);
	const char* mode = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring specific thermal mode table to TOS
	lua_getglobal(L, mode);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Noneq_sound_speed_model::Noneq_sound_speed_model():\n";
	    ost << "Error locating information table for mode: " << mode << endl;
	    input_error(ost);
	}
	
	lua_getfield( L, -1, "components" );
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Noneq_sound_speed_model::Noneq_sound_speed_model()\n";
	    ost << "Error locating 'components' table" << endl;
	    input_error(ost);
	}
	
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    string component_name = luaL_checkstring(L, -1);
	    lua_pop(L,1);
	    if ( component_name=="e_minus-translation" || 
	    	 component_name=="all-translation" ) iTe_ = itm;
	    if ( component_name=="hp-translation" || 
	    	 component_name=="all-translation" ) iT_ = itm;
	}
	
	lua_pop(L, 1);	// pop 'components'
	
	lua_pop(L, 1);	// pop 'mode'
    }
    
    lua_pop(L, 1);	// pop 'thermal_modes'
    
    if ( iT_ < 0 || iTe_ < 0 ) {
	ostringstream ost;
	ost << "Noneq_sound_speed_model::Noneq_sound_speed_model():\n";
	ost << "Heavy-particle and/or free-electron translational temperature\n"
	    << "could not be determined (iT=" << iT_ << ", iTe = " << iTe_ << ")\n";
	input_error(ost);
    }
}

Noneq_sound_speed_model::
~Noneq_sound_speed_model() {}

int
Noneq_sound_speed_model::
s_eval_sound_speed(Gas_data &Q)
{
    // Reference:
    // Cinnella and Grossman (1990)
    
    // "frozen" sound speed, corrected for presence of free electrons
    
    int status;
    double gamma = gm_.gamma(Q, status);
    if ( status != SUCCESS ) return status;
    
    Q.a = sqrt(gamma*gm_.dpdrho_const_T(Q, status) + (gamma - 1.0) * \
    	(Q.T[iT_] / Q.T[iTe_] - 1.0) * Q.p_e/Q.rho);
    return status;
}

