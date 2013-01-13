/** \file spradian_radiator.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Jan-2012
 *  \brief Definitions for class describing a spradian radiator
 *
 **/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>

#include "../../util/source/lua_service.hh"

#include "spradian_radiator.hh"
#include "spectral_model.hh"
#include "radiation_constants.hh"

using namespace std;

/************************** Radiator ***************************/

SpradianRadiator::SpradianRadiator( lua_State * L, string name )
: name( name )
{
    // Basic radiator data
    type = get_string( L, -1, "type" );
    
    m_w = get_positive_number( L, -1, "mol_weight" );

    isp = get_int( L, -1, "isp" );
    
    if ( ECHO_RAD_INPUT > 1 ) {
	cout << "isp = " << isp << endl;
    }
}

SpradianRadiator::
~SpradianRadiator() {}

void SpradianRadiator::set_concentration( double rho_i )
{
    *conc = rho_i / m_w * RC_Na;
}

SpradianRadiator * create_new_spradian_radiator( lua_State * L, const std::string name )
{
    lua_getglobal(L, name.c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "create_new_spradian_radiator()\n";
	ost << "Error locating information table for radiator: " << name << endl;
	input_error(ost);
    }

    cout << "Creating new spradian radiator: " << name << endl;

    SpradianRadiator * new_spradian_radiator = new SpradianRadiator( L, name );
    
    lua_pop(L,1); 	// pop radiator
    
    return new_spradian_radiator;
}
