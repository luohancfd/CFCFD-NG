/** \file parade_radiator.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Jan-2012
 *  \brief Definitions for class describing a parade radiator
 *
 **/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>

#include "../../util/source/lua_service.hh"

#include "parade_radiator.hh"
#include "spectral_model.hh"

using namespace std;

/************************** Radiator ***************************/

ParadeRadiator::ParadeRadiator( lua_State * L, string name )
: name( name )
{
    // Basic radiator data
    type = get_string( L, -1, "type" );
    
    m_w = get_positive_number( L, -1, "mol_weight" );

    isp = get_int( L, -1, "isp" );
    
    if ( type=="diatomic_radiator" || type=="triatomic_radiator" ) {
	iTr = get_int( L, -1, "iTr" );
    
	iTv = get_int( L, -1, "iTv" );
    }
    else {
	iTr = 0;

	iTv = 0;
    }

    if ( ECHO_RAD_INPUT > 1 ) {
	cout << "isp = " << isp << endl
	     << "iTr = " << iTr << endl
	     << "iTv = " << iTv << endl;
    }
}

ParadeRadiator::
~ParadeRadiator() {}

ParadeRadiator * create_new_parade_radiator( lua_State * L, const std::string name )
{
    lua_getglobal(L, name.c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "create_new_parade_radiator()\n";
	ost << "Error locating information table for radiator: " << name << endl;
	input_error(ost);
    }

    cout << "Creating new parade radiator: " << name << endl;

    ParadeRadiator * new_parade_radiator = new ParadeRadiator( L, name );
    
    lua_pop(L,1); 	// pop radiator
    
    return new_parade_radiator;
}
