// Author: Daniel F. Potter
// Date: 18-Nov-2009

#include <iostream>
#include <sstream>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "energy-exchange-rate.hh"

using namespace std;

Energy_exchange_rate::
Energy_exchange_rate(lua_State *L)
{
    int imode = get_int(L, -1, "imode");

    // Now get mechanisms.
    lua_getfield(L, -1, "mechanisms");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Energy_exchange_rate::Energy_exchange_rate()\n";
	ost << "Error interpreting 'mechanisms'; a table of mechanisms is expected.\n";
	input_error(ost);
    }
    int nmechs = lua_objlen(L, -1);
    for ( int i = 1; i <= nmechs; ++i ) {
	lua_rawgeti(L, -1, i);
	ee_mech_.push_back(create_energy_exhange_mechanism(L, imode));
	lua_pop(L, 1);
    }
    lua_pop(L, 1);
}


Energy_exchange_rate::
~Energy_exchange_rate()
{
    for( size_t i = 0; i < ee_mech_.size(); ++i ) delete ee_mech_[i];
}

int
Energy_exchange_rate::
compute_all_relaxation_times(Gas_data &Q, vector<double> &molef)
{
    for( size_t i = 0; i < ee_mech_.size(); ++i ) {
	ee_mech_[i]->compute_relaxation_time(Q, molef);
    }

    return SUCCESS;
}

double
Energy_exchange_rate::
compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
{
    double rate = 0.0;
    for( size_t i = 0; i < ee_mech_.size(); ++i ) {
	rate += ee_mech_[i]->compute_rate(y, Q, molef);
    }
    return rate;
}

Energy_exchange_rate* create_energy_exchange_rate(lua_State *L)
{
    // FIXME - maybe sub-types of Energy_exchange_rate are needed
    return new Energy_exchange_rate(L);
}
