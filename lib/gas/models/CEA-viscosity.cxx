// Author: Rowan J. Gollan
// Date: 05-Nov-2008
// Place: Hampton, Virginia, USA

#include "../../util/source/lua_service.hh"

#include "CEA-viscosity.hh"

CEA_viscosity::
CEA_viscosity(lua_State *L)
    : Viscosity_model()
{    // Assume a table with the model paramaters is TOS.
    int ncurves = lua_objlen(L, -1);
    curve_.resize(ncurves);

    for ( int i = 0; i < ncurves; ++i ) {
	lua_rawgeti(L, -1, i+1); // Lua indexes from 1 not 0
	
	curve_[i].T_low = get_positive_number(L, -1, "T_low");
	curve_[i].T_high = get_positive_number(L, -1, "T_high");
	curve_[i].A = get_number(L, -1, "A");
	curve_[i].B = get_number(L, -1, "B");
	curve_[i].C = get_number(L, -1, "C");
	curve_[i].D = get_number(L, -1, "D");

	lua_pop(L, 1); // pop table of parameters from stack
    }
}

CEA_viscosity::
~CEA_viscosity() {}

double
CEA_viscosity::
s_eval_viscosity(const Gas_data &Q)
{
    double T = Q.T[0];
    double mu = 0.0;

    if ( T < curve_[0].T_low ) {
	// Just evaluate at curve_[0].T_low
	mu = eval_CEA_transport_curve(curve_[0].T_low, curve_[0]);
	return mu*1.0e-7; // convert CEA value in microPoise
	                  // to SI units: kg/(m.s)
    }
    
    if ( T > curve_.back().T_high ) {
	// Just evaluate at final T_high
	mu = eval_CEA_transport_curve(curve_.back().T_high, curve_.back());
	return mu*1.0e-7; // convert CEA value in microPoise
	                  // to SI units: kg/(m.s)
    }

    for ( size_t i = 0; i < curve_.size(); ++i ) {
	if ( check_T_transport_range(T, curve_[i]) ) {
	    mu = eval_CEA_transport_curve(T, curve_[i]);
	    return mu*1.0e-7; // convert CEA value in microPoise
                              // to SI units: kg/(m.s)
	}
    }

    // We should never get here.
    return 0.0;
}

