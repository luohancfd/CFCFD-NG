// Author: Rowan J. Gollan
// Date: 30-July-2008


#include "../../util/source/lua_service.hh"

#include "CEA-thermal-conductivity.hh"

CEA_thermal_conductivity::
CEA_thermal_conductivity(lua_State *L)
    : Thermal_conductivity_model()
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

CEA_thermal_conductivity::
~CEA_thermal_conductivity() {}

double
CEA_thermal_conductivity::
s_eval_thermal_conductivity(const Gas_data &Q)
{
    double T = Q.T[0];
    double k = 0.0;

    if ( T < curve_[0].T_low ) {
	// Just evaluate at curve_[0].T_low
	k = eval_CEA_transport_curve(curve_[0].T_low, curve_[0]);
	return k*1.0e-4; // convert CEA value in microWatts/(cm.K)
	                 // to SI units: W/(m.K) 
    }
    
    if ( T > curve_.back().T_high ) {
	// Just evaluate at final T_high
	k = eval_CEA_transport_curve(curve_.back().T_high, curve_.back());
	return k*1.0e-4; // convert CEA value in microWatts/(cm.K)
	                 // to SI units: W/(m.K)
    }

    for ( size_t i = 0; i < curve_.size(); ++i ) {
	if ( check_T_transport_range(T, curve_[i]) ) {
	    k = eval_CEA_transport_curve(T, curve_[i]);
	    return k*1.0e-4; // convert CEA value in microWatts/(cm.K)
	                     // to SI units: W/(m.K)
	}
    }

    // We should never get here.
    return 0.0;
}

