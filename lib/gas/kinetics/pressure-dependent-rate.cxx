// Author: Rowan J. Gollan
// Date: 16-Apr-2009
// Place: NIA, Hampton, Virginia, USA
//
// This is a port of code written by Brendan O'Flaherty
// which is found in gas_models2.
//

#include <cmath>
#include <sstream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "pressure-dependent-rate.hh"

using namespace std;

Pressure_dependent::
Pressure_dependent(lua_State *L, Gas_model &g, double T_upper, double T_lower)
    : Reaction_rate_coefficient(T_upper, T_lower)
{
    // Initialise k_inf
    lua_getfield(L, -1, "k_inf");
    k_inf_ = new Generalised_Arrhenius(L, g, T_upper, T_lower);
    lua_pop(L, 1);
    
    // Initialise k_0
    lua_getfield(L, -1, "k_0");
    k_0_ = new Generalised_Arrhenius(L, g, T_upper, T_lower);
    lua_pop(L, 1);

    // Pull out efficiencies
    read_table_as_map(L, -1, "efficiencies", efficiencies_);

    // Fill in Troe values, if available...
    Troe_model_ = false;
    lua_getfield(L, -1, "Troe");
    if ( lua_istable(L, -1) ) {
	Troe_model_ = true;
	lua_getfield(L, -1, "F_cent");
	if ( !lua_isnil(L, -1) ) {
	    F_cent_supplied_ = true;
	    F_cent_ = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	}
	else {
	    F_cent_supplied_ = false;
	    lua_pop(L, 1);
	    // We need to get other values.
	    lua_getfield(L, -1, "a"); a_ = luaL_checknumber(L, -1); lua_pop(L, 1);
	    lua_getfield(L, -1, "T1"); T1_ = luaL_checknumber(L, -1); lua_pop(L, 1);
	    lua_getfield(L, -1, "T3"); T3_ = luaL_checknumber(L, -1); lua_pop(L, 1);
	    lua_getfield(L, -1, "T2");
	    if ( !lua_isnumber(L, -1) ) {
		T2_supplied_ = false;
		T2_ = 0.0;
	    }
	    else {
		T2_supplied_ = true;
		T2_ = luaL_checknumber(L, -1);
	    }
	}
    }
    lua_pop(L, 1);

    // Setup array for storage of molecular weights
    M_.resize(g.get_number_of_species());

    // Fill in molecular weights...
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	M_[isp] = g.molecular_weight(isp);
    }

}

Pressure_dependent::
~Pressure_dependent()
{
    delete k_inf_;
    delete k_0_;
}

int
Pressure_dependent::
s_eval(const Gas_data &Q)
{
    // Find value of third body concentration
    double M = compute_third_body_concentration(Q);
    
    // Evaluate the limiting reaction rates (at high and low
    // pressure limits)
    k_inf_->eval(Q);
    k_0_->eval(Q);

    double k_inf = k_inf_->k();
    double k_0 = k_0_->k();

    double p_r = k_0*M/k_inf;
    double small = 1.0e-30;
    double log_p_r = log10(max(p_r, small));

    // Lindemann-Hinshelwood model
    double F = 1.0;

    // Check on temperature limits
    double T = Q.T[0];
    if ( T > T_upper_ )
	T = T_upper_;
    if ( T < T_lower_ )
	T = T_lower_;

    // Troe model
    if ( Troe_model_ ) {
	if ( !F_cent_supplied_ ) {
	    F_cent_ = (1.0 - a_)*exp(-T/T3_) + a_*exp(-T/T1_);
	    if ( T2_supplied_ ) {
		F_cent_ += exp(-T2_/T);
	    }
	}

	double log_F_cent = log10(max(F_cent_, small));
	double c = -0.4 - 0.67*log_F_cent;
	double n = 0.75 - 1.27*log_F_cent;
	double d = 0.14;

	double numer = log_p_r + c; 
	double denom = n - d*numer;
	double frac = numer/denom;
	double log_F = log_F_cent / (1.0 + frac*frac);
	F = pow(10,log_F);
    }

    k_ = F*k_0*k_inf*M/(k_0*M + k_inf);
    return SUCCESS;
}

// double
// Pressure_dependent::
// compute_third_body_concentration(const Gas_data &Q)
// {
//     // First compute concentration of requisite species...
//     double M = 0.0;
//     map<int, double>::const_iterator it;
//     for ( it = efficiencies_.begin(); it != efficiencies_.end(); ++it ) {
// 	int isp = it->first;
// 	double eff = it->second;
// 	M += eff * (Q.massf[isp]*Q.rho / M_[isp]);
//     }

//     return M;
// }

double
compute_third_body_value(const Gas_data &Q, map<int, double> efficiencies, vector<double> M)
{
    // First compute concentration of requisite species...
    double tbv = 0.0;
    map<int, double>::const_iterator it;
    for (it = efficiencies.begin(); it != efficiencies.end(); ++it) {
	int isp = it->first;
	double eff = it->second;
	tbv += eff * (Q.massf[isp]*Q.rho/M[isp]);
    }
    
    return tbv;
}

Reaction_rate_coefficient* create_pressure_dependent_coefficient(lua_State *L, Gas_model &g,
								 double T_upper, double T_lower)
{
    return new Pressure_dependent(L, g, T_upper, T_lower);
}
