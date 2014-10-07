// Author: Elise J. Fahy
// Date: 24-Sep-2014

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "gas-model.hh"
#include "Armaly-Sutton-mixing-rule.hh"
#include "Wilke-mixing-rule.hh"
#include "Sutherland-viscosity.hh"
#include "Blottner-viscosity.hh"
#include "CEA-viscosity.hh"
#include "Sutherland-thermal-conductivity.hh"
#include "CEA-thermal-conductivity.hh"
#include "constant-Prandtl-thermal-conductivity.hh"
#include "physical_constants.hh"

using namespace std;

Armaly_Sutton_mixing_rule::
Armaly_Sutton_mixing_rule(lua_State *L)
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Armaly_Sutton_mixing_rule::Armaly_Sutton_mixing_rule():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, -1);
    M_.resize(nsp, 0.0);
    x_.resize(nsp, 0.0);
    mu_.resize(nsp, 0.0);

    phi_.resize(nsp);
    A_.resize(nsp);
    B_.resize(nsp);
    F_.resize(nsp);
    for( int i = 0; i < nsp; ++i ) {
	phi_[i].resize(nsp, 0.0);
	A_[i].resize(nsp, 0.0);
	B_[i].resize(nsp, 0.0);
	F_[i].resize(nsp, 0.0);
    }

    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Armaly_Sutton_mixing_rule::Armaly_Sutton_mixing_rule()\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	// Get molecular weight
	double M = get_positive_value(L, -1, "M");
	M_[isp] = M;

	// Initialise a viscosity model for component 'isp'

	//cout << "initialising viscosity model...\n";
	lua_getfield(L, -1, "viscosity");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Armaly_Sutton_mixing_rule::Armaly_Sutton_mixing_rule()\n";
	    ost << "Error setting viscosity table for species: " << sp << endl;
	    input_error(ost);
	}

	string vmodel = get_string(L, -1, "model");

	if ( vmodel == "Sutherland" ) {
	    lua_getfield(L, -1, "parameters");
	    VM_.push_back(new Sutherland_viscosity(L));
	    lua_pop(L, 1);
	}
	else if ( vmodel == "CEA" ) {
	    lua_getfield(L, -1, "parameters");
	    VM_.push_back(new CEA_viscosity(L));
	    lua_pop(L, 1);
	}
	else if ( vmodel == "Blottner" ) {
	    lua_getfield(L, -1, "parameters");
	    VM_.push_back(new Blottner_viscosity(L));
	    lua_pop(L, 1);
	}
	else {
	    ostringstream ost;
	    ost << "The viscosity model set for species: " << sp << endl;
	    ost << vmodel << " is not known.\n";
	    input_error(ost);
	}
	lua_pop(L, 1); // Pop 'viscosity' table off stack

	lua_pop(L, 1); // Pops "sp" off stack
    }
    lua_pop(L, 1); // Pops "species" off stack

    lua_getglobal(L, "ignore_mole_fraction");
    if ( lua_isnil(L, -1) ) {
	cout << "In Armaly-Sutton's mixing rule, setting ignore_mole_fraction\n";
	cout << "to default value: " << DEFAULT_MIN_MOLE_FRACTION << endl;
	ignore_mole_fraction_ = DEFAULT_MIN_MOLE_FRACTION;
    }
    else if ( lua_isnumber(L, -1) ) {
	ignore_mole_fraction_ = lua_tonumber(L, -1);
	if ( ignore_mole_fraction_ <= 0.0 ) {
	    ostringstream ost;
	    ost << "The ignore_mole_fraction is set to a value less than or equal to zero: " << ignore_mole_fraction_ << endl;
	    ost << "This must be a value between 0 and 1.\n";
	    input_error(ost); 
	}
	if ( ignore_mole_fraction_ > 1.0 ) {
	    ostringstream ost;
	    ost << "The ignore_mole_fraction is set to a value greater than 1.0: " << ignore_mole_fraction_ << endl;
	    ost << "This must be a value between 0 and 1.\n";
	    input_error(ost);
	}
    }
    else {
	ostringstream ost;
	ost << "The type ignore_mole_fraction is not a numeric type.\n";
	ost << "This must be a value between 0 and 1.\n";
	input_error(ost);
    }
    lua_pop(L, 1); // Pops "ignore_mole_fraction" off stack.

    // Now populate A_, B_ and F_.
    lua_getglobal(L, "ArmalySutton_params");
    lua_getfield(L, -1, "A");
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1);
	for ( int jsp = 0; jsp < nsp; ++jsp ) {
	    lua_rawgeti(L, -1, jsp+1);
	    A_[isp][jsp] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);
    }
    lua_pop(L, 1);
    lua_getfield(L, -1, "B");
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1);
	for ( int jsp = 0; jsp < nsp; ++jsp ) {
	    lua_rawgeti(L, -1, jsp+1);
	    B_[isp][jsp] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);
    }
    lua_pop(L, 1);
    lua_getfield(L, -1, "F");
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1);
	for ( int jsp = 0; jsp < nsp; ++jsp ) {
	    lua_rawgeti(L, -1, jsp+1);
	    F_[isp][jsp] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);
    }
    lua_pop(L, 1);
    lua_getfield(L, -1, "Pr");
    Pr_ = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    // Pop Armaly-Sutton parameters off stack.
    lua_pop(L, 1);

}

Armaly_Sutton_mixing_rule::
~Armaly_Sutton_mixing_rule()
{
    for ( size_t isp = 0; isp < VM_.size(); ++isp ) {
	delete VM_[isp];
    }
}

int
Armaly_Sutton_mixing_rule::
s_eval_transport_coefficients(Gas_data &Q, Gas_model *gmodel)
{
    // Reference:
    // Palmer and Wright (2003)
    // Comparison of Methods to Compute High-Temperature Gas Viscosity
    // Journal of Thermophysics and Heat Transfer, 17:2, pp.232-239
    // Implementation of equations (29) and (35).


    // Presently, we read in A, B and F values.
    // A_ik = 1.25 for all interactions except an atom with its own ion
    //      = 1.1 for an atom with its own ion, except air species!
    //      = 0.21 for N-N+ and O-O+ interactions
    // B_ik = 0.78 for neutral-neutral interactions
    //      = 0.15 for neutral-ion interactions
    //      = 0.2 for neutral-electron interactions
    //      = 1.0 for ion-ion or ion-electron interactions
    // F_ik = 1.0 for all interactions.

    if ( VM_.size() == 1 ) {
	// Assume there is only one species.
	Q.mu = VM_[0]->eval_viscosity(Q);
	return SUCCESS;
    }
    else {
	// Set up values for calculation before mixing
	convert_massf2molef(Q.massf, M_, x_);
	for ( size_t i = 0; i < Q.massf.size(); ++i ) {
	    mu_[i] = VM_[i]->eval_viscosity(Q);
	}

	// Calculate interaction potentials
	for( size_t i = 0; i < Q.massf.size(); ++i ) {
	    for ( size_t j = 0; j < Q.massf.size(); ++j ) {
		double numer = ( 5.0/(3.0*A_[i][j]) + M_[j]/M_[i] ) * pow(( F_[i][j] + B_[i][j]*sqrt(mu_[i]/mu_[j])*pow(M_[j]/M_[i], 0.25) ), 2.0);
		double denom = ( (1.0 + M_[j]/M_[i]) * sqrt(8.0*(1.0 + M_[j]/M_[i])) );
		phi_[i][j] = numer/denom;
	    }
	}

	// Now apply mixing formula - same formulation as Wilke.
	Q.mu = 0.0;
	double sum_a = 0.0;
	for ( size_t i = 0; i < Q.massf.size(); ++i ) {
	    if( x_[i] < ignore_mole_fraction_ ) continue;
	    sum_a = 0.0;
	    for ( size_t j = 0; j < Q.massf.size(); ++j ) {
		if ( i == j ) continue;
		if( x_[i] < ignore_mole_fraction_ ) continue; 
		sum_a += x_[j]*phi_[i][j];
	    }
	    Q.mu += mu_[i]/(1.0 + (1.0/x_[i])*sum_a);
	}
    }

    double Cp = gmodel->Cp(Q);
    // ELISE: Check this expression please.
    Q.k[0] = Q.mu*Cp/Pr_;

    return SUCCESS;
}

