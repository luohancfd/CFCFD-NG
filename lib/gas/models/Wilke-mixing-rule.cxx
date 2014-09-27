// Author: Rowan J. Gollan
// Date: 10-Jul-2008

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "gas-model.hh"
#include "Wilke-mixing-rule.hh"
#include "Sutherland-viscosity.hh"
#include "CEA-viscosity.hh"
#include "Sutherland-thermal-conductivity.hh"
#include "CEA-thermal-conductivity.hh"
#include "physical_constants.hh"

using namespace std;

Wilke_mixing_rule::
Wilke_mixing_rule(lua_State *L)
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Wilke_mixing_rule::Wilke_mixing_rule():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, -1);
    M_.resize(nsp, 0.0);
    x_.resize(nsp, 0.0);
    mu_.resize(nsp, 0.0);
    k_.resize(nsp, 0.0);

    phi_.resize(nsp);
    for( int i = 0; i < nsp; ++i ) phi_[i].resize(nsp, 0.0);
    psi_.resize(nsp);
    for( int i = 0; i < nsp; ++i ) psi_[i].resize(nsp, 0.0);

    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Wilke_mixing_rule::Wilke_mixing_rule()\n";
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
	    ost << "Wilke_mixing_rule::Wilke_mixing_rule()\n";
	    ost << "Error setting viscosity table for species: " << sp << endl;
	    input_error(ost);
	}

	string vmodel = get_string(L, -1, "model");

	//cout << "vmodel= " << vmodel << endl;
	
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
	else {
	    ostringstream ost;
	    ost << "The viscosity model set for species: " << sp << endl;
	    ost << vmodel << " is not known.\n";
	    input_error(ost);
	}
	lua_pop(L, 1); // Pop 'viscosity' table off stack


	//cout << "Initialising therm conductivity model...\n";
	// Initialise a thermal conductivity model for component 'isp'
	lua_getfield(L, -1, "thermal_conductivity");
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Wilke_mixing_rule::Wilke_mixing_rule()\n";
	    ost << "Error setting thermal_conductivity table for species: " << sp << endl;
	    input_error(ost);
	}

	string kmodel = get_string(L, -1, "model");
	
	if ( kmodel == "Sutherland" ) {
	    lua_getfield(L, -1, "parameters");
	    TCM_.push_back(new Sutherland_thermal_conductivity(L));
	    lua_pop(L, 1);
	}
	else if ( kmodel == "CEA" ) {
	    lua_getfield(L, -1, "parameters");
	    TCM_.push_back(new CEA_thermal_conductivity(L));
	    lua_pop(L, 1);
	}
	else {
	    ostringstream ost;
	    ost << "The thermal conductivity model set for species: " << sp << endl;
	    ost << vmodel << " is not known.\n";
	    input_error(ost);
	}
	lua_pop(L, 1); // Pop 'thermal conductivity' table off stack

	lua_pop(L, 1); // Pops "sp" off stack
    }
    lua_pop(L, 1); // Pops "species" off stack

    lua_getglobal(L, "ignore_mole_fraction");
    if ( lua_isnil(L, -1) ) {
	cout << "In Wilke's mixing rule, setting ignore_mole_fraction\n";
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
    
}

Wilke_mixing_rule::
~Wilke_mixing_rule()
{
    for ( size_t isp = 0; isp < VM_.size(); ++isp ) {
	delete VM_[isp];
    }

    for ( size_t isp = 0; isp < TCM_.size(); ++isp ) {
	delete TCM_[isp];
    }
}

int
Wilke_mixing_rule::
s_eval_transport_coefficients(Gas_data &Q, Gas_model *gmodel)
{
    // Reference:
    // Wilke, C.R. (1950)
    // A Viscosity Equation for Gas Mixtures
    // Journal of Chemical Physics, 18:4 pp. 517--519
    // Implementation of equations (12), (13) and (14)
    // p. 519
    //
    // Note: White states that this mixing rule may also be
    //       used for thermal conductivities (p. 34).
    //
    // White, F.M. (2006)
    // Viscous Fluid Flow, 3rd edition
    // McGraw Hill International, New York 
    //

    if ( VM_.size() == 1 ) {
	// Assume there is only one species.
	Q.mu = VM_[0]->eval_viscosity(Q);
	Q.k[0] = TCM_[0]->eval_thermal_conductivity(Q);
	return SUCCESS;
    }
    else {
	// Set up values for calculation before mixing
	convert_massf2molef(Q.massf, M_, x_);
	for ( size_t i = 0; i < Q.massf.size(); ++i ) {
	    mu_[i] = VM_[i]->eval_viscosity(Q);
	    k_[i] = TCM_[i]->eval_thermal_conductivity(Q);
	}

	// Calculate interaction potentials
	for( size_t i = 0; i < Q.massf.size(); ++i ) {
	    for ( size_t j = 0; j < Q.massf.size(); ++j ) {
		double numer = pow((1.0 + sqrt(mu_[i]/mu_[j])*pow(M_[j]/M_[i], 0.25)), 2.0);
		double denom = (4.0/sqrt(2.0))*sqrt(1.0 + (M_[i]/M_[j]));
		phi_[i][j] = numer/denom;
		numer = pow((1.0 + sqrt(k_[i]/k_[j])*pow(M_[j]/M_[i], 0.25)), 2.0);
		psi_[i][j] = numer/denom;
	    }
	}

	// Now apply mixing formula
	Q.mu = 0.0;
	Q.k[0] = 0.0;
	double sum_a = 0.0;
	double sum_b = 0.0;
	for ( size_t i = 0; i < Q.massf.size(); ++i ) {
	    if( x_[i] < ignore_mole_fraction_ ) continue;
	    sum_a = 0.0;
	    sum_b = 0.0;
	    for ( size_t j = 0; j < Q.massf.size(); ++j ) {
		if ( i == j ) continue;
		if( x_[i] < ignore_mole_fraction_ ) continue; 
		sum_a += x_[j]*phi_[i][j];
		sum_b += x_[j]*psi_[i][j];
	    }
	    Q.mu += mu_[i]/(1.0 + (1.0/x_[i])*sum_a);
	    Q.k[0] += k_[i]/(1.0 + (1.0/x_[i])*sum_b);
	}
    }
    return SUCCESS;
}

