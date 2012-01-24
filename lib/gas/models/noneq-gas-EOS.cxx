// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#include <sstream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "physical_constants.hh"
#include "noneq-gas-EOS.hh"

using namespace std;

Noneq_gas::
Noneq_gas(lua_State *L)
    : Equation_of_state()
{
    // Get required data (R_[] and iel) from species
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Noneq_gas::Noneq_gas():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    // Use this loop over the species to also determine the electron index
    iel_ = -1;
    
    int nsp = lua_objlen(L, -1);
    
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);
	
	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Noneq_gas::Noneq_gas()\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	// NOTE: M_ vector is not used internally
	double M = get_positive_value(L, -1, "M");
	M_.push_back(M);
	R_.push_back(PC_R_u/M);
	
	// Check for electron
	string type = get_string(L, -1, "species_type");
	if ( type=="free electron" ) iel_ = isp;

	lua_pop(L, 1); // pop "sp" off stack
    }
    lua_pop(L, 1); // pop "species" off stack
    
    if ( iel_ < 0 ) {
	ostringstream ost;
	ost << "Noneq_gas::Noneq_gas():\n";
	ost << "A 'free electron' species was not found.\n";
	input_error(ost);
    }
    
    // Now search thermal modes for heavy-particle and free-electron translational
    // temperature indices
    lua_getglobal(L, "thermal_modes");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Noneq_gas::Noneq_gas():\n";
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
	    ost << "Noneq_gas::Noneq_gas():\n";
	    ost << "Error locating information table for mode: " << mode << endl;
	    input_error(ost);
	}
	
	lua_getfield( L, -1, "components" );
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Noneq_gas::Noneq_gas()\n";
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
	ost << "Noneq_gas::Noneq_gas():\n";
	ost << "Heavy-particle and/or free-electron translational temperature\n"
	    << "could not be determined (iT=" << iT_ << ", iTe = " << iTe_ << ")\n";
	input_error(ost);
    }

}

Noneq_gas::
~Noneq_gas() {}

int
Noneq_gas::
s_eval_pressure(Gas_data &Q)
{
    // 1. heavy-particle contribution
    double R_hp = hp_gas_constant(Q);
    Q.p  = Q.rho * Q.T[iT_] * R_hp;
    
    // 2. Free-electron contribution
    double R_e = Q.massf[iel_] * R_[iel_];
    Q.p_e = Q.rho * Q.T[iTe_] * R_e;
    Q.p += Q.p_e;
    
    return SUCCESS;
}

int
Noneq_gas::
s_eval_temperature(Gas_data &Q)
{
    // 1. Heavy-particle temperature    
    double R_hp = hp_gas_constant(Q);
    Q.T[iT_] = ( Q.p - Q.p_e ) / ( Q.rho * R_hp );
    
    // 2. Free-electron temperature
    double R_e = Q.massf[iel_] * R_[iel_];
    Q.T[iTe_] = Q.p_e / ( Q.rho * R_e );
    
    return SUCCESS;
}

int
Noneq_gas::
s_eval_density(Gas_data &Q)
{
    double R_hp = hp_gas_constant(Q);
    double R_e = Q.massf[iel_] * R_[iel_];

    Q.rho = Q.p / ( R_hp * Q.T[iT_] + R_e * Q.T[iTe_] );
    
    return SUCCESS;
}

double
Noneq_gas::
s_gas_constant(const Gas_data &Q, int &status)
{
    double R = 0.0;
    for ( size_t isp = 0; isp < R_.size(); ++isp )
	R += Q.massf[isp]*R_[isp];
    
    return R;
}

double 
Noneq_gas::
s_prho_ratio(const Gas_data &Q, int isp)
{
    if ( isp!=iel_ )
	return R_[isp]*Q.T[iT_];
    else
	return R_[iel_]*Q.T[iTe_];
}

// Just use Perfect_gas relations for the derivates as they should not be used

double
Noneq_gas::
s_dTdp_const_rho(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q,status);
    return 1.0/(Q.rho*R);
}

double
Noneq_gas::
s_dTdrho_const_p(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q,status);
    return (-1.0*Q.p)/(R*Q.rho*Q.rho);
}

double
Noneq_gas::
s_dpdrho_const_T(const Gas_data &Q, int &status)
{
    double R = s_gas_constant(Q,status);
    return R*Q.T[iT_];
}

double
Noneq_gas::
s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
{
    return ( isp !=  iel_ ) ? R_[isp] * Q.T[iT_] : R_[isp] * Q.T[iTe_];
}

double
Noneq_gas::
s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
{
    if ( iT_==iTe_ ) 
    	return s_gas_constant(Q,status) * Q.rho;
    else if ( itm==iTe_ )
    	return R_[iel_] * Q.massf[iel_] * Q.rho;
    else
    	return hp_gas_constant(Q) * Q.rho;
}

double
Noneq_gas::
hp_gas_constant(const Gas_data &Q)
{
    double R = 0.0;
    for ( int isp = 0; isp < iel_; ++isp )
	R += Q.massf[isp]*R_[isp];
    
    return R;
}

