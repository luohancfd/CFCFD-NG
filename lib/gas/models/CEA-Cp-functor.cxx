// Author: Rowan J. Gollan
// Date: 04-Nov-2008
// Place: Hampton, Virginia, USA

#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "CEA-Cp-functor.hh"

using namespace std;

CEA_Cp_functor::
CEA_Cp_functor(lua_State *L, double R)
    : R_(R)
{
    T_low_ = get_positive_number(L, -1, "T_low");
    T_high_ = get_positive_number(L, -1, "T_high");

    if ( T_high_ <= T_low_ ) {
	cout << "CEA_Cp_functor::CEA_Cp_functor()\n";
	cout << "Error initialising CEA_Cp_functor object.\n";
	cout << "T_high should be greater than T_low\n";
	cout << "T_high= " << T_high_ << " T_low= " << T_low_ << endl;
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    lua_getfield(L, -1, "coeffs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "CEA_Cp_functor::CEA_Cp_functor()\n";
	ost << "A 'coeffs' table is expected with the 'parameters' table, but it's not found.\n";
	input_error(ost);
    }

    int n = lua_objlen(L, -1);
    if ( n < 7 ) {
	cout << "CEA_Cp_functor::CEA_Cp_functor()\n";
	cout << "7 coefficients are required to specify the\n";
	cout << "curve for Cp in the 'coeffs' table.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    a_.resize(7);
    
    for ( size_t i = 0; i < a_.size(); ++i ) {
	lua_rawgeti(L, -1, i+1);
	a_[i] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }

    lua_pop(L, 1); // pops 'coeffs'
}

CEA_Cp_functor::
CEA_Cp_functor(const CEA_Cp_functor &c)
    : R_(c.R_), T_low_(c.T_low_), T_high_(c.T_high_), a_(c.a_)  {}

CEA_Cp_functor::
~CEA_Cp_functor() {}

CEA_Cp_functor*
CEA_Cp_functor::
clone() const
{
    return new CEA_Cp_functor(*this);
}

double
CEA_Cp_functor::
operator()(double T)
{
    // These catches give constant extrapolation
    // beyond range of validity.

    if ( T < T_low_ )
	T = T_low_;

    if ( T > T_high_ )
	T = T_high_;

    double Cp = a_[0]/(T*T) + a_[1]/T + a_[2] + a_[3]*T;
    Cp += a_[4]*T*T + a_[5]*T*T*T + a_[6]*pow(T, 4);
    return Cp*R_;
}
