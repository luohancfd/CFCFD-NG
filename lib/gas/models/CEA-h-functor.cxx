// Author: Rowan J. Gollan
// Date: 04-Nov-2008
// Place: Hampton, Virginia, USA

#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "CEA-h-functor.hh"

using namespace std;

CEA_h_functor::
CEA_h_functor(lua_State *L, double Cp_low, double Cp_high, double R)
    : R_(R), Cp_low_(Cp_low), Cp_high_(Cp_high)
{
    T_low_ = get_positive_number(L, -1, "T_low");
    T_high_ = get_positive_number(L, -1, "T_high");

    if ( T_high_ <= T_low_ ) {
	cout << "CEA_h_functor::CEA_h_functor()\n";
	cout << "Error initialising CEA_h_functor object.\n";
	cout << "T_high should be greater than T_low\n";
	cout << "T_high= " << T_high_ << " T_low= " << T_low_ << endl;
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    lua_getfield(L, -1, "coeffs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "CEA_h_functor::CEA_h_functor()\n";
	ost << "A 'coeffs' table is expected with the 'parameters' table, but it's not found.\n";
	input_error(ost);
    }

    int n = lua_objlen(L, -1);
    if ( n < 8 ) {
	cout << "CEA_h_functor::CEA_h_functor()\n";
	cout << "8 coefficients are required to specify the\n";
	cout << "curve for h in the 'coeffs' table.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    a_.resize(8);
    
    for ( size_t i = 0; i < a_.size(); ++i ) {
	lua_rawgeti(L, -1, i+1);
	a_[i] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }

    lua_pop(L, 1); // pops 'coeffs' 
}

CEA_h_functor::
CEA_h_functor(const CEA_h_functor &h)
    : R_(h.R_), T_low_(h.T_low_), T_high_(h.T_high_),
      Cp_low_(h.Cp_low_), Cp_high_(h.Cp_high_), a_(h.a_)  {}

CEA_h_functor::
~CEA_h_functor() {}

CEA_h_functor*
CEA_h_functor::
clone() const
{
    return new CEA_h_functor(*this);
}

double
CEA_h_functor::
operator()(double T)
{
    if ( T < T_low_ ) {
	double h_low = eval(T_low_);
	double h = h_low - Cp_low_*(T_low_ - T);
	return h;
    }

    if ( T > T_high_ ) {
	double h_high = eval(T_high_);
	double h = h_high  + Cp_high_*(T - T_high_);
	return h;
    }
    
    return eval(T);
}

double
CEA_h_functor::
eval(double T)
{
    double h = -a_[0]/(T*T) + a_[1]/T * log(T) + a_[2];
    h += a_[3]*T/2.0 + a_[4]*T*T/3.0 + a_[5]*T*T*T/4.0;
    h += a_[6]*pow(T, 4)/5.0 + a_[7]/T;
    return h*R_*T;
}
