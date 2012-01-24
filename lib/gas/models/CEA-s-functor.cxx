// Author: Rowan J. Gollan
// Date: 04-Nov-2008
// Place: Hampton, Virginia, USA

#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "CEA-s-functor.hh"

using namespace std;

CEA_s_functor::
CEA_s_functor(lua_State *L, double Cp_low, double Cp_high, double R)
    : R_(R), Cp_low_(Cp_low), Cp_high_(Cp_high)
{

    T_low_ = get_positive_number(L, -1, "T_low");
    T_high_ = get_positive_number(L, -1, "T_high");

    if ( T_high_ <= T_low_ ) {
	cout << "CEA_s_functor::CEA_s_functor()\n";
	cout << "Error initialising CEA_s_functor object.\n";
	cout << "T_high should be greater than T_low\n";
	cout << "T_high= " << T_high_ << " T_low= " << T_low_ << endl;
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    lua_getfield(L, -1, "coeffs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "CEA_s_functor::CEA_s_functor()\n";
	ost << "A 'coeffs' table is expected with the 'parameters' table, but it's not found.\n";
	input_error(ost);
    }

    int n = lua_objlen(L, -1);
    if ( n < 9 ) {
	cout << "CEA_s_functor::CEA_s_functor()\n";
	cout << "9 coefficients are required to specify the\n";
	cout << "curve for s in the 'coeffs' table.\n";
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }
    
    a_.resize(9);
    
    for ( size_t i = 0; i < a_.size(); ++i ) {
	lua_rawgeti(L, -1, i+1);
	a_[i] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }

    lua_pop(L, 1); // pops 'coeffs'
}

CEA_s_functor::
CEA_s_functor(const CEA_s_functor &s)
    : R_(s.R_), T_low_(s.T_low_), T_high_(s.T_high_),
      Cp_low_(s.Cp_low_), Cp_high_(s.Cp_high_), a_(s.a_)  {}

CEA_s_functor::
~CEA_s_functor() {}

CEA_s_functor*
CEA_s_functor::
clone() const
{
    return new CEA_s_functor(*this);
}

double
CEA_s_functor::
operator()(double T)
{
    if ( T < T_low_ ) {
	double s_low = eval(T_low_);
	double s = s_low - Cp_low_ * log(T_low_/T);
	return s;
    }
    
    if ( T > T_high_ ) {
	double s_high = eval(T_high_);
	double s = s_high + Cp_high_ * log(T/T_high_);
	return s;
    }

    return eval(T);
}

double
CEA_s_functor::
eval(double T)
{
    double s = -a_[0]/(2.0*T*T) - a_[1]/T + a_[2]*log(T);
    s += a_[3]*T + a_[4]*T*T/2.0 + a_[5]*T*T*T/3.0;
    s += a_[6]*pow(T, 4)/4.0 + a_[8];
    return s*R_;
}
