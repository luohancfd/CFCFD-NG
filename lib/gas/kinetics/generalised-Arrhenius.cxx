// Author: Rowan J. Gollan
// Date: 10-Oct-2008

#include <cmath>
#include <sstream>
#include <iostream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "generalised-Arrhenius.hh"
#include "../models/physical_constants.hh"

using namespace std;

Generalised_Arrhenius::
Generalised_Arrhenius(lua_State *L, Gas_model &g, double T_upper, double T_lower)
    : Reaction_rate_coefficient(T_upper, T_lower)
{
    A_ = get_positive_number(L, -1, "A");
    n_ = get_number(L, -1, "n");
    E_a_ = get_number(L, -1, "E_a");
}

Generalised_Arrhenius::
Generalised_Arrhenius( double A, double n, double E_a, double T_upper, double T_lower )
    : Reaction_rate_coefficient(T_upper, T_lower), A_(A), n_(n), E_a_(E_a) {}

Generalised_Arrhenius::
~Generalised_Arrhenius() {}

int
Generalised_Arrhenius::
s_eval(const Gas_data &Q)
{
    double T = Q.T[0];
    // Check T doesn't exceed defined limits
    if ( T > T_upper_ )
	T = T_upper_;
    if ( T < T_lower_ )
	T = T_lower_;
    
    k_ = A_ * pow(T, n_) * exp(-E_a_ / (PC_k_SI * T) );
    return SUCCESS;
}

int
Generalised_Arrhenius::
s_eval_from_T(const double T)
{
    k_ = A_ * pow(T, n_) * exp(-E_a_ / (PC_k_SI * T) );
    return SUCCESS;
}

Reaction_rate_coefficient* create_Generalised_Arrhenius_coefficient(lua_State *L, Gas_model &g, double T_upper, double T_lower)
{
    return new Generalised_Arrhenius(L, g, T_upper, T_lower);
}
