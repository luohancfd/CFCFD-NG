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
Generalised_Arrhenius(lua_State *L, Gas_model &g)
    : Reaction_rate_coefficient()
{
    A_ = get_positive_number(L, -1, "A");
    n_ = get_number(L, -1, "n");
    E_a_ = get_number(L, -1, "E_a");
}

Generalised_Arrhenius::
Generalised_Arrhenius( double A, double n, double E_a )
    : Reaction_rate_coefficient(), A_(A), n_(n), E_a_(E_a) {}

Generalised_Arrhenius::
~Generalised_Arrhenius() {}

int
Generalised_Arrhenius::
s_eval(const Gas_data &Q)
{
    k_ = A_ * pow(Q.T[0], n_) * exp(-E_a_ / (PC_k_SI * Q.T[0]) );
    return SUCCESS;
}

int
Generalised_Arrhenius::
s_eval_from_T(const double T)
{
    k_ = A_ * pow(T, n_) * exp(-E_a_ / (PC_k_SI * T) );
    return SUCCESS;
}

Reaction_rate_coefficient* create_Generalised_Arrhenius_coefficient(lua_State *L, Gas_model &g)
{
    return new Generalised_Arrhenius(L, g);
}
