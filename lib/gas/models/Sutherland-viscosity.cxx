// Author: Rowan J. Gollan
// Date: 10-Jul-2008

#include <iostream>
#include <cmath>

#include "../../util/source/lua_service.hh"
#include "Sutherland-viscosity.hh"

using namespace std;

Sutherland_viscosity::
Sutherland_viscosity(lua_State *L)
    : Viscosity_model()
{
    // Assume a table with the model parameters is TOS.
    mu0_ = get_positive_number(L, -1, "mu_ref");
    T0_ = get_positive_number(L, -1, "T_ref");
    S_ = get_positive_number(L, -1, "S");
}

Sutherland_viscosity::
~Sutherland_viscosity() {}

double
Sutherland_viscosity::
s_eval_viscosity(const Gas_data &Q)
{
    return S_viscosity(Q.T[0], mu0_, T0_, S_);
}

double S_viscosity(double T, double mu0, double T0, double S)
{
    // Reference:
    // White, F.M. (2006)
    // Viscous Fluid Flow, 3rd edition
    // McGraw Hill International, New York
    // Equation 1-36, p. 28
    return mu0 * pow(T/T0, 3.0/2.0) * (T0 + S)/(T + S);
}
