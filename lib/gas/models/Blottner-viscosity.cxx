// Author: Elise J. Fahy
// Date: 24-Sep-2014

#include <iostream>
#include <cmath>

#include "../../util/source/lua_service.hh"
#include "Blottner-viscosity.hh"

using namespace std;

Blottner_viscosity::
Blottner_viscosity(lua_State *L)
    : Viscosity_model()
{
    // Assume a table with the model parameters is TOS.
    A_mu_ = get_positive_number(L, -1, "A_mu");
    B_mu_ = get_positive_number(L, -1, "B_mu");
    C_mu_ = get_positive_number(L, -1, "C_mu");
}

Blottner_viscosity::
~Blottner_viscosity() {}

double
Blottner_viscosity::
s_eval_viscosity(const Gas_data &Q)
{
    return S_viscosity(Q.T[0], A_mu_, B_mu_, C_mu_);
}

double S_viscosity(double T, double A_mu, double B_mu, double C_mu)
{
    // Reference:
    // Blottner, Johnson and Ellis (1971)
    // Chemically Reacting Viscous Flow Program for Multicomponent
    // Gas Mixtures
    // Sandia Laboratories Report SC-RR-70-754

    return 0.1 * exp((A_mu*log(T) + B_mu)*log(T) + C_mu);
}
