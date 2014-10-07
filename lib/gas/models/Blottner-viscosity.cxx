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
    A_mu_ = get_number(L, -1, "A_mu");
    B_mu_ = get_number(L, -1, "B_mu");
    C_mu_ = get_number(L, -1, "C_mu");
}

Blottner_viscosity::
~Blottner_viscosity() {}

double
Blottner_viscosity::
s_eval_viscosity(const Gas_data &Q)
{
    // Reference:
    // Blottner, Johnson and Ellis (1971)
    // Chemically Reacting Viscous Flow Program for Multicomponent
    // Gas Mixtures
    // Sandia Laboratories Report SC-RR-70-754

    return 0.1 * exp((A_mu_*log(Q.T[0]) + B_mu_)*log(Q.T[0]) + C_mu_);
}
