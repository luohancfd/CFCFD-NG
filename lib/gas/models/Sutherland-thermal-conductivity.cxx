// Author: Rowan J. Gollan
// Date: 10-Jul-2008

#include <cmath>
#include <stdexcept>

#include "../../util/source/lua_service.hh"
#include "Sutherland-thermal-conductivity.hh"

using namespace std;

Sutherland_thermal_conductivity::
Sutherland_thermal_conductivity(lua_State *L)
    : Thermal_conductivity_model()
{
    // Assume a table with the model paramaters is TOS.
    k0_ = get_positive_number(L, -1, "k_ref");
    T0_ = get_positive_number(L, -1, "T_ref");
    S_ = get_positive_number(L, -1, "S");

}

Sutherland_thermal_conductivity::
~Sutherland_thermal_conductivity() {}

double
Sutherland_thermal_conductivity::
s_eval_thermal_conductivity(const Gas_data &Q)
{
    return S_thermal_conductivity(Q.T[0], k0_, T0_, S_);
}

double S_thermal_conductivity(double T, double k0, double T0, double S)
{
    // Reference:
    // White, F.M. (2006)
    // Viscous Fluid Flow, 3rd edition
    // McGraw Hill International, New York
    // Equation 1-44b, p. 30
    return k0 * pow(T/T0, 3.0/2.0) * (T0 + S)/(T + S);
}
