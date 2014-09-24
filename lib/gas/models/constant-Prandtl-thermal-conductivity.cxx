// Author: Elise J. Fahy
// Date: 24-Sep-2014

#include <cmath>
#include <stdexcept>

#include "../../util/source/lua_service.hh"
#include "constant-Prandtl-thermal-conductivity.hh"

using namespace std;

constant_Prandtl_thermal_conductivity::
constant_Prandtl_thermal_conductivity(lua_State *L)
    : Thermal_conductivity_model()
{
    // set constant Prandtl number here for now...
    Pr_ = 0.66

}

constant_Prandtl_thermal_conductivity::
~constant_Prandtl_thermal_conductivity() {}

double
constant_Prandtl_thermal_conductivity::
s_eval_thermal_conductivity(const Gas_data &Q)
{
    return S_thermal_conductivity(Q.mu, Cp, Pr_);
    // I don't know where Cp lives! But I don't think it's in Q...
}

double S_thermal_conductivity(double mu, double Cp, double Pr_)
{
    // using Prandtl number expression, rearranged for k
    return (mu * Cp) / Pr_
}
