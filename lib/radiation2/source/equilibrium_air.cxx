/** \file equilibrium_air.cxx
 *  \ingroup radiation
 *
 *  \author Rowan J. Gollan
 *  \version 13-Jan-07 : ported from work of Sep-Oct 2003
 *           06-Jul-09 : ported from lib/radiation
 *
 **/

#include <cmath>
#include <iostream>

#include "equilibrium_air.hh"
#include "radiation_constants.hh"

using namespace std;

EquilibriumAir::
~EquilibriumAir() {}

string EquilibriumAir::str() const
{
    return "EquilibriumAir";
}


double
EquilibriumAir::
integrated_emission_for_gas_state( Gas_data &Q, bool spectrally_resolved )
{
    UNUSED_VARIABLE(spectrally_resolved);
    
    double T = Q.T[0];		  // Assume T[0] is the gas temperature
    
    double kappa_P = Planck_absorption(Q);
    double j_total = kappa_P * RC_sigma_SI * pow(T, 4) / M_PI;

    return j_total;
}


double
EquilibriumAir::
Planck_absorption(Gas_data &Q)
{
    double T = Q.T[0];		  // Assume T[0] is the gas temperature
    
    const double rho_0 = 1.225; // kg/m^3 (sea-level density)
    double kappa_P = 7.94 * pow(Q.rho/rho_0, 1.10) * pow(T/1.0e4, 6.95);
    return kappa_P;
}

void
EquilibriumAir::
spectra_for_gas_state( Gas_data &Q, CoeffSpectra &X )
{
    // NOTE: - storing INTEGRATED j and kappa at single frequency interval
    //       - this allows this model to be used with all the radiation transport models
    
    // 1. resize spectral vectors to size of 1
    X.nu.resize(1);
    X.j_nu.resize(1);
    X.kappa_nu.resize(1);
    
    // 2. calculate INTEGRATED emission and absorption coefficients for this gas_state
    double j_total = integrated_emission_for_gas_state(Q, false);
    double kappa_P = Planck_absorption(Q);
    
    // 3. store in spectral vectors
    X.nu[0] = 1.0;
    X.j_nu[0] = j_total;
    X.kappa_nu[0] = kappa_P;
    
    return;
}

void
EquilibriumAir::
spectral_distribution_for_gas_state(Gas_data &Q, vector<double> &nus)
{
    UNUSED_VARIABLE(Q);
    nus.resize(1);
    nus[0] = 1.0;
    
    return;
}

void
EquilibriumAir::
write_line_widths( Gas_data &Q )
{
    cout << "EquilibriumAir::write_line_widths()" << endl
         << "This function is not applicable for this spectral model!" << endl;
    return;
}
