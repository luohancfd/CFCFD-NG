// Author: Daniel F. Potter
// Version: 24-Mar-2010
//          Ported from lib/radiation/source/diatomic_radiator.cxx

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "polyatom-electronic-level.hh"
#include "physical_constants.hh"

#include "../../util/source/useful.h"

using namespace std;

Polyatom_vibrational_mode::
Polyatom_vibrational_mode( double omega, int Vmax, std::vector<int> Jmax )
 : omega( omega ), Vmax( Vmax ), Jmax( Jmax )
{}

Polyatom_electronic_level::
Polyatom_electronic_level( vector<double> lev_data )
{
    E = lev_data[0] * PC_c * PC_h_SI;
    r_e = lev_data[1];
    g = int(lev_data[2]);
    D = lev_data[3] * PC_c * PC_h_SI;

    A0 = lev_data[4] * PC_c * PC_h_SI;
    B0 = lev_data[5] * PC_c * PC_h_SI;
    C0 = lev_data[6] * PC_c * PC_h_SI;

    sigma = lev_data[7];
    sigma_r = lev_data[8];

    for ( size_t i=9; i<lev_data.size(); ++i )
         omega_e_vec.push_back( lev_data[i] * PC_c * PC_h_SI );

    for ( size_t ivm=0; ivm<omega_e_vec.size(); ++ivm ) {
        int Vmax = this->calculate_V_max(ivm);
        vector<int> Jmax(Vmax);
        for ( int iV=0; iV<Vmax; ++iV ) {
            Jmax[iV] = this->calculate_J_max(ivm,iV);
        }
        vmodes.push_back( new Polyatom_vibrational_mode(omega_e_vec[ivm],Vmax,Jmax) );
    }
}

Polyatom_electronic_level::
~Polyatom_electronic_level()
{
    for ( size_t ivm=0; ivm<vmodes.size(); ++ivm ) {
        delete vmodes[ivm];
    }
}

int
Polyatom_electronic_level::
calculate_V_max( int ivm )
{
    return 10;
}

int
Polyatom_electronic_level::
calculate_J_max( int ivm, int iV )
{
    return 100;
}

double
Polyatom_electronic_level::
eval_E_vib( int ivm, int iV )
{
    double omega = vmodes[ivm]->omega;

    /* implement the Dunham expansion to solve for the anharmonic vibrational energy
     * NOTE: currently only using the harmonic term
     *
     */
    double E_v = omega * ( 0.5 + double(iV) );

    return E_v;
}


double
Polyatom_electronic_level::
eval_E_rot( int ivm, int iV, int iJ )
{
    /* Expression for a spherical top from Capitelli 2005 */

    // Coupling terms

    double B_v = B0;    // FIXME: is this right?
    // Energy
    double E_r = B_v*double(iJ)*(double(iJ)+1.0);

    return E_r;
}
