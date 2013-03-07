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

static int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

Polyatom_vibrational_mode::
Polyatom_vibrational_mode( int g, double omega, int Vmax )
 : g( g ), omega( omega ), Vmax( Vmax ) {}

int
Polyatom_vibrational_mode::
calculate_p( int iV )
{
    // compute statistical weight
    // Eq 40 p12 Capitelli (2005)
    int p = factorial( iV + g - 1 ) / factorial( iV ) / factorial( g - 1 );

    return p;
}

double
Polyatom_vibrational_mode::
eval_E_vib( int iV )
{
    /* implement the Dunham expansion to solve for the anharmonic vibrational energy
     * NOTE: currently only using the harmonic term
     *
     */
    double Gv = omega * ( double(g)/2.0 + double(iV) );

    return Gv;
}

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

    /* create the vmodes with degeneracy considered */
    vector<double> tmp;
    for ( size_t ivm=0; ivm<omega_e_vec.size(); ++ivm ) {
	// calculate level degeneracy
	int g_vib = 0;
	for ( size_t jvm=0; jvm<omega_e_vec.size(); ++jvm ) {
	    if ( fabs(omega_e_vec[jvm]-omega_e_vec[ivm])/omega_e_vec[ivm] < 1.0e-6 )
	        g_vib++;
	}
	// check that this frequency has not yet been added
	bool not_yet_considered = true;
	for ( size_t jvm=0; ivm<tmp.size(); ++jvm ) {
	    if ( fabs(omega_e_vec[ivm]-tmp[jvm])/tmp[jvm] < 1.0e-6 ) {
		not_yet_considered = false;
	        break;
	    }
	}
	if ( !not_yet_considered ) continue;
        int Vmax = this->calculate_V_max(ivm);
        vmodes.push_back( new Polyatom_vibrational_mode(g_vib,omega_e_vec[ivm],Vmax) );
	tmp.push_back( omega_e_vec[ivm] );
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

double
Polyatom_electronic_level::
eval_E_vib( vector<int> iV )
{
    double Gv = 0.0;
    for ( size_t im=0; im<iV.size(); ++im )
        Gv += vmodes[im]->eval_E_vib(iV[im]);

    return Gv;
}

double
Polyatom_electronic_level::
eval_E_rot( int ivm, int iV, int iJ )
{
    /* Expression for a spherical top from Capitelli 2005 */

    // Coupling terms
    double B_v = B0;

    // Energy
    double E_r = B_v*double(iJ)*(double(iJ)+1.0);

    return E_r;
}

void
Polyatom_electronic_level::
vib_loop( int im, vector<int> vib_state, vector< vector<int> > &results )
{
    for ( int iV=0; iV<vmodes[im]->get_Vmax(); ++iV ) {
	vib_state[im] = iV;
        if ( im==((int)vmodes.size()-1) ) {
            // FIXME: check if this state is permitted
            results.push_back( vib_state );
        }
        else
            vib_loop( im+1, vib_state, results );
    }
}
