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

double
Polyatom_vibrational_mode::
eval_Q_from_T( double T )
{
    /* Partition function assuming harmonic oscillator model */
    double Q_HO = 1.0 / ( 1.0 - exp( -omega / PC_k_SI / T ) );

    return Q_HO;
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
    return 3;
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

void
Polyatom_electronic_level::
vib_loop( int im, vector<int> &vib_state, vector< vector<int> > &results )
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

Spherical_top_polyatom_electronic_level::
Spherical_top_polyatom_electronic_level( vector<double> lev_data )
 : Polyatom_electronic_level( lev_data )
{}

Spherical_top_polyatom_electronic_level::
~Spherical_top_polyatom_electronic_level()
{}

double
Spherical_top_polyatom_electronic_level::
eval_E_rot( int iJ, int iK )
{
    /* Expression for a spherical top from Capitelli 2005 */

    // Quantum numbers
    double J = double(iJ);
    UNUSED_VARIABLE(iK);

    // Energy
    double E_r = B0*J*(J+1.0);

    return E_r;
}

double
Spherical_top_polyatom_electronic_level::
eval_Q_from_T( double T )
{
    double QvQr;

#   if FULL_ROVIBRATIONAL_SUMMATION
    // Sum over all rovibrational states
    // FIXME: this currently cannot be done correctly as V_max and J_max are just guesses
    QvQr = 0.0;
    vector< vector<int> > vib_states;
    vector<int> vib_state(vmodes.size(),0);
    vib_loop(0,vib_state,vib_states);
    for ( size_t i=0; i<vib_states.size(); ++i) {
	double J_max = 100; // compute_J_max( vib_states[i] );
	double E_vib = eval_E_vib( vib_states[i] ); // - E_vib_0 ?
	for ( int iJ=0; iJ<=J_max; ++iJ ) {
            double E_rot = eval_E_rot( iJ );
            double E_state = E + E_vib + E_rot;
            QvQr += g * exp( - E_state / PC_k_SI / T );
	}
    }
#   else
    // use rigid-rotator and infinite harmonic oscillator approximations
    // assume a linear molecule for the moment...
    QvQr = 1.0 / sigma * PC_k_SI * T / B0;
    for ( size_t im=0; im<vmodes.size(); ++im )
        QvQr *= vmodes[im]->eval_Q_from_T(T);
#   endif

    return QvQr;
}

Symmetrical_top_polyatom_electronic_level::
Symmetrical_top_polyatom_electronic_level( vector<double> lev_data )
 : Polyatom_electronic_level( lev_data )
{}

Symmetrical_top_polyatom_electronic_level::
~Symmetrical_top_polyatom_electronic_level()
{}

double
Symmetrical_top_polyatom_electronic_level::
eval_E_rot( int iJ, int iK )
{
    /* See Eq. 43 in Capitelli 2005 */

    // Quantum numbers
    double J = double(iJ);
    double K = double(iK);

    // Energy
    double E_r = B0*J*(J+1.0) + ( A0 - B0 )*K*K;

    return E_r;
}

double
Symmetrical_top_polyatom_electronic_level::
eval_Q_from_T( double T )
{
    double QvQr;

#   if FULL_ROVIBRATIONAL_SUMMATION
    // Sum over all rovibrational states
    // FIXME: this currently cannot be done correctly as V_max and J_max are just guesses
    QvQr = 0.0;
    vector< vector<int> > vib_states;
    vector<int> vib_state(vmodes.size(),0);
    vib_loop(0,vib_state,vib_states);
    for ( size_t i=0; i<vib_states.size(); ++i) {
	double J_max = 100; // compute_J_max( vib_states[i] );
	double E_vib = eval_E_vib( vib_states[i] ); // - E_vib_0 ?
	for ( int iJ=0; iJ<=J_max; ++iJ ) {
	    for ( int iK=-iJ; iK<=iJ; ++ik ) {
                double E_rot = eval_E_rot( iJ, iK );
                double E_state = E + E_vib + E_rot;
                QvQr += g * exp( - E_state / PC_k_SI / T );
	    }
	}
    }
#   else
    // use rigid-rotator and infinite harmonic oscillator approximations
    // assume a linear molecule for the moment...
    QvQr = sqrt( PC_k_SI*T*PC_k_SI*T*PC_k_SI*T * M_PI / ( A0 * B0 * C0 ) );
    for ( size_t im=0; im<vmodes.size(); ++im )
        QvQr *= vmodes[im]->eval_Q_from_T(T);
#   endif

    return QvQr;
}

Asymmetric_top_polyatom_electronic_level::
Asymmetric_top_polyatom_electronic_level( vector<double> lev_data )
 : Polyatom_electronic_level( lev_data )
{}

Asymmetric_top_polyatom_electronic_level::
~Asymmetric_top_polyatom_electronic_level()
{}

double
Asymmetric_top_polyatom_electronic_level::
eval_E_rot( int iJ, int iK )
{
    /* See Eq. 44 in Capitelli 2005 */

    // Quantum numbers
    double J = double(iJ);
    double K = double(iK);

    // Energy
    double E_r = 0.5*(B0+C0)*J*(J+1.0) + ( A0 - B0 )*K*K;

    return E_r;
}

double
Asymmetric_top_polyatom_electronic_level::
eval_Q_from_T( double T )
{
    double QvQr;

#   if FULL_ROVIBRATIONAL_SUMMATION
    // Sum over all rovibrational states
    // FIXME: this currently cannot be done correctly as V_max and J_max are just guesses
    QvQr = 0.0;
    vector< vector<int> > vib_states;
    vector<int> vib_state(vmodes.size(),0);
    vib_loop(0,vib_state,vib_states);
    for ( size_t i=0; i<vib_states.size(); ++i) {
	double J_max = 100; // compute_J_max( vib_states[i] );
	double E_vib = eval_E_vib( vib_states[i] ); // - E_vib_0 ?
	for ( int iJ=0; iJ<=J_max; ++iJ ) {
	    for ( int iK=-iJ; iK<=iJ; ++ik ) {
                double E_rot = eval_E_rot( iJ, iK );
                double E_state = E + E_vib + E_rot;
                QvQr += g * exp( - E_state / PC_k_SI / T );
	    }
	}
    }
#   else
    // use rigid-rotator and infinite harmonic oscillator approximations
    // assume a linear molecule for the moment...
    QvQr = sqrt( PC_k_SI*T*PC_k_SI*T*PC_k_SI*T * M_PI / ( A0 * B0 * C0 ) );
    for ( size_t im=0; im<vmodes.size(); ++im )
        QvQr *= vmodes[im]->eval_Q_from_T(T);
#   endif

    return QvQr;
}
