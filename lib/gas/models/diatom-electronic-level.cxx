// Author: Daniel F. Potter
// Version: 24-Mar-2010
//          Ported from lib/radiation/source/diatomic_radiator.cxx

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "diatom-electronic-level.hh"
#include "physical_constants.hh"

#include "../../util/source/useful.h"

using namespace std;

Diatom_electronic_level::Diatom_electronic_level( vector<double> lev_data )
{
    E = lev_data[0] * PC_c * PC_h_SI;		// Convert cm^-1 to J
    r_e = lev_data[1];
    g = int(lev_data[2]);
    D = lev_data[3] * PC_c * PC_h_SI;
    omega_e = lev_data[4] * PC_c * PC_h_SI;
    xomega_e = lev_data[5] * PC_c * PC_h_SI;
    yomega_e = lev_data[6] * PC_c * PC_h_SI;
    zomega_e = lev_data[7] * PC_c * PC_h_SI;
    B_e = lev_data[8] * PC_c * PC_h_SI;
    alpha_e = lev_data[9] * PC_c * PC_h_SI;
    D_e = lev_data[10] * PC_c * PC_h_SI;
    beta_e = lev_data[11] * PC_c * PC_h_SI;
    A_spin = lev_data[12] * PC_c * PC_h_SI;
    Lambda = (int) lev_data[13];
    spin = (int) lev_data[14];
    // Check that a usable dissociation energy is present
    if ( D < 0.0 && xomega_e > 0.0 ) {
    	// use D ~ omega_e**2 / ( 4*xomega_e )
    	D = omega_e*omega_e / ( 4*xomega_e );
    }
    else if ( D < 0.0 ) {
        // set to 164,000K in Joules - E
        D = PC_k_SI * 1.64e5 - E;
    }    	    
    // Now calculate V_max
    calculate_V_max();
    // .. and calculate the J_max vector
    calculate_J_max();
}

void
Diatom_electronic_level::
calculate_V_max()
{
    // Calculate the maximum vibrational quantum number that is just below the 
    // dissociation limit, or when d G_v / d v goes negative
    // See Babou et al 2009
    V_max = -1;
    double G_v_im1 = eval_E_vib( 0 );
    for ( int iV=1; iV<9999; ++iV ) {
    	double G_v_i = eval_E_vib( iV );
    	if ( eval_E_vib( iV ) > D || G_v_i < G_v_im1 ) break;
    	V_max = iV;
    	G_v_im1 = G_v_i;
    }
    
    // Check if this value will cause any problems
    if ( test_V_max() ) {
    	cout << "Diatom_electronic_level::calculate_V_max()" << endl
    	     << "Fix this before going any further!" << endl;
    	exit( FAILURE );
    }
    
    return;
}

bool
Diatom_electronic_level::
test_V_max()
{
    // 1. Check if between 0 and 99
    if ( V_max < 0 || V_max > 99 ) {
    	cout << "Diatom_electronic_level::test_V_max()" << endl
    	     << "V_max = " << V_max << " which is outside the expected range" << endl;
    	return true;
    }
    
    // If we get here this V_max is okay
    return false;
}

double
Diatom_electronic_level::
eval_potential_curve( int iJ, double r )
{
    // Combined Morse and centrifugal potential
    // See Babou et al 2009 IJT 30:416-438 for details
    
    double J = double(iJ);
    double beta = omega_e / ( 4.0 * sqrt( B_e * D ) );
    double tmpC = - 2.0 * beta * ( r - r_e ) / r_e;
    double tmpA = 1.0 - exp( tmpC );
    
    return D * tmpA * tmpA + B_e * ( r_e / r ) * ( r_e / r ) * J * ( J + 1.0 );
}

double
Diatom_electronic_level::
eval_potential_curve_first_derivative( int iJ, double r )
{
    // NOTE: this function has been numerically verified
    
    double J = double(iJ);
    double beta = omega_e / ( 4.0 * sqrt( B_e * D ) );
    double tmpC = - 2.0 * beta * ( r - r_e ) / r_e;
    double dtmpC_dr = - 2.0 * beta / r_e;
    double tmpA = 1.0 - exp( tmpC );
    double dtmpA_dr = - exp( tmpC ) * dtmpC_dr;
    
    return 2.0 * D * tmpA * dtmpA_dr + 2.0 * B_e * r_e * r_e / r / r / r * J * ( J + 1.0 );
}

double
Diatom_electronic_level::
eval_potential_curve_second_derivative( int iJ, double r )
{
    double delta_r = r * 0.0001;
    double r_u = r + 0.5 * delta_r;
    double r_l = r - 0.5 * delta_r;
    
    return ( eval_potential_curve_first_derivative( iJ, r_u ) - eval_potential_curve_first_derivative( iJ, r_l ) ) / delta_r;
}

double
Diatom_electronic_level::
solve_for_potential_curve_r_max( int iJ )
{
    // The task: Find where dVdr = 0 and d2Vdr2 < 0
    // We will cheat as we know that there will be a mininum at r_e
    double r_lower = 1.001 * r_e;
    double r_upper = 5 * r_e;
    double r_mid = 0.5 * ( r_upper + r_lower );
    double tol = 1.0e-6;

    int max_iterations = 10000;
    for ( int i=0; i<=max_iterations; ++i ) {
    	if ( eval_potential_curve_first_derivative(iJ,r_lower) * eval_potential_curve_first_derivative(iJ,r_mid) > 0.0 ) r_lower = r_mid;
    	else r_upper = r_mid;
    	// cout << "iteration = " << i << ", r_mid = " << r_mid << ", fabs(r_upper - r_lower)/r_mid = " << fabs(r_upper - r_lower)/r_mid << endl;
    	r_mid = 0.5 * ( r_upper + r_lower );
    	if ( fabs(r_upper - r_lower)/r_mid < tol ) break;
    	if ( i==max_iterations ) {
    	    cout << "Diatom_electronic_level::solve_for_potential_curve_r_max()" << endl
    	         << "Reached maximum iterations with iJ = " << iJ << endl
    	         << "Bailing out!" << endl;
    	    exit( FAILURE );
    	}
    }
    
    // ensure that the second derivative is less than zero
    if ( eval_potential_curve_second_derivative( iJ, r_mid ) >= 0.0 ) {
    	cout << "Diatom_electronic_level::solve_for_potential_curve_r_max()" << endl
    	     << "r_max = " << r_mid << ", dVdr = " << eval_potential_curve_first_derivative(iJ,r_mid) << ", d2Vdr2 = " << eval_potential_curve_second_derivative( iJ, r_mid ) << endl
    	     << "The secondary derivative must be zero!" << endl
    	     << "Bailing out!" << endl;
    	exit( FAILURE );
    }
    
    return r_mid;
}

void
Diatom_electronic_level::
calculate_J_max()
{
    J_max.resize(V_max+1);
    for ( int iV=0; iV<V_max; ++iV ) {
    	J_max[iV] = -1;
    	for ( int iJ=0; iJ<999; ++iJ ) {
    	    // 1. Find r_max for this J value
    	    double r_max = solve_for_potential_curve_r_max( iJ );
    	    // 2. See if this rotational level exceeds the maximised potential energy curve
    	    double E_vr = eval_E_vib( iV ) + eval_E_rot( iV, iJ );
    	    double V_rmax = eval_potential_curve( iJ, r_max );
    	    if ( E_vr > V_rmax ) break;
    	    // else:
    	    J_max[iV] = iJ;
    	}
    	if ( J_max[iV] < 0 ) {
    	    J_max[iV] = J_max[iV-1];
    	}
    	if ( iV > 0 ) {
    	    if ( J_max[iV] > J_max[iV-1]) {
    	        J_max[iV] = J_max[iV-1];
    	    }
    	}
    	if ( test_J_max(iV) ) {
    	    ostringstream oss;
    	    oss << "Diatom_electronic_level::calculate_J_max()" << endl
    	        << "J_max[" << iV << "] = " << J_max[iV] << " failed the test." << endl;
    	}
    }
    
    return;
}

bool
Diatom_electronic_level::
test_J_max(int iV)
{
    // 1. Check if between 0 and 999
    if ( J_max[iV] < 0 || J_max[iV] > 999 ) {
    	cout << "Diatom_electronic_level::test_J_max()" << endl
    	     << "J_max[" << iV << "] = " << J_max[iV] << " which is outside the expected range" << endl;
    	return true;
    }
    
    // If we get here this J_max is okay
    return false;
}

double
Diatom_electronic_level::
eval_E_vib( int iV )
{
    double E_v;
    
    /* implement the Dunham expansion to solve for the anharmonic vibrational energy */
    E_v = omega_e * ( 0.5 + double(iV) ) - \
          xomega_e * pow ( 0.5 + double(iV), 2 ) + \
          yomega_e * pow ( 0.5 + double(iV), 3 ) + \
          zomega_e * pow ( 0.5 + double(iV), 4 );
    
    return E_v;
}


double
Diatom_electronic_level::
eval_E_rot( int iV, int iJ )
{
    /* Expression for a (singlet) non-rigid, vibrating symmetric top from ref:
       Huber and Herzberg (1950) Constants of Diatomic Molecules p 118, 169 */
    
    // Coupling terms
    
    double D_v = D_e + (double(iV) + 0.5) * beta_e;
    double B_v = B_e - (double(iV) + 0.5) * alpha_e;
    // Energy
    double E_r = B_v*iJ*(iJ+1.0) + (A_spin-B_v)*double(Lambda*Lambda) - D_v*iJ*iJ*(iJ+1.0)*(iJ+1.0);
    
    return E_r;
}

double
Diatom_electronic_level::
eval_B_v( int iV )
{
    double B_v = B_e - alpha_e * ( double(iV) + 0.5 );
    
    return B_v;
}

double
Diatom_electronic_level::
eval_Qr( double T_rot )
{
    double Q_rot = 0.0;
    for ( int iV=0; iV<V_max; ++iV ) {
        for ( int iJ=0; iJ<J_max[iV]; ++iJ ) {
            double E_rot = eval_E_rot( iV, iJ ) - eval_E_rot( iV, 0 );
            Q_rot += double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
        }
    }

    return Q_rot;
}

double
Diatom_electronic_level::
eval_QvQr( double T_vib, double T_rot )
{
    double Q_vib = 0.0;
    for ( int iV=0; iV<V_max; ++iV ) {
        double Q_rot = 0.0;
        for ( int iJ=0; iJ<J_max[iV]; ++iJ ) {
            double E_rot = eval_E_rot( iV, iJ );
            Q_rot += double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
        }
        double E_vib = eval_E_vib( iV );
        Q_vib += exp( - E_vib / PC_k_SI / T_vib ) * Q_rot;
    }
    return Q_vib;
}

double
Diatom_electronic_level::
eval_QeQvQr( double T_el, double T_vib, double T_rot )
{
    return double(g) * exp( - E / PC_k_SI / T_el ) * eval_QvQr( T_vib, T_rot );
}
