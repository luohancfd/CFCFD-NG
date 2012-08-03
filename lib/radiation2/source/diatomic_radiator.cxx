/** \file diatomic_radiator.cxx
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: New radiation2 version
 *
 *  \brief Definitions for diatomic radiator classes
 *
 **/

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "../../util/source/lua_service.hh"
#include "../../gas/models/gas_data.hh"
#include "diatomic_radiator.hh"
#include "diatomic_system.hh"
#include "radiation_constants.hh"
#include "photaura.hh"

using namespace std;

DiatomicElecLev::DiatomicElecLev( int ilev, vector<double> lev_data, bool homonuclear, double I_spin )
 : ElecLev( ilev, lev_data[0]*RC_c*RC_h_SI, int(lev_data[2]) ),
   homonuclear( homonuclear ), I_spin( I_spin )
{
    r_e = lev_data[1];				// Angstroms?
    D = lev_data[3] * RC_c * RC_h_SI;		// Convert cm^-1 to J
    omega_e = lev_data[4] * RC_c * RC_h_SI;
    xomega_e = lev_data[5] * RC_c * RC_h_SI;
    yomega_e = lev_data[6] * RC_c * RC_h_SI;
    zomega_e = lev_data[7] * RC_c * RC_h_SI;
    B_e = lev_data[8] * RC_c * RC_h_SI;
    alpha_e = lev_data[9] * RC_c * RC_h_SI;
    D_e = lev_data[10] * RC_c * RC_h_SI;
    beta_e = lev_data[11] * RC_c * RC_h_SI;
    A_spin = lev_data[12] * RC_c * RC_h_SI;
    lambda = (int) lev_data[13];
    spin = (int) lev_data[14];
    if ( homonuclear ) { 
	P_pm = int(lev_data[15]);
	P_gu = int(lev_data[16]);
    }
    // Now calculate V_max
    calculate_V_max();
    // .. and calculate the J_max vector
    calculate_J_max();
    // Read in optional Sigma - Sigma spin splitting parameters
    if ( int(lev_data.size()) > 999917 ) {
    	gamma_V.resize(V_max+1,0.0);
    	for ( int iV=0; iV<(int(lev_data.size()) - 17); ++iV ) {
    	    gamma_V[iV] = double(lev_data[17+iV])*1.0e-3 * RC_c * RC_h_SI;	// Convert 1e3 cm^-1 to J
    	    if ( ECHO_RAD_INPUT > 1 ) {
    	    	cout << "gamma[" << iV << "] = " << lev_data[17+iV] << endl;
    	    }
    	}
    }
}

void
DiatomicElecLev::
calculate_V_max()
{
    // Calculate the maximum vibrational quantum number that is just below the 
    // dissociation limit, or when d G_v / d v goes negative
    // See Babou et al 2009
    V_max = -1;
    double G_v_im1 = calculate_E_vib( 0 );
    for ( int iV=1; iV<9999; ++iV ) {
    	double G_v_i = calculate_E_vib( iV );
    	if ( calculate_E_vib( iV ) > D || G_v_i < G_v_im1 ) break;
    	V_max = iV;
    	G_v_im1 = G_v_i;
    }
    
    // Check if this value will cause any problems
    if ( test_V_max() ) {
    	cout << "DiatomicElecLev::calculate_V_max()" << endl
    	     << "Fix this before going any further!" << endl;
    	exit( FAILURE );
    }
    
    if ( ECHO_RAD_INPUT > 1 ) cout << "V_max = " << V_max << endl;

    return;
}

bool
DiatomicElecLev::
test_V_max()
{
    // 1. Check if between 0 and 99
    if ( V_max < 0 || V_max > 99 ) {
    	cout << "DiatomicElecLev::test_V_max()" << endl
    	     << "V_max = " << V_max << " which is outside the expected range" << endl;
    	return true;
    }
    
    // 2. Check QvQr behaviour at low temperature
    double E_v = calculate_E_vib(V_max);
    double QvQr_test = calculate_Q_rot(200.0, V_max) * exp( - E_v / ( RC_k_SI * 200.0 ) );
    if ( isinf( QvQr_test ) ) {
    	cout << "DiatomicElecLev::test_V_max()" << endl
    	     << "V_max = " << V_max << " leads to infinite QvQr term at low temperatures" << endl;
    	return true;
    }
    
    // If we get here this V_max is okay
    return false;
}

double
DiatomicElecLev::
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
DiatomicElecLev::
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
DiatomicElecLev::
eval_potential_curve_second_derivative( int iJ, double r )
{
    double delta_r = r * 0.0001;
    double r_u = r + 0.5 * delta_r;
    double r_l = r - 0.5 * delta_r;
    
    return ( eval_potential_curve_first_derivative( iJ, r_u ) - eval_potential_curve_first_derivative( iJ, r_l ) ) / delta_r;
}

double
DiatomicElecLev::
solve_for_potential_curve_r_max( int iJ )
{
    // The task: Find where dVdr = 0 and d2Vdr2 < 0
    // We will cheat as we know that there will be a mininum at r_e
    double r_lower = 1.001 * r_e;
    double r_upper = 5 * r_e;
    double r_mid = 0.5 * ( r_upper + r_lower );
    
    int max_iterations = 1000;
    for ( int i=0; i<=max_iterations; ++i ) {
    	if ( eval_potential_curve_first_derivative(iJ,r_lower) * eval_potential_curve_first_derivative(iJ,r_mid) > 0.0 ) r_lower = r_mid;
    	else r_upper = r_mid;
    	// cout << "iteration = " << i << ", r_mid = " << r_mid << ", fabs(r_upper - r_lower)/r_mid = " << fabs(r_upper - r_lower)/r_mid << endl;
    	r_mid = 0.5 * ( r_upper + r_lower );
    	if ( fabs(r_upper - r_lower)/r_mid < 1.0e-6 ) break;
    	if ( i==max_iterations ) {
    	    cout << "DiatomicElecLev::solve_for_potential_curve_r_max()" << endl
    	         << "Reached maximum iterations with iJ = " << iJ << endl
    	         << "Bailing out!" << endl;
    	    exit( FAILURE );
    	}
    }
    
    // ensure that the second derivative is less than zero
    if ( eval_potential_curve_second_derivative( iJ, r_mid ) >= 0.0 ) {
    	cout << "DiatomicElecLev::solve_for_potential_curve_r_max()" << endl
    	     << "r_max = " << r_mid << ", dVdr = " << eval_potential_curve_first_derivative(iJ,r_mid) << ", d2Vdr2 = " << eval_potential_curve_second_derivative( iJ, r_mid ) << endl
    	     << "The secondary derivative must be zero!" << endl
    	     << "Bailing out!" << endl;
    	exit( FAILURE );
    }
    
    return r_mid;
}

void
DiatomicElecLev::
calculate_J_max()
{
    J_max.resize(V_max+1);
    for ( int iV=0; iV<V_max; ++iV ) {
    	J_max[iV] = -1;
    	for ( int iJ=0; iJ<999; ++iJ ) {
    	    // 1. Find r_max for this J value
    	    double r_max = solve_for_potential_curve_r_max( iJ );
    	    // 2. See if this rotational level exceeds the maximised potential energy curve
    	    double E_vr = calculate_E_vib( iV ) + calculate_E_rot_HundA( iV, 2*iJ );
    	    double V_rmax = eval_potential_curve( iJ, r_max );
    	    if ( E_vr > V_rmax ) break;
    	    // else:
    	    J_max[iV] = iJ;
    	}
    	if ( J_max[iV] < 0 ) {
    	    if ( DEBUG_RAD > 0 ) {
		cout << "DiatomicElecLev::calculate_J_max()" << endl
		     << "The calculation failed for iV = " << iV << " (J_max = " << J_max[iV] << ")" << endl
		     << "Setting J_max to the value for the previous vibrational level (J_max = " << J_max[iV-1] << ")" << endl;
	    }
    	    J_max[iV] = J_max[iV-1];
    	}
    	if ( iV > 0 ) {
    	    if ( J_max[iV] > J_max[iV-1]) {
    	    	if ( DEBUG_RAD > 0 ) {
		    cout << "DiatomicElecLev::calculate_J_max()" << endl
			 << "The calculation failed for iV = " << iV << " (J_max = " << J_max[iV] << ")" << endl
			 << "Setting J_max to the value for the previous vibrational level (J_max = " << J_max[iV-1] << ")" << endl;
		}
    	        J_max[iV] = J_max[iV-1];
    	    }
    	}
    	if ( ECHO_RAD_INPUT > 1 ) cout << "J_max[" << iV << "] = " << J_max[iV] << endl;
    	if ( test_J_max(iV) ) {
    	    ostringstream oss;
    	    oss << "DiatomicElecLev::calculate_J_max()" << endl
    	        << "J_max[" << iV << "] = " << J_max[iV] << " failed the test." << endl;
    	}
    }
    
    return;
}

bool
DiatomicElecLev::
test_J_max(int iV)
{
    // 1. Check if between 0 and 999
    if ( J_max[iV] < 0 || J_max[iV] > 999 ) {
    	cout << "DiatomicElecLev::test_J_max()" << endl
    	     << "J_max[" << iV << "] = " << J_max[iV] << " which is outside the expected range" << endl;
    	return true;
    }
    
    // If we get here this J_max is okay
    return false;
}

void
DiatomicElecLev::
write_potential_curves( std::string species_name, int ilev )
{
    for ( int iJ=0; iJ < J_max[0]; iJ+=10 ) {
    	ostringstream oss;
    	oss << species_name << "_ilev-" << ilev << "_iJ-" << iJ << ".txt";
    	ofstream ofile;
    	ofile.open(oss.str().c_str());
        ofile << "# Column 1: r (Ang)" << endl
              << "# Column 2: r/r_e" << endl
              << "# Column 3: V(r) (J)" << endl;
        ofile << setprecision(12) << scientific << showpoint;
        double delta_f = 0.001;
        for ( double f=delta_f; f<=5.0; f+=delta_f ) {
            double r = f * r_e;
            double V = eval_potential_curve( iJ, r );
            ofile << setw(20) << r << setw(20) << r / r_e << setw(20) << V << endl;
        }
        ofile.close();
    }
    
    return;
}

string
DiatomicElecLev::
string()
{
    ostringstream ost;
    ost << setw(15) << E / ( RC_c * RC_h_SI )    
        << setw(5)  << g 
        << setw(5)  << spin 
        << setw(5)  << lambda 
        << setw(10) << V_max            
        << setw(18) << omega_e / ( RC_c * RC_h_SI ) 
        << setw(18) << xomega_e / ( RC_c * RC_h_SI ) 
        << setw(18) << yomega_e / ( RC_c * RC_h_SI )
        << setw(18) << zomega_e / ( RC_c * RC_h_SI ) 
        << setw(15) << B_e / ( RC_c * RC_h_SI )     
        << setw(18) << alpha_e / ( RC_c * RC_h_SI )  
        << setw(18) << D_e / ( RC_c * RC_h_SI )     
        << setw(18) << beta_e / ( RC_c * RC_h_SI )
        << setw(10) << A_spin / ( RC_c * RC_h_SI );
	
    return ost.str();
}

double
DiatomicElecLev::
calculate_equilibrium_Q_total( double T )
{
    return calculate_Q_el( T ) * calculate_Q_vib( T, T );
}

double
DiatomicElecLev::
calculate_Q_vib( double T_vib, double T_rot )
{
    /* Full summation expression for vibrational partition function */
    double E_v, Q_rot;
    double QvibQrot = 0.0;
    
    for (int iV=0; iV<V_max; iV++) {
    	// Vibrational energy
        E_v = calculate_E_vib(iV);
        // Rotational partition function
        Q_rot = calculate_Q_rot(T_rot, iV);
        QvibQrot += Q_rot * exp( - E_v / ( RC_k_SI * T_vib ) );
    }
    
    return QvibQrot;
}

double
DiatomicElecLev::
calculate_and_store_Q_vib( double T_vib, double T_rot )
{
    /* Full summation expression for vibrational partition function */
    double E_v, Q_rot;
    QvQr = 0.0;
    
    for (int iV=0; iV<V_max; iV++) {
    	// Vibrational energy
        E_v = calculate_E_vib(iV);
        // Rotational partition function
        Q_rot = calculate_Q_rot(T_rot, iV);
        QvQr += Q_rot * exp( - E_v / ( RC_k_SI * T_vib ) );
    }
    
    return QvQr;
}

double
DiatomicElecLev::
calculate_E_vib( int iV )
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
DiatomicElecLev::
calculate_Q_rot(double T, int iV)
{
    double Q_rot=0.0;
    
#   if EXACT_Q_ROT
    // Full summation over rotational levels (for testing only)
    // NOTE: - using singlet energy expression
    //       - initialisation should have ensured iV < J_max.size()
    double E_r;
    int itJ_min = 2 * lambda;
    int itJ_max = 2 * (int) J_max[iV];
    
    for (int itJ=itJ_min; itJ<itJ_max; itJ+=2) {
	E_r = calculate_E_rot_HundA(iV,itJ);
	Q_rot += double(itJ + 1) * exp( - E_r / ( RC_k_SI * T ) );
    }
#   else
    // Approximation expression from [ h.c.Bv << k.Tr ] case
    double B_v = B_e - (double(iV) + 0.5) * alpha_e;
    // Vibrational and rotational constants have already been converted to Joules 
    Q_rot = RC_k_SI * T / B_v;
#   endif

    // Apply the homonuclear 'sigma' factor
    if ( homonuclear ) Q_rot *= 0.5;
    
    return Q_rot;
}

double
DiatomicElecLev::
calculate_E_rot_HundA( int iV, int itJ )
{
    /* Expression for a (singlet) non-rigid, vibrating symmetric top from ref:
       Huber and Herzberg (1950) Constants of Diatomic Molecules p 118, 169 */
    
    // Quantum numbers
    double iJ = double(itJ)/2.0;
    
    // Coupling terms
    double B_v = B_e - (double(iV) + 0.5) * alpha_e;
    double D_v = D_e + (double(iV) + 0.5) * beta_e;
    
    // Energy
    // NOTE: the A - B_v term is not needed as we use energies referenced from the ground state
    // double E_r = B_v*iJ*(iJ+1.0) + (A_spin-B_v)*double(lambda*lambda) - D_v*iJ*iJ*(iJ+1.0)*(iJ+1.0);
    double E_r = B_v*iJ*(iJ+1.0) - D_v*iJ*iJ*(iJ+1.0)*(iJ+1.0);
    
    return E_r;
}

double
DiatomicElecLev::
calculate_E_rot_SigmaTriplet( int iV, int itJ )
{
    /* Expression for the effective centre of a triplet from ref:
       Arnold, Whiting and Lyle (1969) JQRST Vol. 9 pp 775-798         */
    
    // Quantum numbers
    double iJ = double(itJ)/2.0;
    
    // Coupling terms
    double B_v = B_e - (double(iV) + 0.5) * alpha_e;
    double D_v = D_e + (double(iV) + 0.5) * beta_e;
    double Y = A_spin / B_v;
    double Z_2 = (double(lambda*lambda)*Y*(Y-1.0)-4.0/9.0-2.0*(iJ*(iJ+1.0)))
		/(3.0*(double(lambda*lambda)*Y*(Y-4.0)+4.0/3.0+4.0*iJ*(iJ+1.0)));
    
    // Energy
    double E_r = B_v * ( iJ*(iJ+1.0) + 4.0 * Z_2 )
		-  D_v * (iJ+0.5)*(iJ+0.5)*(iJ+0.5)*(iJ+0.5);
    
    return E_r;
}

double
DiatomicElecLev::
calculate_E_rot_HundBDoublet( int iV, int itJ, int itSig )
{
    /* Expression for Hund case b (multiplet Sigma state) from ref:
       Huber and Herzberg (1950) Constants of Diatomic Molecules p 222 */
    
    // Quantum numbers
    double iJ = double(itJ)/2.0;
    double iSig = double(itSig)/2.0;
    double iK = iJ - iSig;
    
    // Coupling terms
    double B_v = B_e - (double(iV) + 0.5) * alpha_e;
    double D_v = D_e + (double(iV) + 0.5) * beta_e;
    
    // Energy
    double E_r = B_v*iK*(iK+1.0) - D_v*iK*iK*(iK+1.0)*(iK+1.0) + iSig * gamma_V[iV] * ( iK + 0.5 - iSig );
    
    return E_r;
}

// Todo:
// calculate_E_rot_HundBTriplet
/* Huber and Herzberg (1950) Constants of Diatomic Molecules p 223 */

double
DiatomicElecLev::
calculate_E_rot_HundABDoublet( int iV, int itJ, int itSig )
{
    /* Expression for Doublet Hund a/b intermediate case from ref:
       Huber and Herzberg (1950) Constants of Diatomic Molecules p 232 */
    
    // Quantum numbers
    double iJ = double(itJ)/2.0;
    double iSig = double(itSig)/2.0;
    
    // Coupling terms
    double B_v = B_e - (double(iV) + 0.5) * alpha_e;
    double D_v = D_e + (double(iV) + 0.5) * beta_e;
    double Y = A_spin / B_v;
    
    // Energy
    double E_r = B_v * ( (iJ+0.5)*(iJ+0.5) - double(lambda*lambda)
	    	- iSig * sqrt( 4.0*(iJ+0.5)*(iJ+0.5)
		+ Y * ( Y - 4.0 ) * double(lambda*lambda) ) ) - D_v*pow(iJ+0.5-iSig,4);
    
    return E_r;
}

double
DiatomicElecLev::
calculate_E_rot_HundABTriplet( int iV, int itJ, int itSig )
{
    /* Expression for Triplet Hund a/b intermediate case from ref:
       Huber and Herzberg (1950) Constants of Diatomic Molecules p 235 */
    
    // Quantum numbers
    double iJ = double(itJ)/2.0;
    double iSig = double(itSig)/2.0;
    
    // Coupling terms
    double B_v = B_e - (double(iV) + 0.5) * alpha_e;
    double D_v = D_e + (double(iV) + 0.5) * beta_e;
    double Y = A_spin / B_v;
    double Z_1 = double(lambda)*double(lambda)*Y*( Y - 4.0 ) + 4.0/3.0 + 4.0*iJ*(iJ+1.0);
    double Z_2 = 1.0/(3.0*Z_1)*(double(lambda)*double(lambda)*Y*(Y-1.0)-4.0/9.0-2.0*iJ*(iJ+1));
    
    // Energy
    double E_r = B_v*(iJ*(iJ+1.0) - iSig*sqrt(Z_1) + Z_2 * ( 4.0 - fabs(iSig)*6.0 ) ) - D_v*(iJ+0.5-iSig);
    
    return E_r;
}

double
DiatomicElecLev::
calculate_N_vib( double T_vib, double T_rot, int iV )
{
    double E_vib = calculate_E_vib( iV );
    double Q_rot = calculate_Q_rot( T_rot, iV );
    
    double N_vib = N * Q_rot * exp( - E_vib / ( RC_k_SI * T_vib ) ) / QvQr;
    
    return N_vib;
}

double
DiatomicElecLev::
calculate_N_rot_HundA( double T_vib, double T_rot, int iV, int itJ, bool apply_LAF )
{
    double E_vib = calculate_E_vib( iV );
    double E_rot = calculate_E_rot_HundA( iV, itJ );
    double g_rot = double(itJ+1);
    
    double N_rot = N * g_rot * exp( - E_vib / ( RC_k_SI * T_vib ) - E_rot / ( RC_k_SI * T_rot ) ) / QvQr;
    
    if ( homonuclear ) {
    	N_rot *= 0.5;
    	if ( apply_LAF ) N_rot *= line_alternation_factor( itJ );
    }
    
    return N_rot;
}

double
DiatomicElecLev::
calculate_N_rot_SigmaTriplet( double T_vib, double T_rot, int iV, int itJ, bool apply_LAF )
{
    double E_vib = calculate_E_vib( iV );
    double E_rot = calculate_E_rot_SigmaTriplet( iV, itJ );
    double g_rot = double(itJ+1);
    
    double N_rot = N * g_rot * exp( - E_vib / ( RC_k_SI * T_vib ) - E_rot / ( RC_k_SI * T_rot ) ) / QvQr;
    
    if ( homonuclear ) {
    	N_rot *= 0.5;
    	if ( apply_LAF ) N_rot *= line_alternation_factor( itJ );
    }
    
    return N_rot;
}

double
DiatomicElecLev::
calculate_N_rot_HundBDoublet( double T_vib, double T_rot, int iV, int itJ, int itSig, bool apply_LAF )
{
    double E_vib = calculate_E_vib( iV );
    double E_rot = calculate_E_rot_HundBDoublet( iV, itJ, itSig );
    double g_rot = double(itJ+1);
    
    double N_rot = N * g_rot * exp( - E_vib / ( RC_k_SI * T_vib ) - E_rot / ( RC_k_SI * T_rot ) ) / QvQr;
    
    if ( homonuclear ) {
    	N_rot *= 0.5;
    	if ( apply_LAF ) N_rot *= line_alternation_factor( itJ, itSig );
    }
    
    return N_rot;
}

double
DiatomicElecLev::
calculate_N_rot_HundABDoublet( double T_vib, double T_rot, int iV, int itJ, int itSig, bool apply_LAF )
{
    double E_vib = calculate_E_vib( iV );
    double E_rot = calculate_E_rot_HundABDoublet( iV, itJ, itSig );
    double g_rot = double(itJ+1);
    
    double N_rot = N * g_rot * exp( - E_vib / ( RC_k_SI * T_vib ) - E_rot / ( RC_k_SI * T_rot ) ) / QvQr;
    
    if ( homonuclear ) {
    	N_rot *= 0.5;
    	if ( apply_LAF ) N_rot *= line_alternation_factor( itJ, itSig );
    }
    
    return N_rot;
}

double
DiatomicElecLev::
calculate_N_rot_HundABTriplet( double T_vib, double T_rot, int iV, int itJ, int itSig, bool apply_LAF )
{
    double E_vib = calculate_E_vib( iV );
    double E_rot = calculate_E_rot_HundABTriplet( iV, itJ, itSig );
    double g_rot = double(itJ+1);
    
    double N_rot = N * g_rot * exp( - E_vib / ( RC_k_SI * T_vib ) - E_rot / ( RC_k_SI * T_rot ) ) / QvQr;
    
    if ( homonuclear ) {
    	N_rot *= 0.5;
    	if ( apply_LAF ) N_rot *= line_alternation_factor( itJ, itSig );
    }
    
    return N_rot;
}

double
DiatomicElecLev::
line_alternation_factor( int itJ, int itSig )
{
    // 0. If lambda is not equal to 0 then return 1.0
    if ( lambda!=0 ) return 1.0;
    
    // 1. Calculate K and J_star
    int K = (itJ - itSig)/2;
    int J_star = itJ / 2;	// This should be rounded down from J
    int m1_pow_J_star = ( J_star%2==0 ) ? 1 : -1;
    
    // 2. Calculate total parity of the rovibronic level
    int P_pm_rot = 0;
    if ( lambda==0 ) {
    	// Sigma state
	if ( P_pm==1 ) {
	    if ( K%2 == 0 ) P_pm_rot = 1;
	    else P_pm_rot = -1;
	}
	else {
	    if ( K%2 == 0 ) P_pm_rot = -1;
	    else P_pm_rot = 1;
	}
    }
    else {
    	cout << "DiatomicElecLev::line_alternation_factor()" << endl
    	     << "Only implemented for Sigma states!" << endl
    	     << "Bailing out!" << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    // 3. Calculate P_ef (see Laux p 13)
    int P_ef = 0;
    if ( m1_pow_J_star * P_pm * P_pm_rot > 0 ) P_ef = 1;
    else P_ef = -1;
    
    // 4. Calculate the line alternation factor (see Laux p 19)
    // NOTE: multiplying LAF x 2 as we apply homonuclear factor for Q_rot_i
    
    if ( P_ef * P_gu * m1_pow_J_star == 1 ) {
    	if ( I_spin==0.5 ) return 2.0 * I_spin / ( 2.0 * I_spin + 1.0 );
    	else return 2.0 * ( I_spin + 1.0 ) / ( 2.0 * I_spin + 1.0 );
    }
    else {
    	if ( I_spin==0.5 ) return 2.0 * ( I_spin + 1.0 ) / ( 2.0 * I_spin + 1.0 );
    	else return 2.0 * I_spin / ( 2.0 * I_spin + 1.0 );
    }
}


double
DiatomicElecLev::
calculate_B_v( int iV )
{
    double B_v = B_e - alpha_e * ( double(iV) + 0.5 );
    
    return B_v;
}

int
DiatomicElecLev::
get_J_max( int iV )
{
    return J_max[iV];
}

int
DiatomicElecLev::
get_lambda()
{
    return lambda;
}

int
DiatomicElecLev::
get_spin()
{
    return spin;
}

double
DiatomicElecLev::
get_A_spin()
{
    return A_spin;
}

bool
DiatomicElecLev::
has_gamma_Vs()
{
    if ( gamma_V.size() > 0 ) return true;
    else return false;
}

double
DiatomicElecLev::
get_gamma_V( int iV )
{
    return gamma_V[iV];
}

double
DiatomicElecLev::
get_D_e()
{
    return D_e;
}

int
DiatomicElecLev::
get_V_max()
{
    return V_max;
}

int
DiatomicElecLev::
eval_kronecker_delta()
{
    if ( lambda==0 ) return 1;
    else return 0;
}

bool
DiatomicElecLev::
get_homonuclear()
{
    return homonuclear;
}

double
DiatomicElecLev::
get_D()
{
    return D;
}

double
DiatomicElecLev::
get_B_e()
{
    return B_e;
}

DiatomicRadiator::
DiatomicRadiator( lua_State * L, string name )
 : Radiator(L, name)
{
    iTv = get_int( L, -1, "iTv" );
    
    iTr = get_int( L, -1, "iTr" );
    
    I_spin = get_int( L, -1, "I_spin" );
    
    D = get_number( L, -1, "eta_D" );
    D *= RC_c * RC_h_SI;		// Convert cm**-1 -> J
    
    homonuclear = false;
    if ( name.find("2")!=string::npos ) homonuclear = true;
    
    read_elevel_data( L );
    
    read_system_data( L );
    
    read_photoionization_data( L );
}

DiatomicRadiator::
~DiatomicRadiator()
{
    for ( int ilev=0; ilev<nlevs; ++ilev )
    	delete elevs[ilev];
    
    for ( size_t isys=0; isys<systems.size(); ++isys )
	delete systems[isys];
}

void
DiatomicRadiator::
set_e_index( int iel )
{
    Radiator::set_e_index( iel );
    
    for ( size_t isys=0; isys<systems.size(); ++isys ) {
    	systems[isys]->e_index = iel;
    }
}


ElecLev *
DiatomicRadiator::
get_elev_pointer( int ie )
{
    return elevs[ie];
}

DiatomicElecLev *
DiatomicRadiator::
get_diatomic_elev_pointer( int ie )
{
    return elevs[ie];
}

string
DiatomicRadiator::
get_system_name( int isys )
{
    return systems[isys]->name;
}

DiatomicSystem *
DiatomicRadiator::
get_system_pointer( int isys )
{
    return systems[isys];
}

double
DiatomicRadiator::
get_D()
{
    return D;
}

double
DiatomicRadiator::
calculate_vibronic_wavenumber( int ie_u, int ie_l, int Vu, int Vl )
{
    for ( size_t isys=0; isys<systems.size(); ++isys ) {
    	if ( systems[isys]->ie_l == ie_l && systems[isys]->ie_u == ie_u ) {
    	    return systems[isys]->band_pointer( Vu, Vl )->calculate_average_frequency() / RC_c;
    	}
    }
    
    cout << "DiatomicRadiator::calculate_vibronic_wavenumber()\n";
    cout << "System[ie_u=" << ie_u << ", ie_l=" << ie_l << "] not found.\n";
    cout << "Bailing out!\n";
    exit( BAD_INPUT_ERROR );
}

double
DiatomicRadiator::
calculate_vibronic_wavelength( std::string sys_name, int Vu, int Vl )
{
    for ( size_t isys=0; isys<systems.size(); ++isys ) {
    	if ( systems[isys]->name == sys_name ) {
    	    double nu_ul = systems[isys]->band_pointer( Vu, Vl )->calculate_average_frequency();
    	    return nu2lambda(nu_ul);
    	}
    }
    
    cout << "DiatomicRadiator::calculate_vibronic_wavelength()\n";
    cout << "System: " << sys_name << " not found.\n";
    cout << "Bailing out!\n";
    exit( BAD_INPUT_ERROR );
}

int
DiatomicRadiator::
calculate_kronecker_delta_for_elev( int ie )
{
    return dynamic_cast<DiatomicElecLev*>(elevs[ie])->eval_kronecker_delta();
}

double
DiatomicRadiator::
calculate_system_transition_probability( int isys, double Tv )
{
    return systems[isys]->calculate_transition_probability( Tv );
}

void
DiatomicRadiator::
read_elevel_data( lua_State * L )
{
    lua_getfield(L, -1, "level_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "DiatomicRadiator::read_elevel_data()\n";
	ost << "Error locating 'elec_levels' table" << endl;
	input_error(ost);
    }
    
    nlevs = get_positive_int(L, -1, "n_levels");
    elevs.resize( nlevs );
    
    lua_getfield(L, -1, "isp_list");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "DiatomicRadiator::read_elevel_data()\n";
	ost << "Error locating isp_list table" << endl;
	input_error(ost);
    }
    vector<double> isp_list;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	isp_list.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L,1);	// pop isp_list
    
    // FIXME: do something with the isp_list!
    
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
	ostringstream lev_oss;
	lev_oss << "ilev_" << ilev;
	lua_getfield(L, -1, lev_oss.str().c_str());
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "DiatomicRadiator::read_elevel_data()\n";
	    ost << "Error locating " << lev_oss << " table" << endl;
	    input_error(ost);
	}
	vector<double> lev_data;
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    lev_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop ilev
	// Create the electronic level
	elevs[ilev] = new DiatomicElecLev(ilev,lev_data,homonuclear,I_spin);
#       if 0
        elevs[ilev]->write_potential_curves( name, ilev );
#       endif
    }
    lua_pop(L,1); 	// pop elec_levels
    
    if ( ECHO_RAD_INPUT > 1 ) {
	cout << "nlevs = " << nlevs << endl;
	cout << setw(10) << "[ilev]"
	     << setw(15) << "[E_el (1/cm)]"
	     << setw(5)  << "[g_el]"
	     << setw(5)  << "[spin]"
	     << setw(5)  << "[lambda]"
	     << setw(10) << "[V_max]"        
	     << setw(18) << "[omega_e (1/cm)]"
	     << setw(18) << "[xomega_e (1/cm)]"
	     << setw(18) << "[yomega_e (1/cm)]"
	     << setw(18) << "[zomega_e (1/cm)]"
	     << setw(15) << "[B_e (1/cm)]"
	     << setw(18) << "[alpha_e (1/cm)]"
	     << setw(18) << "[D_e (1/cm)]"
	     << setw(18) << "[beta_e (1/cm)]" 
	     << setw(10) << "[A_spin (1/cm)]" << endl;
	for ( int ilev=0; ilev<nlevs; ++ilev )
	    cout << "ilev_" << ilev << " = " << elevs[ilev]->string() << endl;
    }
    
    return;
}

void
DiatomicRadiator::
read_system_data( lua_State * L )
{
    lua_getfield(L, -1, "systems");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "DiatomicRadiator::DiatomicRadiator()\n";
	ost << "Error locating 'systems' table" << endl;
	input_error(ost);
    }
    
    lua_getfield(L, -1, "systems_list" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "DiatomicRadiator::DiatomicRadiator():\n";
	ost << "Error in the declaration of 'systems_list': a table is expected.\n";
	input_error(ost);
    }
    
    nsys = lua_objlen(L, -1);
    systems.resize(nsys);
    
    for ( int isys = 0; isys < nsys; ++isys ) {
	lua_rawgeti(L, -1, isys+1); // A Lua list is offset one from the C++ vector index
	const char* sys = luaL_checkstring(L, -1);
	lua_pop(L, 1);		// pop system
	
	// system_list is top-of-stack, reach one back to get 'systems' field
	lua_getfield(L, -2, sys);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "read_system_data():\n";
	    ost << "Error in the declaration of system: " << sys << " - a table is expected.\n";
	    input_error(ost);
	}
	// Electronic level pointers
	int ie_u = get_int( L, -1, "ie_u" );
	int ie_l = get_int( L, -1, "ie_l" );
	DiatomicElecLev * elev_u = dynamic_cast<DiatomicElecLev*>(elevs[ie_u]);
	DiatomicElecLev * elev_l = dynamic_cast<DiatomicElecLev*>(elevs[ie_l]);
	
	// Now construct the system
	systems[isys] = create_new_diatomic_system( L, sys, elev_u, elev_l, iT, iTe, iTv, iTr, m_w, I );
	
	lua_pop(L,1); 	// pop system
    }
    
    lua_pop(L,1); 		// pop system_list
    
    lua_pop(L,1); 		// pop systems field
    
    return;
}

void
DiatomicRadiator::
calculate_Q_int( Gas_data &Q )
{
    double T_el = Q.T[iTe];
    double T_vib = Q.T[iTv];
    double T_rot = Q.T[iTr];
    double Q_el_i, QvQr;
    Q_int = 0.0;
    Q_el = 0.0;
    
    // Q_int = sum ( Q_el * QvQr )
    for (int ilev=0; ilev<nlevs; ilev++) {
    	QvQr = elevs[ilev]->calculate_and_store_Q_vib(T_vib, T_rot);
    	Q_el_i = elevs[ilev]->calculate_and_store_Q_el(T_el);
    	elevs[ilev]->set_Q_int( Q_el_i * QvQr );
	Q_el += Q_el_i;
	Q_int += QvQr * Q_el_i;
    }
    
    return;
}

double
DiatomicRadiator::
calculate_total_equil_partition_function( double T )
{
    // 1. Translational contribution
    double Q_tr = pow( 2.0 * M_PI * m_w / RC_Na * RC_k_SI * T / RC_h_SI / RC_h_SI, 1.5 );
    
    // 2. Rovibronic (coupled electronic) contribution
    double Q_elec = 0.0;
    for (int ilev=0; ilev<nlevs; ++ilev) {
	Q_elec += elevs[ilev]->calculate_Q_el(T) * elevs[ilev]->calculate_Q_vib(T, T);
    }
    
    // 3. Return the product of the modal contributions
    return Q_tr * Q_elec;
}

void
DiatomicRadiator::
initialise_mechanisms( Gas_data &Q )
{
    // 0. Pre-calculate n_hvy, n_elecs, mw_av
    // NOTE: we are assuming Q.p includes the electron contribution, and that 
    //       Q.p_e is present and correct
    //       [if a 1T model is being used p_e will be zero!]
    double N_elecs = 0.0;
    if ( e_index >= 0 )
    	N_elecs = Q.p_e * RC_Na/(RC_R_u*Q.T[iTe]) * 1.0e-6;
    double N_hvy = (Q.p - Q.p_e ) * RC_Na/(RC_R_u*Q.T[iT]) * 1.0e-6;
    double mw_av = ( Q.rho * RC_Na ) / ( ( N_hvy + N_elecs ) * 1.0e6 ) * 1.0e3;
    
    for (size_t isys=0; isys<systems.size(); ++isys) {
	systems[isys]->initialize( Q.T[iT], Q.T[iTe], Q.p, N_hvy, N_elecs, mw_av );
    }
    
    return;
}

double
DiatomicRadiator::
calculate_unresolved_emission_coefficient( Gas_data &Q )
{
    double j_ul = 0.0;
    
    /* Pre-initialised, just sum the contributions */

    for (size_t isys=0; isys<systems.size(); ++isys)
	j_ul += systems[isys]->calculate_j_ul(Q.T[iTv],Q.T[iTr]);

    return j_ul;
}

double
DiatomicRadiator::
calculate_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u )
{
    double j_ul = 0.0;
    
    /* Pre-initialised, just sum the contributions */

    for (size_t isys=0; isys<systems.size(); ++isys)
	j_ul += systems[isys]->calculate_OV_j_ul(Q.T[iTv],Q.T[iTr],wavel_switch,Lambda_l,Lambda_u);

    return j_ul;
}

void
DiatomicRadiator::
spectral_distribution( vector<double> &nus )
{
    return;
}

void
DiatomicRadiator::
calculate_spectrum( Gas_data &Q, CoeffSpectra &X )
{
    // Calculate the unresolved total emission to measure the importance of each band
    double j_av = 0.0;
    for (size_t isys=0; isys<systems.size(); ++isys) {
    	j_av += systems[isys]->calculate_j_ul(Q.T[iTv],Q.T[iTr]) /
    			double(systems[isys]->lRe_dim * systems[isys]->uRe_dim);
    }
    
    // Loop over the systems and add the contributions
    for (size_t isys=0; isys<systems.size(); ++isys) {
   	systems[isys]->calculate_spectrum(X,Q.T[iTv],Q.T[iTr],j_av);
    }
    
    return;
}

string
DiatomicRadiator::
line_width_string( Gas_data &Q )
{
    // 0. Pre-calculate n_hvy, n_elecs, mw_av
    // NOTE: we are assuming Q.p includes the electron contribution, and that 
    //       Q.p_e is present and correct
    //       [if a 1T model is being used p_e will be zero!]
    double N_elecs = 0.0;
    if ( e_index >= 0 )
    	N_elecs = Q.p_e * RC_Na/(RC_R_u*Q.T[iTe]) * 1.0e-6;
    double N_hvy = (Q.p - Q.p_e ) * RC_Na/(RC_R_u*Q.T[iT]) * 1.0e-6;
    double mw_av = ( Q.rho * RC_Na ) / ( ( N_hvy + N_elecs ) * 1.0e6 ) * 1.0e3;
    
    // 1. Loop through systems, call line_width_string function
    string lws = "";
    for (size_t isys=0; isys<systems.size(); ++isys) {
	lws += systems[isys]->line_width_string( Q.T[iT], Q.T[iTe], Q.p, N_hvy, N_elecs, mw_av );
    }
    
    return lws;
}

/************************** BoltzDiatomicRadiator **************************/

BoltzDiatomicRadiator::
BoltzDiatomicRadiator( lua_State * L, std::string name )
 : DiatomicRadiator(L, name) {}
 
BoltzDiatomicRadiator::
~BoltzDiatomicRadiator() {}

void
BoltzDiatomicRadiator::
calculate_n_e( Gas_data &Q )
{
    double n_total = Q.massf[isp] * Q.rho / m_w * RC_Na;	// convert kg/m**3 -> particles/m**3
    
    for (int ilev=0; ilev<nlevs; ++ilev) {
	elevs[ilev]->set_N( n_total * elevs[ilev]->get_Q_int() / Q_int );
#       if DEBUG_RAD > 0
	cout << "N_el[" << ilev << "] = " << elevs[ilev]->get_N() << endl;
#	endif
    }
}

/*************************** QSSDiatomicRadiator ***************************/

QSSDiatomicRadiator::
QSSDiatomicRadiator( lua_State * L, std::string name )
 : DiatomicRadiator(L, name)
{
    lua_getfield(L, -1, "QSS_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "QSSDiatomicRadiator::QSSDiatomicRadiator()\n";
	ost << "Error locating 'QSS_model' table" << endl;
	input_error(ost);
    }
    
    // 0. Read-in the mininum temperature to apply the CR model
    T_lower = get_number( L, -1, "T_lower" );
    
    // 1. Create the nonequilibrium electronic levels
    lua_getfield(L, -1, "noneq_elevs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "QSSDiatomicRadiator::QSSDiatomicRadiator()\n";
	ost << "Error locating 'noneq_elevs' table" << endl;
	input_error(ost);
    }
    
    vector<int> noneq_ilevs;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	noneq_ilevs.push_back( (int) luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop noneq_elevs table
    
    lua_getfield(L, -1, "noneq_elev_labels");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "QSSDiatomicRadiator::QSSDiatomicRadiator()\n";
	ost << "Error locating 'noneq_elev_labels' table" << endl;
	input_error(ost);
    }
    
    vector<string> noneq_elev_labels;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	noneq_elev_labels.push_back( (string) luaL_checkstring(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop noneq_labels table
    
    int inc_eq_elevs = get_int(L,-1,"inc_eq_elevs");
    
    for ( size_t ne_ilev=0; ne_ilev<noneq_ilevs.size(); ++ne_ilev ) {
    	int ilev = noneq_ilevs[ne_ilev];
    	// a. create vector of pointers to equilibriated levels
    	vector<ElecLev*> eq_elevs;
    	if ( inc_eq_elevs ) {
	    int next_ilev = (int) elevs.size();				// if this is the last noneq level
	    if ( ne_ilev!=noneq_ilevs.size()-1 ) next_ilev = noneq_ilevs[ne_ilev+1];	// if not
	    for ( int eq_ilev=ilev+1; eq_ilev<next_ilev; ++eq_ilev ) 
		eq_elevs.push_back( elevs[eq_ilev] );
	}
    	// b. create the noneq level
    	noneq_elevs.push_back( new NoneqElecLev( ilev, ne_ilev, noneq_elev_labels[ne_ilev], elevs[ilev], eq_elevs ) );
    }
    
    // 2. Create the collision-radiative reaction mechanisms
    int nreactions = create_reactions( L );
    cout << " - Created " << nreactions << " for QSSDiatomicRadiator: " << name << endl;
    
    lua_pop(L, 1);	// pop QSS_model table
    
    // 3. Initialise the working matrices and valarrays
    dGdy = new Valmatrix();
    dGdy->resize( noneq_elevs.size(), noneq_elevs.size() );
    C.resize( noneq_elevs.size(), 0.0 );
    y_out.resize( noneq_elevs.size(), 0.0 );
}
 
QSSDiatomicRadiator::
~QSSDiatomicRadiator()
{
    for ( size_t i=0; i<noneq_elevs.size(); ++i )
    	delete noneq_elevs[i];
    
    for ( size_t i=0; i<reactions.size(); ++i )
    	delete reactions[i];
    
    delete dGdy;
}

void
QSSDiatomicRadiator::
set_radiator_pointers( std::vector<Radiator*> radiators )
{
    if ( name.size()==3 ) {
    	ostringstream oss;
    	oss << "QSSDiatomicRadiator::set_radiator_pointers()" << endl
    	    << "Cannot handle diatoms with three character names!" << endl;
    	input_error( oss );
    }
    
    string atom_A_name = name.substr(0,1);
    string atom_B_name = name.substr(1,1);
    if ( atom_B_name=="2" && name.size()==2 ) atom_B_name = atom_A_name;
    else if ( atom_B_name=="2" && name.size()==7 ) atom_B_name = atom_A_name + "_plus";
    
    atom_A = get_radiator_pointer_from_name( radiators, atom_A_name );
    atom_B = get_radiator_pointer_from_name( radiators, atom_B_name );
    elec = get_radiator_pointer_from_name( radiators, "e_minus" );
    
    for ( size_t ir=0; ir<reactions.size(); ++ir ) {
    	if ( reactions[ir]->get_type()=="ElectronImpactExcitation" ) {
    	    dynamic_cast<ElectronImpactExcitation*>(reactions[ir])->set_electron_pointer( elec );
    	}
    	else if ( reactions[ir]->get_type()=="ElectronImpactDissociation" ) {
    	    dynamic_cast<ElectronImpactDissociation*>(reactions[ir])->set_atom_A_pointer( atom_A );
    	    dynamic_cast<ElectronImpactDissociation*>(reactions[ir])->set_atom_B_pointer( atom_B );
    	    dynamic_cast<ElectronImpactDissociation*>(reactions[ir])->set_electron_pointer( elec );
    	}
    	else if ( reactions[ir]->get_type()=="HeavyParticleImpactDissociation" ) {
    	    string M_name = dynamic_cast<HeavyParticleImpactDissociation*>(reactions[ir])->M_name;
    	    Radiator * M = 0;
    	    for ( size_t irad=0; irad<radiators.size(); ++irad ) {
		if ( radiators[irad]->name==M_name ) {
		    M = radiators[irad];
		    break;
		}
	    }
    	    dynamic_cast<HeavyParticleImpactDissociation*>(reactions[ir])->set_atom_A_pointer( atom_A );
    	    dynamic_cast<HeavyParticleImpactDissociation*>(reactions[ir])->set_atom_B_pointer( atom_B );
    	    dynamic_cast<HeavyParticleImpactDissociation*>(reactions[ir])->set_heavy_particle_pointer( M );
    	}
    	else if ( reactions[ir]->get_type()=="HeavyParticleImpactExcitation" ) {
    	    string M_name = dynamic_cast<HeavyParticleImpactExcitation*>(reactions[ir])->M_name;
    	    Radiator * M = 0;
    	    for ( size_t irad=0; irad<radiators.size(); ++irad ) {
		if ( radiators[irad]->name==M_name ) {
		    M = radiators[irad];
		    break;
		}
	    }
    	    dynamic_cast<HeavyParticleImpactExcitation*>(reactions[ir])->set_heavy_particle_pointer( M );
    	}
    }
    
    return;
}

void
QSSDiatomicRadiator::
level_population_file( Gas_data &Q, int index )
{
    atom_A->calc_partition_functions(Q);
    atom_B->calc_partition_functions(Q);
    calculate_Q_int(Q);
    calculate_n_e(Q);
    
    ofstream ofile;
    ostringstream oss;
    oss << name << "_level_populations-" << index << ".txt";
    string fname = oss.str();
    ofile.open(fname.c_str(),ios::out);
    ofile << setprecision(12) << scientific << showpoint;
    ofile << name + "_level_populations.txt" << endl
    	  << "# Column 1: Level" << endl
    	  << "# Column 2: Level energy (eV)" << endl
	  << "# Column 3: Number density (particles/m**3)" << endl
	  << "# Column 4: Number density / degeneracy (particles/m**3)" << endl
	  << "# Column 5: Boltzmann population / degeneracy" << endl
	  << "# Column 6: Dissociation equilibrium population / degeneracy" << endl;
double delta_N_DE_N_boltz = 0.0;
double delta_N_QSS_N_boltz = 0.0;
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
    	ofile << setw(20) << ilev
    	      << setw(20) << get_elev_pointer(ilev)->get_E() / RC_e_SI 
    	      << setw(20) << get_elev_pointer(ilev)->get_N()
    	      << setw(20) << get_elev_pointer(ilev)->get_N() / get_elev_pointer(ilev)->get_g()
    	      << setw(20) << eval_Boltzmann_population_for_level(Q,ilev) / get_elev_pointer(ilev)->get_g()
    	      << setw(20) << eval_DE_population_for_level(Q,ilev) / get_elev_pointer(ilev)->get_g()
    	      << endl;
delta_N_DE_N_boltz += eval_DE_population_for_level(Q,ilev) / eval_Boltzmann_population_for_level(Q,ilev) - 1.0;
delta_N_QSS_N_boltz += get_elev_pointer(ilev)->get_N() / eval_Boltzmann_population_for_level(Q,ilev) -1.0;
    }
    ofile.close();
    
cout << name << ": delta N_DE / N_boltz = " << delta_N_DE_N_boltz / double(nlevs) << ", delta N_QSS / N_boltz = " << delta_N_QSS_N_boltz / double(nlevs) << endl;
    return;
}

double
QSSDiatomicRadiator::
eval_Boltzmann_population_for_level( Gas_data &Q, int ilev )
{
    double N_total = Q.massf[isp] * Q.rho / m_w * RC_Na;	// convert kg/m**3 -> particles/m**3
    
    return N_total * elevs[ilev]->get_Q_int() / Q_int;
}

double
QSSDiatomicRadiator::
eval_DE_population_for_level( Gas_data &Q, int ilev )
{
    // 1. Number densities
    double N_A = Q.massf[atom_A->isp] * Q.rho / atom_A->m_w * RC_Na;		// convert kg/m**3 -> particles/m**3
    double N_B = Q.massf[atom_B->isp] * Q.rho / atom_B->m_w * RC_Na;
    
    // 2. Partition functions
    double Q_level = this->eval_translational_partition_function(Q) * elevs[ilev]->get_Q_int() * exp( D / RC_k_SI / Q.T[iT] );
    double Q_A = atom_A->eval_translational_partition_function(Q) * atom_A->Q_int;
    double Q_B = atom_B->eval_translational_partition_function(Q) * atom_B->Q_int;
    
    // 3. Evaluate the DE equation
    return N_A * N_B * Q_level / Q_A / Q_B;
}

string
QSSDiatomicRadiator::
CR_model_latex_string()
{
    ostringstream oss;
    oss << "\\begin{table}[!p]" << endl
        << " \\centering" << endl
        << " \\begin{threeparttable}" << endl
        << " \\label{tab:CR-excite}" << endl
        << " \\begin{tabular*}{1.0\\textwidth}%" << endl
        << "     {@{\\extracolsep{\\fill}}lcccc}" << endl
        << " \\hline \\hline Reaction                                                            & $A$ (cm$^{3}$/s)                   &  $n$    & $E_{a}$ (K) &  Source \\\\" << endl;
        
    oss << " \\hline  \\multicolumn{1}{l}{\\emph{Electron impact excitation}} \\\\" << endl;
    
    for ( size_t ir=0; ir<reactions.size(); ++ir ) {
    	if ( reactions[ir]->get_type()=="ElectronImpactExcitation" ) {
    	    cout << "reaction " << ir << " is ElectronImpactExcitation" << endl;
    	    cout << "latex string = " << reactions[ir]->get_latex_string() << endl;
    	    oss << reactions[ir]->get_latex_string();
    	}
    }
    
    oss << " \\hline  \\multicolumn{1}{l}{\\emph{Electron impact dissociation}} \\\\" << endl;
    
    for ( size_t ir=0; ir<reactions.size(); ++ir ) {
    	if ( reactions[ir]->get_type()=="ElectronImpactDissociation" ) {
    	    cout << "reaction " << ir << " is ElectronImpactDissociation" << endl;
    	    cout << "latex string = " << reactions[ir]->get_latex_string() << endl;
    	    oss << reactions[ir]->get_latex_string();
    	}
    }
    
    oss << " \\hline  \\multicolumn{1}{l}{\\emph{Heavy-particle impact excitation}} \\\\" << endl;
    
    for ( size_t ir=0; ir<reactions.size(); ++ir ) {
    	if ( reactions[ir]->get_type()=="HeavyParticleImpactExcitation" ) {
    	    oss << reactions[ir]->get_latex_string();
    	}
    }
    
    oss << " \\hline  \\multicolumn{1}{l}{\\emph{Heavy-particle impact dissociation}} \\\\" << endl;
    
    for ( size_t ir=0; ir<reactions.size(); ++ir ) {
    	if ( reactions[ir]->get_type()=="HeavyParticleImpactDissociation" ) {
    	    oss << reactions[ir]->get_latex_string();
    	}
    }
    
    oss << " \\hline  \\multicolumn{1}{l}{\\emph{Radiative transitions}} \\\\" << endl;
    oss << "           \\hline Reaction                                                                      & $A$ (s$^{-1}$)             &              &                         & Source \\\\" << endl;
    
    
    for ( size_t ir=0; ir<reactions.size(); ++ir ) {
    	if ( reactions[ir]->get_type()=="RadiativeTransition" ) {
    	    cout << "reaction " << ir << " is RadiativeTransition" << endl;
    	    cout << "latex string = " << reactions[ir]->get_latex_string() << endl;
    	    oss << reactions[ir]->get_latex_string();
    	}
    }
    
    string tmp;
    string lrad_name;
    get_latex_species_piecies( name, lrad_name, tmp );
    
    oss << " \\hline" << endl
        << " \\end{tabular*}" << endl
        << " \\end{threeparttable}" << endl
        << " \\caption{Collisional-radiative model for " << lrad_name << ".}" << endl
        << " \\label{tab:" << name << "_CR_model}" << endl
        << "\\end{table}" << endl;
    
    
    return oss.str();
}

int
QSSDiatomicRadiator::
create_reactions( lua_State * L )
{
    lua_getfield(L, -1, "reactions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "QSSDiatomicRadiator::QSSDiatomicRadiator()\n";
	ost << "Error locating 'reactions' table" << endl;
	input_error(ost);
    }
    
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	string type = get_string( L, -1, "type" );
	if ( type=="heavy_particle_impact_excitation" ) {
	    reactions.push_back( new HeavyParticleImpactExcitation( L, this ) );
	}
	else if ( type=="electron_impact_excitation" ) {
	    reactions.push_back( new ElectronImpactExcitation( L, this ) );
	}
	else if ( type=="heavy_particle_impact_dissociation" ) {
	    reactions.push_back( new HeavyParticleImpactDissociation( L, this ) );
	}
	else if ( type=="electron_impact_dissociation" ) {
	    reactions.push_back( new ElectronImpactDissociation( L, this ) );
	}
	else if ( type=="radiative_transition" ) {
	    reactions.push_back( new RadiativeTransition( L, this ) );
    	}
    	else {
    	    ostringstream oss;
    	    oss << "QSSDiatomicRadiator::create_reactions()" << endl
    	        << "Reaction type: " << type << " not recognised." << endl;
    	    input_error(oss);
    	}
	lua_pop(L, 1 );		// pop reaction 'i'
    }
    lua_pop(L, 1);	// pop reactions table
    
    return reactions.size();
}

void
QSSDiatomicRadiator::
calculate_n_e( Gas_data &Q )
{
    /* Linear approach, solve system via Gaussian Elimination */
    
    // Use Boltzmann populations if temperature is low
    if ( Q.T.back() < T_lower ) {
    	for ( int ilev=0; ilev<nlevs; ++ilev ) {
	    elevs[ilev]->set_N( eval_Boltzmann_population_for_level(Q,ilev) );
    	}
    	return;
    }
    
    // 0. Reset Jacobian matrix and source and solution vectors
    for ( size_t i=0; i<noneq_elevs.size(); ++i ) {
	C[i] = 0.0;
	y_out[i] = 0.0;
	for ( size_t j=0; j<noneq_elevs.size(); ++j ) {
	    dGdy->set( i, j, 0.0 );
	}
    }
    
    // 1.  Population summations, first matrix row
    for ( size_t ne_ilev=0; ne_ilev<noneq_elevs.size(); ++ne_ilev ) {
    	double bf_acc = 0.0;	// accumulated boltzmann fractions
    	for ( size_t eq_ilev=0; eq_ilev<noneq_elevs[ne_ilev]->eq_elevs.size(); ++eq_ilev ) {
    	    bf_acc += noneq_elevs[ne_ilev]->eq_elevs[eq_ilev]->get_Q_int() / noneq_elevs[ne_ilev]->elev->get_Q_int();
    	}
    	double tmp = dGdy->get(0,ne_ilev) + 1.0 + bf_acc;
    	dGdy->set(0,ne_ilev,tmp);
    }
    
    // 2. Contributions from reactions
    for ( size_t ir=0; ir<reactions.size(); ++ir )
	reactions[ir]->add_jacobian_contributions( Q, *dGdy );
    
    // 3. Construct source vector
    C[0] = Q.massf[isp] * Q.rho / m_w / 1.0e6;		// Convert moles/m**3 to moles/cm**3
    // NOTE: some reactions such as dissociation have terms that will not be a function of the unkown
    //       populations, therefore they need to go in the source vector for this method
    for ( size_t ir=0; ir<reactions.size(); ++ir )
    	reactions[ir]->add_source_vector_contributions( Q, C );
    
    // 4. Solve the system
    if( dGdy->gaussian_elimination( y_out, C ) ) {
        cout << "QSSDiatomicRadiator::calculate_n_e()" << endl
             << "Gaussian elimination failed for QSSDiatomicRadiator: " << name << endl
             << "The gas-state was: " << endl
             << "Q.T[0] = " << Q.T[0] << endl
             << "Q.p = " << Q.p << endl
             << "Q.p_e = " << Q.p_e << endl;
        exit( FAILURE );
    }
   
    // 5.  Map results back onto radiator
    for ( size_t ne_ilev=0; ne_ilev<noneq_elevs.size(); ++ne_ilev ) {
    	// 5a. Firstly noneq levels
    	double N_ne_ilev = y_out[ne_ilev] * 1.0e6 * RC_Na;	// Convert moles/cm**3 -> particles/m**3
    	if ( N_ne_ilev < 0.0 ) {
    	    cout << "QSSDiatomicRadiator::calculate_n_e()" << endl
    	         << name << ".N[" << noneq_elevs[ne_ilev]->elev->i << "] = " << N_ne_ilev << endl
    	         << "Bailing out!" << endl;
    	    exit( FAILURE );
    	}
    	noneq_elevs[ne_ilev]->elev->set_N( N_ne_ilev );
	// 5b. Now equilibriated levels
	for ( size_t eq_ilev=0; eq_ilev<noneq_elevs[ne_ilev]->eq_elevs.size(); ++eq_ilev ) {
	    double bf = noneq_elevs[ne_ilev]->eq_elevs[eq_ilev]->get_Q_int() / noneq_elevs[ne_ilev]->elev->get_Q_int();
	    noneq_elevs[ne_ilev]->eq_elevs[eq_ilev]->set_N( N_ne_ilev * bf );
	}
    }
    
#   if ENFORCE_DIATOMIC_BOLTZ_QSS_LIMIT
    // 6. Enforce the Boltzmann limit
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
    	double N_boltz = eval_Boltzmann_population_for_level(Q,ilev);
    	// Boltzmann as the upper limit
    	if (  !finite(elevs[ilev]->get_N()) || elevs[ilev]->get_N() > N_boltz ) elevs[ilev]->set_N(N_boltz);
    	// Check for finiteness
    	/*
    	if ( !finite( elevs[ilev]->get_N() ) ) {
    	    cout << "QSSDiatomicRadiator::calculate_n_e()" << endl
    	         << "Radiator: " << name << " QSS calculation failed." << endl;
    	    cout << "Source vector: \n"; print_valarray(C);
    	    cout << "y_out: \n"; print_valarray(y_out);
    	    cout << "Bailing out!" << endl;
    	    exit( FAILURE );
    	}
    	*/
    }
#   endif
    
#   if DEBUG_RAD > 0
    // 7. Check that the total species density has been conserved
    double N_total = 0.0;
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
	N_total += elevs[ilev]->get_N();
	cout << "N_el[" << ilev << "] = " << elevs[ilev]->get_N() << endl;
    }
    cout << "N_total(QSS) = " << N_total << ", N_total(CFD) = " << Q.massf[isp] * Q.rho / m_w * RC_Na << endl;
#   endif
}

/*************************** Helper functions ******************************/

int get_diatomic_transition_type( DiatomicElecLev * elev_u, DiatomicElecLev * elev_l )
{
    // Will usually be type 0 (unsplit parallel or perpendicular transition)
    int transition_type = HUND_A;
    
    // Remove the following clause once we are sure the input file is correct
    int delta_L = elev_u->get_lambda()-elev_l->get_lambda();
    int delta_S = elev_u->get_spin()-elev_l->get_spin();
    if ( abs(delta_L)>1 || delta_S!=0 ) {
    	cout << "get_diatomic_transition_type()" << endl
	     << "delta_L = " << delta_L << ", delta_S = " << delta_S << endl
	     << "This is a forbidden transition, need to check input file." << endl;
	cout << "Bailing out!" << endl;
	exit( BAD_INPUT_ERROR );
    }
    
    if ( elev_u->has_gamma_Vs() && elev_u->get_lambda()==0 && elev_u->get_spin()==3
      && elev_l->has_gamma_Vs() && elev_l->get_lambda()==0 && elev_l->get_spin()==3 ) {
    	// 3Sigma<->3Sigma transition with spin splitting
    	transition_type = HUND_B_TRIPLET;
    	ostringstream oss;
    	oss << "get_diatomic_transition_type()" << endl
    	    << "Triplet transitions belonging to Hund's case (b) are not yet available." << endl
    	    << "Remove the spin splitting data and this transition will be modeled by" << endl
    	    << "Hund's case (a) with an effective line-center" << endl;
    	input_error(oss);
    }
    else if ( elev_u->get_lambda()==0 && elev_u->get_spin()==3
           && elev_l->get_lambda()==0 && elev_l->get_spin()==3 ) {
	// Parallel (Sigma - Sigma) triplet
	transition_type = SIGMA_TRIPLET;
    }
    else if ( elev_u->has_gamma_Vs() && elev_u->get_lambda()==0 && elev_u->get_spin()==2
           && elev_l->has_gamma_Vs() && elev_l->get_lambda()==0 && elev_l->get_spin()==2 ) {
    	// 2Sigma<->2Sigma transition with spin splitting (eg. CN Violet)
    	transition_type = HUND_B_DOUBLET;
    }
#   if HUND_AB_DOUBLET_HLF_METHOD==0
    else if ( ( elev_u->get_spin()==2 && elev_l->get_spin()==2 ) &&
              ( ( elev_u->get_lambda()==1 && elev_l->get_lambda()==0 ) ||
              	( elev_u->get_lambda()==0 && elev_l->get_lambda()==1 ) ) ) {
    	// 2Sigma <-> 2Pi transition
    	transition_type = HUND_AB_DOUBLET;
    }
#   else
    else if ( ( elev_u->get_spin()==2 && elev_l->get_spin()==2 ) &&
              ( elev_u->get_lambda()!=0 || elev_l->get_lambda()!=0 ) ) {
    	// Doublet intermediate (a)-(b) case
    	transition_type = HUND_AB_DOUBLET;
    }
#   endif
#   if MODEL_HUND_AB_TRIPLETS
    else if ( ( elev_u->get_spin()==3 && elev_l->get_spin()==3 ) &&
               ( ( elev_u->get_lambda()!=0 || elev_l->get_lambda()!=0 ) ) ) {
    	// Triplet intermediate (a)-(b) case
    	transition_type = HUND_AB_TRIPLET;
    }
#   endif
	       
    return transition_type;
}

