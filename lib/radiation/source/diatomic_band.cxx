/** \file diatomic_band.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 11-Mar-08: Original version
 *           11-Aug-09: Photaura version
 *  \brief Definitions for Diatomic_band class.
 *
 **/

#include <stdlib.h>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "../../util/source/useful.h"

#include "radiation_constants.hh"
#include "diatomic_band.hh"
#include "diatomic_radiator.hh"

using namespace std;

// This constant is used for converting transition moments to probabilities
const double A_ul_const = ( 64.0 / 3.0 ) * pow( M_PI, 4 ) / RC_h * pow( RC_a0*RC_e, 2 );

DiatomicLine::
DiatomicLine( double nu_00, double E_u, double m_w, double I, double sigma_nm )
: nu_00( nu_00 ), E_u( E_u ), m_w( m_w ), I( I ), sigma_nm( sigma_nm ) {}

DiatomicLine::
~DiatomicLine() {}

void
DiatomicLine::
initialize(  double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    // 1. Lorentz half-width
    gamma_L = calculate_Lorentz_width( T, Te, p, N_hvy, N_elecs, mw_av );
    
    // 2. Doppler half-width
    gamma_D = calculate_Doppler_width( T );
    
    // 3. Voigt half-width
    gamma_V = calculate_Voigt_width();
    
#   if DEBUG_RAD > 1
    cout << "gamma_L = " << nu2lambda(nu_00) - nu2lambda(gamma_L+nu_00) << " nm" << endl
	 << "gamma_D = " << nu2lambda(nu_00) - nu2lambda(gamma_D+nu_00) << " nm" << endl
	 << "gamma_V = " << nu2lambda(nu_00) - nu2lambda(gamma_V+nu_00) << " nm" << endl;
#   endif

    return;
}

double
DiatomicLine::
calculate_Lorentz_width( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    // 1. Collision broadening
    double gamma_C = calculate_collision_width( T, p, N_hvy, mw_av );
    
    // 2. Stark broadening
    double gamma_S = calculate_Stark_width( Te, N_elecs );
    
    // 3. Natural line broadening
    double gamma_N = calculate_natural_width();
    
    // Sum the contributions
    return gamma_C + gamma_S + gamma_N;
}

double
DiatomicLine::
calculate_collision_width( double T, double p, double N_hvy, double mw_av )
{
    // Two different methods are available (give very similiar results, it seems):
    
#   if DIATOMIC_COLLISION_WIDTH_METHOD == 0
    // Expression A-19 from Johnston 2006 CN paper, with reference:
    // Colket, M.B., Sectroscopic Absorption Model for CN(X2Sig-B2Sig): Comparison of Experiments and Theory, 
    // Journal of Quantitative Spectroscopy and Radiative Transfer, Vol. 31, 1984, pp. 7-13
    
    double lambda_cl = nu2lambda( nu_00 ) * 1.0e-7;		// cm
    double p_torr = p / 133.33;
    
    double gamma_C_cm = 0.5*lambda_cl*lambda_cl*0.04687*sigma_nm*sigma_nm*
			    p_torr * sqrt( 2.0 / ( T * m_w * 1.0e3 ) );
    
    // Convert delta cm into delta Hz
    double gamma_C = gamma_C_cm / lambda_cl * nu_00;
    
    return gamma_C;
#   elif DIATOMIC_COLLISION_WIDTH_METHOD == 1
    // Spradian07 method, pressure broadening by non-resonant collisions (Traving)
    
    double lambda_cl = nu2lambda( nu_00 ) * 1.0e1;		// Angstroms
    
    double gamma_C_Ang = lambda_cl*lambda_cl*5.85e-30*sqrt(2.0/mw_av)*N_hvy*sqrt(T);
    
    // Convert delta Ang into delta Hz
    double gamma_C = gamma_C_Ang / lambda_cl * nu_00;
    
    return gamma_C;
#   endif
}

double
DiatomicLine::
calculate_Stark_width( double Te, double N_elecs )
{
    /* Implementing the default method from Spradian07:
       - Use default exponential term n of 0.33
       - various constC methods:
          0 : Johnston 2006 curve fit
          1 : Cowley 1971 
          2 : Arnold 1979 curve fit developed for Si and C atoms 
    */
#   if DIATOMIC_STARK_METHOD==0
    double constA = 1.69e10; double constB = 2.623;
#   elif DIATOMIC_STARK_METHOD==1
    double constA = 9.27e07; double constB = 2.000;
#   elif DIATOMIC_STARK_METHOD==2
    double constA = 4.20e07; double constB = 2.000;
#   endif
    // NOTE: - converting energies from J -> cm^-1
    //       - applying fabs() around delta_E as E_u can be > E_ionise
    double constC =  0.5 * constA / pow( fabs(I - E_u)/(RC_c*RC_h_SI), constB );
    double gamma_S0 = constC * RC_c;
    double gamma_S = gamma_S0 * pow( (Te / 1.0e4), 0.33 ) * ( N_elecs / 1.0e16 );
    
    return gamma_S;
}

double
DiatomicLine::
calculate_natural_width()
{
    /* Method A: use constant value as used by Spradian07.
       Method B: use classical radiation theory expression, p.190 in Thorne et al "Spectrophysics..." */
       
#   if DIATOMIC_NAUTRAL_WIDTH_METHOD==0
    return 1.18e-5 * nu_00 / nu2lambda(nu_00);
#   elif DIATOMIC_NAUTRAL_WIDTH_METHOD==1
    return ( 2.0 * M_PI * RC_e_SI*RC_e_SI * nu_00*nu_00 ) /
	          ( 3.0 * RC_eps0_SI * RC_m_SI * RC_c_SI*RC_c_SI*RC_c_SI );
#   endif
}

double
DiatomicLine::
calculate_Doppler_width( double T )
{
    /* Calculate the doppler half-width at half-intensity for a particular line */
    double m_s = m_w * 1000.0 / RC_Na;		/* mass per particle in grams */
    
    /* NOTE: Johnston uses T_e, but T should govern Maxwell distribution... */
    return (nu_00 / RC_c) * sqrt ( (2.0 * RC_k * T * log(2.0)) / m_s);
}

double
DiatomicLine::
calculate_Voigt_width()
{
    double d = ( gamma_L - gamma_D ) / (gamma_L + gamma_D);
    double beta = 0.023665 * exp ( 0.6 * d) + 0.0418 * exp ( -1.9 * d);
    double alpha = 0.18121;
    double R_d = 1.0 - alpha * ( 1.0 - d * d ) - beta * sin(M_PI * d);
    
    return R_d * ( gamma_L + gamma_D);
}

double
DiatomicLine::
get_lorentz_point( double delta_nu )
{
    /* Line shape as a function of delta_nu for a Lorentz profile */
    double bp_nu;
    
    bp_nu = gamma_L / (M_PI * ( delta_nu * delta_nu ) + gamma_L * gamma_L );
    
    return bp_nu;
}

double
DiatomicLine::
get_doppler_point( double delta_nu )
{
    /* Line shape as a function of delta_nu for a Doppler profile */
    double bd_nu;
    
    bd_nu = sqrt (log(2.0) / M_PI ) / gamma_D * exp ( -log(2.0) * pow( delta_nu/gamma_D, 2.0) );
    
    return bd_nu;
}

double
DiatomicLine::
get_voigt_point( double delta_nu )
{
    // Ref: Whiting (1968) JQRST Vol. 8 pp 1379-1384
    double R_l = delta_nu / ( 2.0 * gamma_V );
    double R_d = gamma_L / gamma_V;
    
#   if DIATOMIC_VOIGT_PROFILE_METHOD == 0
    // Accurate expression
    double tmpA = ( 1.0 - R_d ) * exp( -2.772 * R_l * R_l ) + R_d / ( 1.0 + 4.0 * R_l * R_l );
    double tmpB = 0.016 * ( 1.0 - R_d ) * R_d * ( exp( -0.4 * pow( R_l, 2.25 ) ) - 10.0 / ( 10.0 + pow( R_l, 2.25 ) ) );
    double tmpC = 2.0 * gamma_V * (1.065 + 0.447 * R_d + 0.058 * R_d * R_d);
    
    double b_nu = ( tmpA + tmpB ) / tmpC;
    
#   elif DIATOMIC_VOIGT_PROFILE_METHOD == 1
    // Approximate expression
    double tmpA = ( 1.0 - R_d ) * exp( -2.772 * R_l * R_l ) + R_d / ( 1.0 + 4.0 * R_l * R_l );
    double tmpC = 2.0 * gamma_V * (1.065 + 0.447 * R_d + 0.058 * R_d * R_d);
    
    double b_nu = tmpA / tmpC;
    
#   endif
    
    return b_nu;
}

void
DiatomicLine::
calculate_spectrum( CoeffSpectra &X )
{
    // 1. Calculate specific spectral range if we can, otherwise loop over all frequencies
    int inu_start = 0;
    int inu_end = int( X.nu.size() );
#   if DIATOMIC_LIMITED_LINE_EXTENT
    double nu_lower = nu_ul - double(DIATOMIC_LINE_EXTENT) * gamma_V;
    double nu_upper = nu_ul + double(DIATOMIC_LINE_EXTENT) * gamma_V;
    inu_start = get_nu_index(X.nu,nu_lower,X.adaptive) + 1;
    inu_end = get_nu_index(X.nu,nu_upper,X.adaptive) + 1;
#   endif
	
    // 2. Loop over predetermined frequency range,
    //    compute j_nu and kappa_nu
    for ( int inu=inu_start; inu<inu_end; inu++ ) {
	double nu = X.nu[inu];
	double delta_nu = fabs( nu - nu_ul );
	double b_nu = get_voigt_point(delta_nu);
	// cout << "inu = " << inu << ", b_nu = " << b_nu << endl;
	X.j_nu[inu] += j_ul * b_nu;
	X.kappa_nu[inu] += kappa_lu * b_nu;
    }
    
    return;
}

string
DiatomicLine::
line_width_string(  double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    ostringstream ost;
    ost << setprecision(12) << scientific << showpoint;
    
    double lambda_nm_ul = nu2lambda(nu_00);
    double gamma_C_nm = lambda_nm_ul / nu_00 * calculate_collision_width(T,p,N_hvy,mw_av);
    double gamma_S_nm = lambda_nm_ul / nu_00 * calculate_Stark_width(Te,N_elecs);
    double gamma_N_nm = lambda_nm_ul / nu_00 * calculate_natural_width();
    double gamma_D_nm = lambda_nm_ul / nu_00 * calculate_Doppler_width(T);
    ost << setw(20) << lambda_nm_ul << setw(20) << gamma_D_nm << setw(20)
	<< setw(20) << gamma_C_nm   << setw(20) << gamma_S_nm << setw(20) << gamma_N_nm << endl;
    
    return ost.str();
}

DiatomicBand::
DiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, double Re_vib )
 : Vu( Vu ), Vl( Vl ), elev_u( elev_u ), elev_l( elev_l ), Re_vib( Re_vib )
{
    // Set the electronic level related data
    double E_u = elev_u->calculate_E_vib(Vu) + elev_u->E; 
    double E_l = elev_l->calculate_E_vib(Vl) + elev_l->E;
    
    // 0. Transition type
    transition_type = get_diatomic_transition_type( elev_u, elev_l );
    
    // 1. set the reverse_band flag
    if ( E_u < E_l ) {
	reverse_band = true;
#       if DEBUG_RAD > 0
    	cout << "DiatomicBand::set_elev_data()" << endl
    	     << "[Vu=" << Vu << ", Vl=" << Vl << "] is a reverse band as delta_E = " << E_u - E_l << endl;
#       endif
    }
    else {
    	reverse_band = false;
    }
    
    // 2. Spin splitting
    if ( elev_u->get_spin()==1 ) tS = 0;
    else if ( elev_u->get_spin()==2 ) tS = 1;
    else if ( elev_u->get_spin()==3 ) tS = 2;
    else {
    	cout << "DiatomicBand::DiatomicBand()" << endl
    	    << "2S+1 = " << elev_u->get_spin() << " was not expected." << endl
    	    << "Bailing out!" << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    // 3. Pull out lambda and J_max values from the electronic levels
    tJu_max = 2*elev_u->get_J_max(Vu);
    tJl_max = 2*elev_l->get_J_max(Vl);
    
    lambda_u = elev_u->get_lambda();
    lambda_l = elev_l->get_lambda();

    // 4. Calculate characteristic frequencies
    // 4a. 00 transition frequency
    nu_00 = ( E_u - E_l ) / RC_h_SI;
    
    // 4b. Band-head frequency (See Golden (1966) eq 34.)
    double Bv_u  = elev_u->calculate_B_v( Vu );
    double Bv_l  = elev_l->calculate_B_v( Vl );
    double Bv_term = 0.25 * pow( (Bv_u + Bv_l), 2) / ( Bv_u - Bv_l );
    nu_bh = ( E_u - E_l - Bv_term ) / RC_h_SI;
    
    // 4c. Band-average frequency 
    nu_av = calculate_average_frequency();
    
    // 5. Calculate the band average transition probability
    A_ul_av = A_ul_const * pow( Re_vib, 2 ) * pow( nu_av / RC_c, 3 );
}

DiatomicBand::
~DiatomicBand() {}

double DiatomicBand::calculate_average_frequency()
{
    // NOTE: this function should replicate the behaviour of 
    //       DiatomicSystem::calculate_average_frequency()
    double E_u = elev_u->calculate_E_vib(Vu) + elev_u->E; 
    double E_l = elev_l->calculate_E_vib(Vl) + elev_l->E;
    
    double nu_av = 0.0;
#   if BAND_AVERAGE_FREQUENCY_METHOD==0
    // NOTE: count should be equal to 2 Ju + 1 for Sigma transitions and 3 Ju + 1 otherwise
    int count = 0;
    double nu_acc = 0.0;
    
    // J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components -> CHECKME
    int tJu_min = abs(2*lambda_u);
    for ( int tJu = tJu_min; tJu <= tJu_max; tJu += 2 ) {
	// J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components -> CHECKME
	int tJl_min = abs(2*lambda_l);
	for ( int delta_J = -1; delta_J <=1; delta_J+=1 ) {
	    int tJl = tJu + 2*delta_J;
	    /* Check if final rotational state is permitted */
	    if ( tJl>tJl_max || tJl<tJl_min ) continue;
	    /* 0-0 transition is universally prohibited */
	    if ( tJu==0 && tJl==0 ) continue;
	    /* delta_J = 0 transitions are forbidden for Sigma-Sigma transitions */
	    if ( delta_J == 0 && ( lambda_l == 0 && lambda_u == 0 ) ) continue;
	    /* If we have got to here this is an allowed transition */
	    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundA(Vu,tJu) - elev_l->calculate_E_rot_HundA(Vl,tJl) ) / RC_h_SI;
	    nu_acc += fabs(nu_cl);
	    ++count;
	}
    }
    nu_av = nu_acc / double( count );
#   elif BAND_AVERAGE_FREQUENCY_METHOD==1
    // See Lin Hartung Chambers thesis p 32 for an explanation of this expression
    nu_00 = ( E_u - E_l ) / RC_h_SI;
    double Bv_u  = elev_u->calculate_B_v( Vu );
    double Bv_l  = elev_l->calculate_B_v( Vl );
    int J_max_l = (int) elev_l->get_J_max(Vl);
    Bv_u /= RC_h_SI;	// convert J -> Hz
    Bv_l /= RC_h_SI;
    if ( elev_l->get_lambda() == 0 && elev_u->get_lambda() == 0 ) {
    	// Sigma -> Sigma transition
    	double tmpA = 2.0 * ( Bv_u ) * double ( J_max_l + 1) / double ( 2 * J_max_l + 1 );
    	double tmpB = 2.0 * ( Bv_u - Bv_l ) / double ( 2 * J_max_l + 1 );
    	double tmpC = double( J_max_l * ( J_max_l + 1 ) * ( J_max_l + 2 ) ) / 3.0;
    	nu_av = nu_00 + tmpA + tmpB * tmpC;
    }
    else {
    	// Not a Sigma -> Sigma transition
    	double tmpA = 2.0 * ( Bv_u ) * double ( J_max_l + 1) / double ( 3 * J_max_l + 1 );
    	double tmpB = 3.0 * ( Bv_u - Bv_l ) / double ( 3 * J_max_l + 1 );
    	double tmpC = double( J_max_l * ( J_max_l + 1 ) * ( J_max_l + 2 ) ) / 1.0;
    	nu_av = nu_00 + tmpA + tmpB * tmpC;
    }
#   else
    // just use the zero-zero frequency
    nu_av = ( E_u - E_l ) / RC_h_SI;
#   endif
    
    return fabs(nu_av);
}

double DiatomicBand::calculate_j_ul( double N_u )
{
    return N_u * RC_h_SI * nu_av * A_ul_av / ( 4.0 * M_PI );
}

LBLDiatomicBand::
LBLDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		 double Re_vib, double m_w, double I, double sigma_nm  )
 : DiatomicBand( Vu, Vl, elev_u, elev_l, Re_vib )
{
    line = new DiatomicLine( nu_00, elev_u->get_E(), m_w, I, sigma_nm );
}

LBLDiatomicBand::
~LBLDiatomicBand()
{
    delete line;
}

void
LBLDiatomicBand::
initialize(  double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    line->initialize( T, Te, p, N_hvy, N_elecs, mw_av );

    return;
}

std::string
LBLDiatomicBand::
line_width_string(  double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    return line->line_width_string( T, Te, p, N_hvy, N_elecs, mw_av );
}

HundALBLDiatomicBand::
HundALBLDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		 double Re_vib, double m_w, double I, double sigma_nm  )
: LBLDiatomicBand( Vu, Vl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm )
{}

void
HundALBLDiatomicBand::
write_lines_to_file( string fname )
{
    /* Setup the output file */
    ofstream ofile;
    ofile << setprecision(6) << scientific << showpoint;
    ofile.open(fname.c_str());
    ofile << "# " << fname << endl
	  << "# Column 1: J_u" << endl
	  << "# Column 2: J_l" << endl
	  << "# Column 3: S" << endl
	  << "# Column 4: lambda_cl" << endl
	  << "# Column 5: A_ul" << endl;
    
    double S_sum = 0.0;
    int tJu_sum = 10;
    
    // J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
    int tJu_min = abs(2*lambda_u);
    for ( int tJu = tJu_min; tJu <= tJu_max; tJu += 2 ) {
	// J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
	int tJl_min = abs(2*lambda_l);
	for ( int delta_J = -1; delta_J <=1; delta_J+=1 ) {
	    int tJl = tJu + 2*delta_J;
	    /* Check if final rotational state is permitted */
	    if ( tJl>tJl_max || tJl<tJl_min ) continue;
	    /* 0-0 transition is universally prohibited */
	    if ( tJu==0 && tJl==0 ) continue;
	    /* delta_J = 0 transitions are forbidden for Sigma-Sigma transitions */
	    if ( delta_J == 0 && ( lambda_l == 0 && lambda_u == 0 ) ) continue;
	    /* If we have got to here this is an allowed transition */
	    double S = get_HLF( tJu, tJl );
	    if ( tJu==tJu_sum ) S_sum += S;
	    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundA(Vu,tJu) - elev_l->calculate_E_rot_HundA(Vl,tJl) ) / RC_h_SI;
	    double lambda_cl = RC_c/nu_cl;
	    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );
	    ofile << setw(20) << 0.5*double(tJu)
		  << setw(20) << 0.5*double(tJl)
		  << setw(20) << S
		  << setw(20) << lambda_cl*1.0e7 
		  << setw(20) << A_ul
		  << endl;
	}
    }
    
    cout << "S_sum = " << S_sum << ", double(tJu+1) = " << double(tJu_sum+1) << endl;
    
    ofile.close();
    
    return;
}

void
HundALBLDiatomicBand::
calculate_spectrum( CoeffSpectra &X, double Tv, double Tr )
{
    double j_ul_00 = -1.0;
    // J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
    int tJu_min = abs(2*lambda_u);
    for ( int tJu = tJu_min; tJu <= tJu_max; tJu += 2 ) {
	// J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
	int tJl_min = abs(2*lambda_l);
	for ( int delta_J = -1; delta_J <=1; delta_J+=1 ) {
	    int tJl = tJu + 2*delta_J;
	    /* Check if final rotational state is permitted */
	    if ( tJl>tJl_max || tJl<tJl_min ) continue;
	    /* 0-0 transition is universally prohibited */
	    if ( tJu==0 && tJl==0 ) continue;
	    /* delta_J = 0 transitions are forbidden for Sigma-Sigma transitions */
	    if ( delta_J == 0 && ( lambda_l == 0 && lambda_u == 0 ) ) continue;
	    /* If we have got to here this is an allowed transition */
	    calculate_rot_line_spectrum(X,Tv,Tr,tJu,tJl,j_ul_00);
	}
    }
    
    return;
}

void
HundALBLDiatomicBand::
calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, double &j_ul_00 )
{
    // 1. Calculate integrated emission and absorption coefficients
    // 1a. calculate Einstein coefficients
    double S = get_HLF( tJu, tJl );
    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundA(Vu,tJu) - elev_l->calculate_E_rot_HundA(Vl,tJl) ) / RC_h_SI;
    // NOTE: lambda must have cgs units
    double lambda_cl = RC_c/nu_cl;
    // NOTE: Spin/lambda splitting degeneracy should be included in Honl-London factor S
    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );

    // only degeneracy ratio is used, so don't worry about correcting for spin splitting
    int g_u = elev_u->get_g() * ( tJu + 1 );
    int g_l = elev_l->get_g() * ( tJl + 1 );
    
    double B_lu = A_ul * double( g_u ) / double( g_l ) * RC_c_SI * RC_c_SI 
		/ ( 8.0 * M_PI * RC_h_SI ) / ( nu_cl*nu_cl*nu_cl );
		
    double B_ul = B_lu * double( g_l )  / double ( g_u );
    
    // 1b. Upper and lower rotational state number densities
    double n_u = elev_u->calculate_N_rot_HundA(Tv,Tr,Vu,tJu,true);
    double n_l = elev_l->calculate_N_rot_HundA(Tv,Tr,Vl,tJl);
    
    // 1c. Combine the pieces to get integrated coefficients
    double j_ul = n_u * A_ul * RC_h_SI * nu_cl / ( 4.0 * M_PI );
    double kappa_lu = ( n_l * B_lu - n_u * B_ul ) * RC_h_SI * nu_cl;
    
    if ( j_ul_00 < 0.0 ) j_ul_00 = j_ul;
    else if ( j_ul < F_ROT_LINE_LIMIT * j_ul_00 ) return;
    
    // 2. Give coefficients to the line so that the line profile can be calculated
    line->nu_ul = nu_cl;
    line->j_ul = j_ul;
    line->kappa_lu = kappa_lu;
    line->calculate_spectrum( X );
    
    return;
}

double
HundALBLDiatomicBand::
get_HLF( int tJu, int tJl)
{
    int delta_J = (tJl - tJu)/2;
    int delta_L = lambda_l - lambda_u;
    double S = 0.0;
    double Ju = double(tJu)/2.0;
    
    // Calculate Honl-London factor for this transition
    // See JQRST v9 pp 775-798 Arnold et al 1969 for HLF expressions used here
    
    if ( delta_L==0 ) {
    	// 1. A parallel transition
    	if ( lambda_u==0 ) {
    	    // Use Sigma - Sigma expression
    	    if      ( delta_J == -1 ) S = Ju;
    	    else if ( delta_J ==  0 ) S = 0.0;
    	    else  if ( delta_J == 1 ) S = Ju + 1.0;
    	}
    	else {
	// Use full HLF singlet expressions (See Huber and Herzberg p 208)
	    if      ( delta_J == -1 ) S = ( Ju + double(lambda_u) )
					* ( Ju - double(lambda_u) )
					/ Ju;
	    else if ( delta_J ==  0 ) S = ( 2.0 * Ju + 1.0 )
					* double(lambda_u) * double(lambda_u)
					/ ( Ju * ( Ju + 1.0 ) );
	    else if ( delta_J ==  1 ) S = ( Ju + 1.0 + double(lambda_u) )
					* ( Ju + 1.0 - double(lambda_u) )
					/ ( Ju + 1.0 );
	}
    }
    else if ( abs(delta_L)==1 ) {
    	// 2. A perpendicular transition
        // delta_L = -1 and delta_L = 1 have slightly different expressions that are distinguished by sign of delta_L
        double sign = double(delta_L);
        // R, Q and P branches respectively, ignoring spin splitting
        // NOTE: check these expressions as they are different from Huber & Herzberg p. 208 by a factor of 2
        if      ( delta_J == -1 ) S = ( Ju + sign * double(lambda_u) ) 
                                    * ( Ju - 1.0 + sign * double(lambda_u) )
                                    / ( 2.0 * Ju ); 
        else if ( delta_J ==  0 ) S = ( Ju + sign * double(lambda_u) ) 
                                    * ( Ju + 1.0 - sign * double(lambda_u) )
                                    * ( 2.0 * Ju + 1.0 )
                                    / ( 2.0 * Ju * ( Ju + 1.0 ) ); 
        else if ( delta_J ==  1 ) S = ( Ju + 1.0 - sign * double(lambda_u) ) 
                                    * ( Ju + 2.0 - sign * double(lambda_u) )
                                    / ( 2.0 * ( Ju + 1.0 ) );
    }
    
    return S;
}

// ---------------

SigmaTripletLBLDiatomicBand::
SigmaTripletLBLDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		 double Re_vib, double m_w, double I, double sigma_nm  )
: LBLDiatomicBand( Vu, Vl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm )
{}


void
SigmaTripletLBLDiatomicBand::
write_lines_to_file( string fname )
{
    /* Setup the output file */
    ofstream ofile;
    ofile << setprecision(6) << scientific << showpoint;
    ofile.open(fname.c_str());
    ofile << "# " << fname << endl
	  << "# Column 1: J_u" << endl
	  << "# Column 2: J_l" << endl
	  << "# Column 3: S" << endl
	  << "# Column 4: lambda_cl" << endl
	  << "# Column 5: A_ul" << endl;
	  
    double S_sum = 0.0;
    int tJu_sum = 10;
    
    // J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
    int tJu_min = abs(2*lambda_u);
    for ( int tJu = tJu_min; tJu <= tJu_max; tJu += 2 ) {
	// J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
	int tJl_min = abs(2*lambda_l);
	for ( int delta_J = -1; delta_J <=1; delta_J+=1 ) {
	    int tJl = tJu + 2*delta_J;
	    /* Check if final rotational state is permitted */
	    if ( tJl>tJl_max || tJl<tJl_min ) continue;
	    /* 0-0 transition is universally prohibited */
	    if ( tJu==0 && tJl==0 ) continue;
	    /* delta_J = 0 transitions are forbidden for Sigma-Sigma transitions */
	    if ( delta_J == 0 && ( lambda_l == 0 && lambda_u == 0 ) ) continue;
	    /* If we have got to here this is an allowed transition */
	    double S = get_HLF( tJu, tJl );
	    if ( tJu==tJu_sum ) S_sum += S;
	    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_SigmaTriplet(Vu,tJu) - elev_l->calculate_E_rot_SigmaTriplet(Vl,tJl) ) / RC_h_SI;
	    double lambda_cl = RC_c/nu_cl;
	    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );
	    ofile << setw(20) << 0.5*tJu 
		  << setw(20) << 0.5*tJl
		  << setw(20) << S
		  << setw(20) << lambda_cl*1.0e7 
		  << setw(20) << A_ul
		  << endl;
	}
    }
    
    cout << "S_sum = " << S_sum << ", double(tJu+1) = " << double(tJu_sum+1) << endl;
    
    ofile.close();
    
    return;
}

void
SigmaTripletLBLDiatomicBand::
calculate_spectrum( CoeffSpectra &X, double Tv, double Tr )
{
    double j_ul_00 = -1.0;
    // J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
    int tJu_min = abs(2*lambda_u);
    for ( int tJu = tJu_min; tJu <= tJu_max; tJu += 2 ) {
	// J can't be smaller than the (absolute) sum of nuclear and spin (=0 here) components
	int tJl_min = abs(2*lambda_l);
	// NOTE: incrementing delta_J by two to skip over delta_J = 0 which is prohibited for Sigma transitions
	for ( int delta_J = -1; delta_J <=1; delta_J+=2 ) {
	    int tJl = tJu + 2*delta_J;
	    /* Check if final rotational state is permitted */
	    if ( tJl>tJl_max || tJl<tJl_min ) continue;
	    /* 0-0 transition is universally prohibited */
	    if ( tJu==0 && tJl==0 ) continue;
	    /* If we have got to here this is an allowed transition */
	    calculate_rot_line_spectrum(X,Tv,Tr,tJu,tJl,j_ul_00);
	}
    }
    
    return;
}

void
SigmaTripletLBLDiatomicBand::
calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, double &j_ul_00 )
{
    // 1. Calculate integrated emission and absorption coefficients
    // 1a. calculate Einstein coefficients
    double S = get_HLF( tJu, tJl );
    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_SigmaTriplet(Vu,tJu) - elev_l->calculate_E_rot_SigmaTriplet(Vl,tJl) ) / RC_h_SI;
    // NOTE: lambda must have cgs units
    double lambda_cl = RC_c/nu_cl;
    // NOTE: Spin/lambda splitting degeneracy should be included in Honl-London factor S
    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );

    // only degeneracy ratio is used, so don't worry about correcting for spin splitting
    // as 2S+1 must be the same for the upper and lower levels
    int g_u = elev_u->get_g() * ( tJu + 1 );
    int g_l = elev_l->get_g() * ( tJl + 1 );
    
    double B_lu = A_ul * double( g_u ) / double( g_l ) * RC_c_SI * RC_c_SI 
		/ ( 8.0 * M_PI * RC_h_SI ) / ( nu_cl*nu_cl*nu_cl );
		
    double B_ul = B_lu * double( g_l )  / double ( g_u );
    
    // 1b. Upper and lower rotational state number densities
    double n_u = elev_u->calculate_N_rot_SigmaTriplet(Tv,Tr,Vu,tJu,true);
    double n_l = elev_l->calculate_N_rot_SigmaTriplet(Tv,Tr,Vl,tJl);
    
    // 1c. Combine the pieces to get integrated coefficients
    double j_ul = n_u * A_ul * RC_h_SI * nu_cl / ( 4.0 * M_PI );
    double kappa_lu = ( n_l * B_lu - n_u * B_ul ) * RC_h_SI * nu_cl;
    
    if ( j_ul_00 < 0.0 ) j_ul_00 = j_ul;
    else if ( j_ul < F_ROT_LINE_LIMIT * j_ul_00 ) return;
    
    // 2. Give coefficients to the line so that the line profile can be calculated
    line->nu_ul = nu_cl;
    line->j_ul = j_ul;
    line->kappa_lu = kappa_lu;
    line->calculate_spectrum( X );
    
    return;
}

double
SigmaTripletLBLDiatomicBand::
get_HLF( int tJu, int tJl )
{
    // Calculate Honl-London factor for this transition
    // See JQRST v9 pp 775-798 Arnold et al 1969 for HLF expressions used here
    
    int delta_J = (tJl - tJu)/2;
    double S = 0.0;
    double Ju = double(tJu)/2.0;
    
    // CHECKME: not sure that this is correct
    // P, Q and R branches respectively, ignoring spin splitting (for a parallel transition)
    if      ( delta_J == -1 ) S = Ju + 1.0;
    else if ( delta_J ==  0 ) S = 0.0;
    else if ( delta_J ==  1 ) S = Ju;
    
    return S;
}

// -----------------

HundBDoubletLBLDiatomicBand::
HundBDoubletLBLDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		                 double Re_vib, double m_w, double I, double sigma_nm  )
: LBLDiatomicBand( Vu, Vl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm )
{
    // NOTE: we are assuming this is always a Sigma-Sigma transition
}

void
HundBDoubletLBLDiatomicBand::
write_lines_to_file( string fname )
{
    /* Setup the output file */
    ofstream ofile;
    ofile << setprecision(6) << scientific << showpoint;
    ofile.open(fname.c_str());
    ofile << "# " << fname << endl
          << "# Ju_max = " << 0.5*double(tJu_max) << endl
          << "# Jl_max = " << 0.5*double(tJl_max) << endl
	  << "# Column 1: K_u" << endl
	  << "# Column 2: J_u" << endl
	  << "# Column 3: K_l" << endl
	  << "# Column 4: J_l" << endl
	  << "# Column 5: S" << endl
	  << "# Column 6: lambda_cl" << endl
	  << "# Column 7: A_ul" << endl;
    
    // Quantum numbers:
    // K -> total angular momentum apart from spin
    // S -> electronic spin angular momentum
    // J -> total angular momentum
    
    double S_sum = 0.0;
    int tJu_sum = 11;
    
    // 1.  Upper state 
    // 1a. Loop over upper state 'K' quantum number
    for ( int tKu = 2*lambda_u; tKu <= tJu_max + tS; tKu += 2 ) {
    	// 1b. Loop over upper state 'J' quantum number
    	for ( int tJu = abs( tKu - tS ); tJu <= ( tKu + tS ); tJu += 2 ) {
    	    // 1c. Check that this upper state is permitted
    	    if ( tJu > tJu_max ) continue;
    	    // 2. Apply delta_J transition rule 
    	    for ( int delta_J = -1; delta_J <= 1; delta_J+=1 ) {
    	    	// 3. Lower state
    	    	// 3a. calculate lower state 'J' quantum number
    	    	int tJl = tJu + 2*delta_J;
    	    	// 3b. Check that it is in-range
    	    	if ( tJl < 0 || tJl > tJl_max ) continue;
    	    	// 3c. Impose no 0-0 transition restriction
    	    	if ( tJu==0 && tJl==0 ) continue;
    	    	// 3d. Loop over possible lower state 'K' values
    	    	for ( int tKl = (tJl - tS); tKl <= (tJl + tS); tKl += 2 ) {
    	    	    // 3e. Check that Kl is above the required minimum
    	    	    if ( tKl < 2*lambda_l ) continue;
    	    	    // 3f. Impose |delta_K| = 1 transition restriction
    	    	    if ( abs(tKu-tKl)!=2 ) continue;
    	    	    /* If we have got to here this is an allowed transition */
    	    	    // 4. Calculate tSig_u and tSig_l
    	    	    int tSig_u = tJu - tKu;
    	    	    int tSig_l = tJl - tKl;
		    double S = get_HLF( tJu, tJl, tSig_u, tSig_l );
		    if ( tJu==tJu_sum ) S_sum += S;
		    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundBDoublet(Vu,tJu,tSig_u) - elev_l->calculate_E_rot_HundBDoublet(Vl,tJl,tSig_l) ) / RC_h_SI;
		    double lambda_cl = RC_c/nu_cl;
		    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );
		    ofile << setw(20)  << 0.5*tKu << setw(20) << 0.5*tJu 
		          << setw(20)  << 0.5*tKl << setw(20) << 0.5*tJl
		          << setw(20) << S
		          << setw(20) << lambda_cl*1.0e7 
		          << setw(20) << A_ul
		          << endl;
		}
	    }
	}
    }
    
    cout << "S_sum = " << S_sum << ", double(tJu+1) = " << double(tJu_sum+1) << endl;
    
    ofile.close();
    
    return;
}

void
HundBDoubletLBLDiatomicBand::
calculate_spectrum( CoeffSpectra &X, double Tv, double Tr )
{
    // Quantum numbers:
    // K -> total angular momentum apart from spin
    // S -> electronic spin angular momentum
    // J -> total angular momentum
    
    double j_ul_00 = -1.0;
    // 1.  Upper state 
    // 1a. Loop over upper state 'K' quantum number
    for ( int tKu = 2*lambda_u; tKu <= tJu_max + tS; tKu += 2 ) {
    	// 1b. Loop over upper state 'J' quantum number
    	for ( int tJu = abs( tKu - tS ); tJu <= ( tKu + tS ); tJu += 2 ) {
    	    // 1c. Check that this upper state is permitted
    	    if ( tJu > tJu_max ) continue;
    	    // 2. Apply delta_J transition rule 
    	    for ( int delta_J = -1; delta_J <= 1; delta_J+=1 ) {
    	    	// 3. Lower state
    	    	// 3a. calculate lower state 'J' quantum number
    	    	int tJl = tJu + 2*delta_J;
    	    	// 3b. Check that it is in-range
    	    	if ( tJl < 0 || tJl > tJl_max ) continue;
    	    	// 3c. Impose no 0-0 transition restriction
    	    	if ( tJu==0 && tJl==0 ) continue;
    	    	// 3d. Loop over possible lower state 'K' values
    	    	for ( int tKl = (tJl - tS); tKl <= (tJl + tS); tKl += 2 ) {
    	    	    // 3e. Check that Kl is above the required minimum
    	    	    if ( tKl < 2*lambda_l ) continue;
    	    	    // 3f. Impose |delta_K| = 1 transition restriction
    	    	    if ( abs(tKu-tKl)!=2 ) continue;
    	    	    /* If we have got to here this is an allowed transition */
    	    	    // 4. Calculate tSig_u and tSig_l
    	    	    int tSig_u = tJu - tKu;
    	    	    int tSig_l = tJl - tKl;
    	    	    // 5. Calculate the spectrum
    	    	    calculate_rot_line_spectrum(X,Tv,Tr,tJu,tJl,tSig_u,tSig_l,j_ul_00);
    	    	}
    	    }
    	}
    }
    
    return;
}

void
HundBDoubletLBLDiatomicBand::
calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, int tSig_u, int tSig_l, double &j_ul_00 )
{
    // 1. Calculate integrated emission and absorption coefficients
    // 1a. calculate Einstein coefficients
    double S = get_HLF( tJu, tJl, tSig_u, tSig_l );
    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundBDoublet(Vu,tJu,tSig_u) - elev_l->calculate_E_rot_HundBDoublet(Vl,tJl,tSig_l) ) / RC_h_SI;
    // NOTE: lambda must have cgs units
    double lambda_cl = RC_c/nu_cl;
    // NOTE: Spin/lambda splitting degeneracy should be included in Honl-London factor S
    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );

    // only degeneracy ratio is used, so don't worry about correcting for spin splitting
    int g_u = elev_u->get_g() * ( tJu + 1 );
    int g_l = elev_l->get_g() * ( tJl + 1 );
    
    double B_lu = A_ul * double( g_u ) / double( g_l ) * RC_c_SI * RC_c_SI 
		/ ( 8.0 * M_PI * RC_h_SI ) / ( nu_cl*nu_cl*nu_cl );
		
    double B_ul = B_lu * double( g_l )  / double ( g_u );
    
    // 1b. Upper and lower rotational state number densities
    double n_u = elev_u->calculate_N_rot_HundBDoublet(Tv,Tr,Vu,tJu,tSig_u,true);
    double n_l = elev_l->calculate_N_rot_HundBDoublet(Tv,Tr,Vl,tJl,tSig_l);
    
    // 1c. Combine the pieces to get integrated coefficients
    double j_ul = n_u * A_ul * RC_h_SI * nu_cl / ( 4.0 * M_PI );
    double kappa_lu = ( n_l * B_lu - n_u * B_ul ) * RC_h_SI * nu_cl;
    
    if ( j_ul_00 < 0.0 ) j_ul_00 = j_ul;
    else if ( j_ul < F_ROT_LINE_LIMIT * j_ul_00 ) return;
    
    // 2. Give coefficients to the line so that the line profile can be calculated
    line->nu_ul = nu_cl;
    line->j_ul = j_ul;
    line->kappa_lu = kappa_lu;
    line->calculate_spectrum( X );
    
    return;
}

double
HundBDoubletLBLDiatomicBand::
get_HLF( int tJu, int tJl, int tSig_u, int tSig_l )
{
    // Calculate Honl-London factor for this transition
    // Ref: Huber and Herzberg p 250
    
    int delta_J = (tJl - tJu)/2;
    
    double S = 0.0;
    double Jl = double(tJl)/2.0;
    
    // P, Q and R branches respectively
    if      ( delta_J ==  1 ) S = ( Jl * Jl - 0.25 ) / Jl;
    else if ( delta_J ==  0 ) S = ( 2.0*Jl + 1.0 ) / ( 4.0 * Jl * ( Jl + 1.0 ) );
    else if ( delta_J == -1 ) S = ( ( Jl + 1.0 ) * ( Jl + 1.0 ) - 0.25 ) / ( Jl + 1.0 );
    
    // 2. Divide by 2 to ensure summation = 2J+1
    return S / 2.0;
}

// --------------

HundABDoubletLBLDiatomicBand::
HundABDoubletLBLDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		 double Re_vib, double m_w, double I, double sigma_nm  )
: LBLDiatomicBand( Vu, Vl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm )
{}

void
HundABDoubletLBLDiatomicBand::
write_lines_to_file( string fname )
{
    /* Setup the output file */
    ofstream ofile;
    ofile << setprecision(6) << scientific << showpoint;
    ofile.open(fname.c_str());
    ofile << "# " << fname << endl
          << "# Ju_max = " << 0.5*double(tJu_max) << endl
          << "# Jl_max = " << 0.5*double(tJl_max) << endl
	  << "# Column 1: K_u" << endl
	  << "# Column 2: J_u" << endl
	  << "# Column 3: K_l" << endl
	  << "# Column 4: J_l" << endl
	  << "# Column 5: S" << endl
	  << "# Column 6: S_sum" << endl
	  << "# Column 7: 2 J_u + 1" << endl
	  << "# Column 7: lambda_cl" << endl
	  << "# Column 8: A_ul" << endl;
    
    // Quantum numbers:
    // K -> total angular momentum apart from spin
    // S -> electronic spin angular momentum
    // J -> total angular momentum
    
    double S_sum = 0.0;
    int tJu_sum = 11;
    
    // 1.  Upper state 
    // 1a. Loop over upper state 'K' quantum number
    for ( int tKu = 2*lambda_u; tKu <= tJu_max; tKu += 2 ) {
    	// 1b. Loop over upper state 'J' quantum number
    	for ( int tJu = abs( tKu - tS ); tJu <= ( tKu + tS ); tJu += 2 ) {
    	    // 1c. Check that this upper state is permitted
    	    double S_sum_Ju = 0.0;
    	    if ( tJu > tJu_max ) continue;
    	    // 2. Apply delta_J transition rule
    	    for ( int delta_J = -1; delta_J <= 1; ++delta_J ) {
    	    	// 3. Lower state
    	    	// 3a. calculate lower state 'J' quantum number
    	    	int tJl = tJu + 2*delta_J;
    	    	// 3b. Check that it is in-range
    	    	if ( tJl < 0 || tJl > tJl_max ) continue;
    	    	// 3c. Impose no 0-0 transition restriction
    	    	if ( tJu==0 && tJl==0 ) continue;
    	    	// 3d. Loop over possible lower state 'K' values
    	    	for ( int tKl = (tJl - tS); tKl <= (tJl + tS); tKl += 2 ) {
    	    	    // 3e. Check that Kl is above the required minimum
    	    	    if ( tKl < 2*lambda_l ) continue;
    	    	    /* If we have got to here this is an allowed transition */
    	    	    // 4. Calculate tSig_u and tSig_l
    	    	    int tSig_u = tJu - tKu;
    	    	    int tSig_l = tJl - tKl;
		    double S = get_HLF( tJu, tJl, tSig_u, tSig_l );
		    if ( tJu==tJu_sum ) S_sum += S;
		    S_sum_Ju += S;
		    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundABDoublet(Vu,tJu,tSig_u) - elev_l->calculate_E_rot_HundABDoublet(Vl,tJl,tSig_l) ) / RC_h_SI;
		    double lambda_cl = RC_c/nu_cl;
		    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );
		    ofile << setw(20) << 0.5*tKu << setw(20) << 0.5*tJu 
		          << setw(20) << 0.5*tKl << setw(20) << 0.5*tJl
		          << setw(20) << S
		          << setw(20) << S_sum_Ju
		          << setw(20) << tJu + 1
		          << setw(20) << lambda_cl*1.0e7 
		          << setw(20) << A_ul
		          << endl;
		}
	    }
	}
    }
    
    cout << "S_sum = " << S_sum << ", double(tJu+1) = " << double(tJu_sum+1) << endl;
    
    ofile.close();
    
    return;
}

void
HundABDoubletLBLDiatomicBand::
calculate_spectrum( CoeffSpectra &X, double Tv, double Tr )
{
    // Quantum numbers:
    // K -> total angular momentum apart from spin
    // S -> electronic spin angular momentum
    // J -> total angular momentum
    
    double j_ul_00 = -1.0;
    // 1.  Upper state 
    // 1a. Loop over upper state 'K' quantum number
    for ( int tKu = 2*lambda_u; tKu <= tJu_max + tS; tKu += 2 ) {
    	// 1b. Loop over upper state 'J' quantum number
    	for ( int tJu = abs( tKu - tS ); tJu <= ( tKu + tS ); tJu += 2 ) {
    	    // 1c. Check that this upper state is permitted
    	    if ( tJu > tJu_max ) continue;
    	    // 2. Apply delta_J transition rule
    	    for ( int delta_J = -1; delta_J <= 1; ++delta_J ) {
    	    	// 3. Lower state
    	    	// 3a. calculate lower state 'J' quantum number
    	    	int tJl = tJu + 2*delta_J;
    	    	// 3b. Check that it is in-range
    	    	if ( tJl < 0 || tJl > tJl_max ) continue;
    	    	// 3c. Impose no 0-0 transition restriction
    	    	if ( tJu==0 && tJl==0 ) continue;
    	    	// 3d. Loop over possible lower state 'K' values
    	    	for ( int tKl = (tJl - tS); tKl <= (tJl + tS); tKl += 2 ) {
    	    	    // 3e. Check that Kl is above the required minimum
    	    	    if ( tKl < 2*lambda_l ) continue;
    	    	    /* If we have got to here this is an allowed transition */
    	    	    // 4. Calculate tSig_u and tSig_l
    	    	    int tSig_u = tJu - tKu;
    	    	    int tSig_l = tJl - tKl;
    	    	    // 5. Calculate the spectrum
    	    	    calculate_rot_line_spectrum(X,Tv,Tr,tJu,tJl,tSig_u,tSig_l,j_ul_00);
    	    	}
    	    }
    	}
    }
    
    return;
}

void
HundABDoubletLBLDiatomicBand::
calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, int tSig_u, int tSig_l, double &j_ul_00 )
{
    // 1. Calculate integrated emission and absorption coefficients
    // 1a. calculate Einstein coefficients
    double S = get_HLF( tJu, tJl, tSig_u, tSig_l );
    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundABDoublet(Vu,tJu,tSig_u) - elev_l->calculate_E_rot_HundABDoublet(Vl,tJl,tSig_l) ) / RC_h_SI;
    // NOTE: lambda must have cgs units
    double lambda_cl = RC_c/nu_cl;
    // NOTE: Spin/lambda splitting degeneracy should be included in Honl-London factor S
    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );

    // only degeneracy ratio is used, so don't worry about correcting for spin splitting
    int g_u = elev_u->get_g() * ( tJu + 1 );
    int g_l = elev_l->get_g() * ( tJl + 1 );
    
    double B_lu = A_ul * double( g_u ) / double( g_l ) * RC_c_SI * RC_c_SI 
		/ ( 8.0 * M_PI * RC_h_SI ) / ( nu_cl*nu_cl*nu_cl );
		
    double B_ul = B_lu * double( g_l )  / double ( g_u );
    
    // 1b. Upper and lower rotational state number densities
    double n_u = elev_u->calculate_N_rot_HundABDoublet(Tv,Tr,Vu,tJu,tSig_u,true);
    double n_l = elev_l->calculate_N_rot_HundABDoublet(Tv,Tr,Vl,tJl,tSig_l);
    
    // 1c. Combine the pieces to get integrated coefficients
    double j_ul = n_u * A_ul * RC_h_SI * nu_cl / ( 4.0 * M_PI );
    double kappa_lu = ( n_l * B_lu - n_u * B_ul ) * RC_h_SI * nu_cl;
    
    if ( j_ul_00 < 0.0 ) j_ul_00 = j_ul;
    else if ( j_ul < F_ROT_LINE_LIMIT * j_ul_00 ) return;
    
    // 2. Give coefficients to the line so that the line profile can be calculated
    line->nu_ul = nu_cl;
    line->j_ul = j_ul;
    line->kappa_lu = kappa_lu;
    line->calculate_spectrum( X );
    
    return;
}

double
HundABDoubletLBLDiatomicBand::
get_HLF( int tJu, int tJl, int tSig_u, int tSig_l)
{
#   if HUND_AB_DOUBLET_HLF_METHOD==0
    // Calculate Honl-London factor for this transition
    // See JQRST v9 pp 775-798 Arnold et al 1969 for HLF expressions used here
    // NOTE: these expression are for 2Sigma - 2Pi transitions only
    
    int delta_J = (tJl - tJu)/2;
    double S = 0.0;
    double Ju = double(tJu)/2.0;
    
    double sign = 1.0;
    int expr = -1;
    
    if      ( delta_J==1 && ( tSig_u==-1 && tSig_l==-1 ) ) {
	// P2 branch
	sign = 1.0;
	if      ( lambda_u==1 ) expr = 0;
	else if ( lambda_u==0 ) expr = 4;
    }
    else if ( delta_J==1 && ( tSig_u==+1 && tSig_l==-1 ) ) {
	// OP12 branch
	sign = -1.0;
	if      ( lambda_u==1 ) expr = 0;
	else if ( lambda_u==0 ) expr = 5;
    }
    else if ( delta_J==1 && ( tSig_u==-1 && tSig_l==+1 ) ) {
	// QP21 branch
	sign = -1.0;
	if      ( lambda_u==1 ) expr = 1;
	else if ( lambda_u==0 ) expr = 4;
    }
    else if ( delta_J==1 && ( tSig_u==+1 && tSig_l==+1 ) ) {
	// P1 branch
	sign = 1.0;
	if      ( lambda_u==1 ) expr = 1;
	else if ( lambda_u==0 ) expr = 5;
    }
    else if ( delta_J==0 && ( tSig_u==-1 && tSig_l==-1 ) ) {
	// Q2 branch
	sign = 1.0;
	if      ( lambda_u==1 ) expr = 2;
	else if ( lambda_u==0 ) expr = 2;
    }
    else if ( delta_J==0 && ( tSig_u==+1 && tSig_l==-1 ) ) {
	// PQ12 branch
	sign = -1.0;
	if      ( lambda_u==1 ) expr = 2;
	else if ( lambda_u==0 ) expr = 3;
    }
    else if ( delta_J==0 && ( tSig_u==-1 && tSig_l==+1 ) ) {
	// RQ21 branch
	sign = -1.0;
	if      ( lambda_u==1 ) expr = 3;
	else if ( lambda_u==0 ) expr = 2;
    }
    else if ( delta_J==0 && ( tSig_u==+1 && tSig_l==+1 ) ) {
	// Q1 branch
	sign = 1.0;
	if      ( lambda_u==1 ) expr = 3;
	else if ( lambda_u==0 ) expr = 3;
    }
    else if ( delta_J==-1 && ( tSig_u==-1 && tSig_l==-1 ) ) {
	// R2 branch
	sign = 1.0;
	if      ( lambda_u==1 ) expr = 4;
	else if ( lambda_u==0 ) expr = 0;
    }
    else if ( delta_J==-1 && ( tSig_u==+1 && tSig_l==-1 ) ) {
	// QR12 branch
	sign = -1.0;
	if      ( lambda_u==1 ) expr = 4;
	else if ( lambda_u==0 ) expr = 1;
    }
    else if ( delta_J==-1 && ( tSig_u==-1 && tSig_l==+1 ) ) {
	// SR21 branch
	sign = -1.0;
	if      ( lambda_u==1 ) expr = 5;
	else if ( lambda_u==0 ) expr = 0;
    }
    else if ( delta_J==-1 && ( tSig_u==+1 && tSig_l==+1 ) ) {
	// R1 branch
	sign = 1.0;
	if      ( lambda_u==1 ) expr = 5;
	else if ( lambda_u==0 ) expr = 1;
    }
    
    /* Now use the expr and sign values we have just assigned to calculate S */
    
    double Y = elev_u->get_A_spin() / elev_u->calculate_B_v(Vu);
    double U = 1.0 / sqrt( Y*Y - 4.0*Y + (2.0*Ju+1.0)*(2.0*Ju+1.0) );
    
    // Precompute the universally used terms with int->double conversion
    double tmpA = 2.0*Ju + 1.0;
    double tmpB = 4.0*Ju*Ju + 4.0*Ju;
    
    if      ( expr==0 )
	S = ( tmpA*tmpA + sign*tmpA*U* ( tmpB + 1.0 - 2.0*Y ) ) / ( 16.0*( Ju + 1.0 ) );
    else if ( expr==1 )
	S = ( tmpA*tmpA + sign*tmpA*U* ( tmpB - 7.0 + 2.0*Y ) ) / ( 16.0*( Ju + 1.0 ) );
    else if ( expr==2 )
	S = tmpA*( ( tmpB - 1.0 ) + sign*U*( 8.0*Ju*Ju*Ju + 12.0*Ju*Ju - 2.0*Ju + 1.0 - 2.0*Y ) ) / ( 16.0*Ju*(Ju + 1.0) );
    else if ( expr==3 )
	S = tmpA*( ( tmpB - 1.0 ) + sign*U*( 8.0*Ju*Ju*Ju + 12.0*Ju*Ju - 2.0*Ju - 7.0 + 2.0*Y ) ) / ( 16.0*Ju*(Ju + 1.0) );
    else if ( expr==4 )
	S = ( tmpA*tmpA + sign*tmpA*U*( tmpB - 7.0 + 2.0*Y ) ) / ( 16.0*Ju );
    else if ( expr==5 )
	S = ( tmpA*tmpA + sign*tmpA*U*( tmpB + 1.0 - 2.0*Y ) ) / ( 16.0*Ju );
    else {
	cout << "HundABDoubletLBLDiatomicBand::get_HLF()" << endl 
	     << "Failed selecting an expression for 2Pi<->2Sigma transition" << endl
	     << "Check electronic level input data" << endl
	     << "Exiting program." << endl;
	exit( BAD_INPUT_ERROR );
    }
    // Normalise by 2S+1 (always 2) here so it doesn't have to be done in the A_ul calculation
    S /= 2.0;
    
#   else

    // Kovacs (1969)
    int delta_L = lambda_l - lambda_u;
    int delta_J = (tJl - tJu)/2;
    double S = 0.0;
    double J_u = double(tJu)/2.0;
    double J_l = double(tJl)/2.0;
    double L_l = double(lambda_l);
    double A_u = elev_u->get_A_spin();
    double A_l = elev_l->get_A_spin();
    double B_v_u = elev_u->calculate_B_v(Vu);
    double B_v_l = elev_l->calculate_B_v(Vl);
    
    if ( delta_L==0 ) {
    	// Parallel transition
    	if ( delta_J==1 ) {
    	    // P branches
    	    double tmpA = (J_l-L_l-0.5)*(J_l+L_l+0.5);
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// P1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u-1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l-0.5);
    	    	tmpD = 4.0 * J_l * C_minus(lambda_u, A_u, B_v_u, J_u-1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// P21 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u-1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l-0.5);
    	    	tmpD = 4.0 * J_l * C_plus(lambda_u, A_u, B_v_u, J_u-1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// P12 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u-1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l-0.5);
    	    	tmpD = 4.0 * J_l * C_minus(lambda_u, A_u, B_v_u, J_u-1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// P2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u-1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l-0.5);
    	    	tmpD = 4.0 * J_l * C_plus(lambda_u, A_u, B_v_u, J_u-1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    	else if ( delta_J==0 ) {
    	    // Q branches
    	    double tmpA = J_l + 0.5;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// Q1 branch
    	    	tmpB = ( L_l + 0.5 ) * u_minus(lambda_u, A_u, B_v_u, J_u) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0 * ( L_l - 0.5 ) * ( J_l - L_l + 0.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 2.0 * J_l * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u) * C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// Q21 branch
    	    	tmpB = ( L_l + 0.5 ) * u_plus(lambda_u, A_u, B_v_u, J_u) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = - 4.0 * ( L_l - 0.5 ) * ( J_l - L_l + 0.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 2.0 * J_l * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u) * C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// Q12 branch
    	    	tmpB = ( L_l + 0.5 ) * u_minus(lambda_u, A_u, B_v_u, J_u) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0 * ( L_l - 0.5 ) * ( J_l - L_l + 0.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 2.0 * J_l * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u) * C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// Q2 branch
    	    	tmpB = ( L_l + 0.5 ) * u_plus(lambda_u, A_u, B_v_u, J_u) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0 * ( L_l - 0.5 ) * ( J_l - L_l + 0.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 2.0 * J_l * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u) * C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    	else if ( delta_J==-1 ) {
    	    // R branches
    	    double tmpA = ( J_l - L_l + 0.5 ) * ( J_l + L_l + 1.5 );
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// R1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u+1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0 * ( J_l - L_l + 1.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 4.0 * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u+1.0) * C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// R21 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u+1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0 * ( J_l - L_l + 1.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 4.0 * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u+1.0) * C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// R12 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u+1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0 * ( J_l - L_l + 1.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 4.0 * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u+1.0) * C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// R2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u+1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0 * ( J_l - L_l + 1.5 ) * ( J_l + L_l + 0.5 );
    	    	tmpD = 4.0 * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u+1.0) * C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    }
    else if ( delta_L==1 ) {
    	// Perpendicular plus one transition
    	if ( delta_J==1 ) {
    	    // P branches
    	    double tmpA = (J_l-L_l-1.5)*(J_l-L_l-0.5);
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// P1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u-1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_minus(lambda_u, A_u, B_v_u, J_u-1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// P21 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u-1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_plus(lambda_u, A_u, B_v_u, J_u-1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// P12 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u-1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_minus(lambda_u, A_u, B_v_u, J_u-1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// P2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u-1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_plus(lambda_u, A_u, B_v_u, J_u-1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    	else if ( delta_J==0 ) {
    	    // Q branches
    	    double tmpA = (J_l-L_l-0.5)*(J_l+0.5)*(J_l+L_l+1.5);
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// Q1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// Q21 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// Q12 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// Q2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    	else if ( delta_J==-1 ) {
    	    // R branches
    	    double tmpA = (J_l+L_l+1.5)*(J_l+L_l+2.5);
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// R1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u+1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u+1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// R21 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u+1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u+1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// R12 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u+1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u+1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// R2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u+1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u+1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    }
    else if ( delta_L==-1 ) {
    	// Perpendicular minus one transition
    	if ( delta_J==-1 ) {
    	    // R branches -> increment J indices by one so we can use delta_l=1 P branch equations
    	    J_l += 1.0; J_u += 1.0;
    	    double tmpA = (J_l-L_l-1.5)*(J_l-L_l-0.5);
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// R1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u-1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_minus(lambda_u, A_u, B_v_u, J_u-1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// R12 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u-1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_plus(lambda_u, A_u, B_v_u, J_u-1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// R21 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u-1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_minus(lambda_u, A_u, B_v_u, J_u-1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// R2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u-1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * J_l * C_plus(lambda_u, A_u, B_v_u, J_u-1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    	else if ( delta_J==0 ) {
    	    // Q branches
    	    double tmpA = (J_l-L_l-0.5)*(J_l+0.5)*(J_l+L_l+1.5);
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// Q1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// Q12 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// Q21 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// Q2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 4.0 * J_l * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    	else if ( delta_J==1 ) {
    	    // P branches -> decrement J indices by one so we can use delta_l=1 R branch equations
    	    J_l -= 1.0; J_u -= 1.0;
    	    double tmpA = (J_l+L_l+1.5)*(J_l+L_l+2.5);
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    if ( tSig_u==1 && tSig_l==1 ) {
    	    	// P1 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u+1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u+1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==1 && tSig_l==-1 ) {
    	    	// P12 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u+1.0) * u_minus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u+1.0)*C_minus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==1 ) {
    	    	// P21 branch
    	    	tmpB = u_minus(lambda_u, A_u, B_v_u, J_u+1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = -4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_minus(lambda_u, A_u, B_v_u, J_u+1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    else if ( tSig_u==-1 && tSig_l==-1 ) {
    	    	// P2 branch
    	    	tmpB = u_plus(lambda_u, A_u, B_v_u, J_u+1.0) * u_plus(lambda_l, A_l, B_v_l, J_l);
    	    	tmpC = 4.0*(J_l-L_l+0.5)*(J_l+L_l+0.5);
    	    	tmpD = 8.0 * ( J_l + 1.0 ) * C_plus(lambda_u, A_u, B_v_u, J_u+1.0)*C_plus(lambda_l, A_l, B_v_l, J_l);
    	    }
    	    S = tmpA * ( tmpB + tmpC ) * ( tmpB + tmpC ) / tmpD;
    	}
    }
    S /= 2.0;
#   endif
    
    return S;
}

double
HundABDoubletLBLDiatomicBand::
u_plus( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) Y *= -1.0;

    return sqrt( lambda*lambda*Y*( Y - 4.0 ) + 4.0*( J + 0.5 )*( J + 0.5 ) ) + lambda * ( Y - 2.0 );
}

double
HundABDoubletLBLDiatomicBand::
u_minus( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) Y *= -1.0;

    return sqrt( lambda*lambda*Y*( Y - 4.0 ) + 4.0*( J + 0.5 )*( J + 0.5 ) ) - lambda * ( Y - 2.0 );
}

double
HundABDoubletLBLDiatomicBand::
C_plus( double lambda, double A, double B_v, double J )
{
    double u = u_plus( lambda, A, B_v, J );

    return 0.5 * ( u*u + 4.0 * ( ( J + 0.5 ) * ( J + 0.5 ) - lambda*lambda));
}

double
HundABDoubletLBLDiatomicBand::
C_minus( double lambda, double A, double B_v, double J )
{
    double u = u_minus( lambda, A, B_v, J );

    return 0.5 * ( u*u + 4.0 * ( ( J + 0.5 ) * ( J + 0.5 ) - lambda*lambda));
}

// --------------

HundABTripletLBLDiatomicBand::
HundABTripletLBLDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		 double Re_vib, double m_w, double I, double sigma_nm  )
: LBLDiatomicBand( Vu, Vl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm )
{}

void
HundABTripletLBLDiatomicBand::
write_lines_to_file( string fname )
{
    /* Setup the output file */
    ofstream ofile;
    ofile << setprecision(6) << scientific << showpoint;
    ofile.open(fname.c_str());
    ofile << "# " << fname << endl
          << "# Ju_max = " << 0.5*double(tJu_max) << endl
          << "# Jl_max = " << 0.5*double(tJl_max) << endl
	  << "# Column 1: K_u" << endl
	  << "# Column 2: J_u" << endl
	  << "# Column 3: K_l" << endl
	  << "# Column 4: J_l" << endl
	  << "# Column 5: S" << endl
	  << "# Column 6: lambda_cl" << endl
	  << "# Column 7: A_ul" << endl;
    
    // Quantum numbers:
    // K -> total angular momentum apart from spin
    // S -> electronic spin angular momentum
    // J -> total angular momentum
    
    double S_sum = 0.0;
    int tJu_sum = 50;
    
    // 1.  Upper state 
    // 1a. Loop over upper state 'K' quantum number
    for ( int tKu = 2*lambda_u; tKu <= tJu_max; tKu += 2 ) {
    	// 1b. Loop over upper state 'J' quantum number
    	for ( int tJu = abs( tKu - tS ); tJu <= ( tKu + tS ); tJu += 2 ) {
    	    // 1c. Check that this upper state is permitted
    	    if ( tJu > tJu_max ) continue;
    	    // 2. Apply delta_J transition rule
    	    for ( int delta_J = -1; delta_J <= 1; ++delta_J ) {
    	    	// 3. Lower state
    	    	// 3a. calculate lower state 'J' quantum number
    	    	int tJl = tJu + 2*delta_J;
    	    	// 3b. Check that it is in-range
    	    	if ( tJl < 0 || tJl > tJl_max ) continue;
    	    	// 3c. Impose no 0-0 transition restriction
    	    	if ( tJu==0 && tJl==0 ) continue;
    	    	// 3d. Loop over possible lower state 'K' values
    	    	for ( int tKl = (tJl - tS); tKl <= (tJl + tS); tKl += 2 ) {
    	    	    // 3e. Check that Kl is above the required minimum
    	    	    if ( tKl < 2*lambda_l ) continue;
    	    	    /* If we have got to here this is an allowed transition */
    	    	    // 4. Calculate tSig_u and tSig_l
    	    	    int tSig_u = tJu - tKu;
    	    	    int tSig_l = tJl - tKl;
		    double S = get_HLF( tJu, tJl, tSig_u, tSig_l );
		    if ( tJu == tJu_sum ) S_sum += S;
		    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundABTriplet(Vu,tJu,tSig_u) - elev_l->calculate_E_rot_HundABTriplet(Vl,tJl,tSig_l) ) / RC_h_SI;
		    double lambda_cl = RC_c/nu_cl;
		    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );
		    ofile << setw(20)  << 0.5*tKu << setw(20) << 0.5*tJu 
		          << setw(20)  << 0.5*tKl << setw(20) << 0.5*tJl
		          << setw(20) << S
		          << setw(20) << lambda_cl*1.0e7 
		          << setw(20) << A_ul
		          << endl;
		}
	    }
	}
    }
    
    cout << "S_sum = " << S_sum << ", double(tJu+1) = " << double(tJu_sum+1) << endl;
    
    ofile.close();
    
    return;
}

void
HundABTripletLBLDiatomicBand::
calculate_spectrum( CoeffSpectra &X, double Tv, double Tr )
{
    // Quantum numbers:
    // K -> total angular momentum apart from spin
    // S -> electronic spin angular momentum
    // J -> total angular momentum
    
    double j_ul_00 = -1.0;
    // 1.  Upper state 
    // 1a. Loop over upper state 'K' quantum number
    for ( int tKu = 2*lambda_u; tKu <= tJu_max + tS; tKu += 2 ) {
    	// 1b. Loop over upper state 'J' quantum number
    	for ( int tJu = abs( tKu - tS ); tJu <= ( tKu + tS ); tJu += 2 ) {
    	    // 1c. Check that this upper state is permitted
    	    if ( tJu > tJu_max ) continue;
    	    // 2. Apply delta_J transition rule
    	    for ( int delta_J = -1; delta_J <= 1; ++delta_J ) {
    	    	// 3. Lower state
    	    	// 3a. calculate lower state 'J' quantum number
    	    	int tJl = tJu + 2*delta_J;
    	    	// 3b. Check that it is in-range
    	    	if ( tJl < 0 || tJl > tJl_max ) continue;
    	    	// 3c. Impose no 0-0 transition restriction
    	    	if ( tJu==0 && tJl==0 ) continue;
    	    	// 3d. Loop over possible lower state 'K' values
    	    	for ( int tKl = (tJl - tS); tKl <= (tJl + tS); tKl += 2 ) {
    	    	    // 3e. Check that Kl is above the required minimum
    	    	    if ( tKl < 2*lambda_l ) continue;
    	    	    /* If we have got to here this is an allowed transition */
    	    	    // 4. Calculate tSig_u and tSig_l
    	    	    int tSig_u = tJu - tKu;
    	    	    int tSig_l = tJl - tKl;
    	    	    // 5. Calculate the spectrum
    	    	    calculate_rot_line_spectrum(X,Tv,Tr,tJu,tJl,tSig_u,tSig_l,j_ul_00);
    	    	}
    	    }
    	}
    }
    
    return;
}

void
HundABTripletLBLDiatomicBand::
calculate_rot_line_spectrum( CoeffSpectra &X, double Tv, double Tr, int tJu, int tJl, int tSig_u, int tSig_l, double &j_ul_00 )
{
    // 1. Calculate integrated emission and absorption coefficients
    // 1a. calculate Einstein coefficients
    double S = get_HLF( tJu, tJl, tSig_u, tSig_l );
    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_HundABTriplet(Vu,tJu,tSig_u) - elev_l->calculate_E_rot_HundABTriplet(Vl,tJl,tSig_l) ) / RC_h_SI;
    // NOTE: lambda must have cgs units
    double lambda_cl = RC_c/nu_cl;
    // NOTE: Spin/lambda splitting degeneracy should be included in Honl-London factor S
    double A_ul = A_ul_const * Re_vib * Re_vib * S / ( lambda_cl*lambda_cl*lambda_cl * double(tJu+1) );

    // only degeneracy ratio is used, so don't worry about correcting for spin splitting
    int g_u = elev_u->get_g() * ( tJu + 1 );
    int g_l = elev_l->get_g() * ( tJl + 1 );
    
    double B_lu = A_ul * double( g_u ) / double( g_l ) * RC_c_SI * RC_c_SI 
		/ ( 8.0 * M_PI * RC_h_SI ) / ( nu_cl*nu_cl*nu_cl );
		
    double B_ul = B_lu * double( g_l )  / double ( g_u );
    
    // 1b. Upper and lower rotational state number densities
    double n_u = elev_u->calculate_N_rot_HundABTriplet(Tv,Tr,Vu,tJu,tSig_u,true);
    double n_l = elev_l->calculate_N_rot_HundABTriplet(Tv,Tr,Vl,tJl,tSig_l);
    
    // 1c. Combine the pieces to get integrated coefficients
    double j_ul = n_u * A_ul * RC_h_SI * nu_cl / ( 4.0 * M_PI );
    double kappa_lu = ( n_l * B_lu - n_u * B_ul ) * RC_h_SI * nu_cl;
    
    if ( j_ul_00 < 0.0 ) j_ul_00 = j_ul;
    else if ( j_ul < F_ROT_LINE_LIMIT * j_ul_00 ) return;
    
    // 2. Give coefficients to the line so that the line profile can be calculated
    line->nu_ul = nu_cl;
    line->j_ul = j_ul;
    line->kappa_lu = kappa_lu;
    line->calculate_spectrum( X );
    
    return;
}

double
HundABTripletLBLDiatomicBand::
get_HLF( int tJu, int tJl, int tSig_u, int tSig_l)
{
    // Kovacs (1969)
    int delta_L = lambda_l - lambda_u;
    int delta_J = (tJl - tJu)/2;
    double S = 0.0;
    double J_u = double(tJu)/2.0;
    double J_l = double(tJl)/2.0;
    double L_u = double(lambda_u);
    double L_l = double(lambda_l);
    double A_u = elev_u->get_A_spin();
    double A_l = elev_l->get_A_spin();
    double B_v_u = elev_u->calculate_B_v(Vu);
    double B_v_l = elev_l->calculate_B_v(Vl);
    
    if ( delta_L==0 ) {
    	// Parallel transition
    	if ( delta_J==-1 ) {
    	    // P branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// P1 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 16.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l ) * (  J_l + L_l - 1.0 ) * ( J_l + L_l );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// P21 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 2.0 * J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// P31 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 16.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// P12 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 2.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l - 1.0 ) * ( J_l + L_l - 1.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// P2 branch [v]
    	    	tmpA = 4.0 *( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l * L_l * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpE = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// P32 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 2.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l ) * ( J_l + L_l ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// P13 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 16.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_1_plus( L_u, A_u, B_v_u, J_u ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = -8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// P23 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 2.0 * J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// P3 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l );
    	    	tmpB = 16.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l - 1.0 ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    	else if ( delta_J==0 ) {
    	    // Q branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// Q1 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 16.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * L_l * ( J_l - L_l ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// Q21 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 2.0 * J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * L_l * L_l * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// Q31 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 16.0 * J_l * ( J_l + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * L_l * ( J_l - 1.0 ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// Q12 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 2.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u );
    	    	tmpD = - ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u );
    	    	tmpE = - 2.0 * L_l * L_l * ( J_l - L_l ) * ( J_l + L_l ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// Q2 branch [v]
    	    	tmpA = 4.0 *( 2.0 * J_l + 1.0 );
    	    	tmpB = J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l * L_l * L_l * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 );
    	    	tmpE = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// Q32 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 2.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u );
    	    	tmpD = - ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u );
    	    	tmpE = 2.0 * L_l * L_l * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// Q13 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 16.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * L_l * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// Q23 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 2.0 * J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u );
    	    	tmpE = 2.0 * L_l * L_l * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// Q3 branch [v]
    	    	tmpA = 2.0 * J_l + 1.0;
    	    	tmpB = 16.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * L_l * ( J_l - L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    	else if ( delta_J==1 ) {
    	    // R branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// R1 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 16.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 2.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// R21 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 2.0 * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 2.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// R31 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 16.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 2.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l + 2.0 ) * ( J_l + L_l ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// R12 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 2.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 2.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// R2 branch [v]
    	    	tmpA = 4.0 * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l * L_l * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 2.0 );
    	    	tmpE = ( J_l - L_l + 2.0 ) * ( J_l + L_l );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// R32 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 2.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 2.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l + 2.0 ) * ( J_l + L_l + 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// R13 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 16.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 2.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// R23 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 2.0 * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 2.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// R3 branch [v]
    	    	tmpA = ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    	tmpB = 16.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 2.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l + 1.0 ) * ( J_l - L_l + 2.0 ) * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    }
    else if ( delta_L==1 ) {
    	// Perpendicular plus one transition
    	if ( delta_J==-1 ) {
    	    // P branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// P1 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l - 2.0 ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// P21 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// P31 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// P12 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l - 2.0 ) * ( J_l + L_l ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// P2 branch [v]
    	    	tmpA = 2.0 *( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l * ( L_l + 1.0 ) * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l + 1.0 ) * ( J_l + L_l );
    	    	tmpE = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// P32 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// P13 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = -8.0 * ( J_l - L_l - 2.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// P23 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * ( L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// P3 branch [v]
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    	else if ( delta_J==0 ) {
    	    // Q branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// Q1 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// Q21 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// Q31 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// Q12 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// Q2 branch [v]
    	    	tmpA = 2.0 * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l * ( L_l + 1.0 ) * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l + 1.0 ) * ( J_l + L_l );
    	    	tmpE = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// Q32 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// Q13 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// Q23 branch [v]
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * ( L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// Q3 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    	else if ( delta_J==1 ) {
    	    // R branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// R1 branch [v]
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l ) * ( J_l - L_l ) * ( J_l + L_l ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// R21 branch [v]
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// R31 branch [v]
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * ( J_l + L_l + 3.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// R12 branch [v]
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// R2 branch [v]
    	    	tmpA = 2.0 * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l *( L_l + 1.0 ) * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l + 1.0 ) * ( J_l + L_l );
    	    	tmpE = ( J_l - L_l ) * ( J_l + L_l + 3.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// R32 branch [v]
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 3.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// R13 branch [v]
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// R23 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * ( L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// R3 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 3.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    }
    else if ( delta_L==-1 ) {
    	// Perpendicular minus one transition
    	if ( delta_J==1 ) {
    	    // R branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    J_l += 1; J_u += 1;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// R1 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l - 2.0 ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// R12 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// R13 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// R21 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l - 2.0 ) * ( J_l + L_l ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// R2 branch
    	    	tmpA = 2.0 *( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l * ( L_l + 1.0 ) * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l + 1.0 ) * ( J_l + L_l );
    	    	tmpE = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// R23 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// R31 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_1( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = -8.0 * ( J_l - L_l - 2.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// R32 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 4.0 * J_l * C_2( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * ( L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// R3 branch
    	    	tmpA = ( J_l - L_l - 1.0 ) * ( J_l - L_l );
    	    	tmpB = 32.0 * J_l * C_3( L_u, A_u, B_v_u, J_u - 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 2.0 ) * ( J_l + L_l + 1.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u - 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    	else if ( delta_J==0 ) {
    	    // Q branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// Q1 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// Q12 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// Q13 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l ) * (  J_l + L_l ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// Q21 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l - 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// Q2 branch
    	    	tmpA = 2.0 * ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l * ( L_l + 1.0 ) * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l + 1.0 ) * ( J_l + L_l );
    	    	tmpE = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// Q23 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// Q31 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l - 1.0 ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l + 1.0 ) * ( J_l + L_l + 1.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// Q32 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 4.0 * J_l * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * ( L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// Q3 branch
    	    	tmpA = ( J_l - L_l ) * ( J_l + L_l + 1.0 ) * ( 2.0 * J_l + 1.0 );
    	    	tmpB = 32.0 * J_l * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l - 1.0 ) * ( J_l + L_l + 2.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * (  J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    	else if ( delta_J==1 ) {
    	    // P branches
    	    double tmpA = 0.0;
    	    double tmpB = 0.0;
    	    double tmpC = 0.0;
    	    double tmpD = 0.0;
    	    double tmpE = 0.0;
    	    J_l -= 1.0; J_u -= 1.0;
    	    if ( tSig_u==2 && tSig_l==2 ) {
    	    	// P1 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l ) * ( J_l - L_l ) * ( J_l + L_l ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==0 ) {
    	    	// P12 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 2.0 * ( L_l + 1.0 ) * ( J_l - L_l ) * ( J_l + L_l ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==2 && tSig_l==-2 ) {
    	    	// P13 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_1( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_1_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * ( J_l + L_l + 3.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==2 ) {
    	    	// P21 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpE = - 2.0 * L_l * ( J_l - L_l ) * ( J_l + L_l + 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==0 ) {
    	    	// P2 branch
    	    	tmpA = 2.0 * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = 0.5 * L_l *( L_l + 1.0 ) * ( A_u / B_v_u - 2.0 ) * ( A_l / B_v_l - 2.0 );
    	    	tmpD = ( J_l - L_l + 1.0 ) * ( J_l + L_l );
    	    	tmpE = ( J_l - L_l ) * ( J_l + L_l + 3.0 );
    	    }
    	    else if ( tSig_u==0 && tSig_l==-2 ) {
    	    	// P23 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_2( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 );
    	    	tmpE = 2.0 * L_l * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 3.0 ) * ( A_l / B_v_l - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==2 ) {
    	    	// P31 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_1( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( L_l - L_l + 1.0 ) * ( J_l + L_l ) * u_1_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_1_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = - 8.0 * ( J_l - L_l ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==0 ) {
    	    	// P32 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 4.0 * ( J_l + 1.0 ) * C_2( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = - ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 2.0 * ( L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( A_u / B_v_u - 2.0 );
    	    }
    	    else if ( tSig_u==-2 && tSig_l==-2 ) {
    	    	// P3 branch
    	    	tmpA = ( J_l + L_l + 1.0 ) * ( J_l + L_l + 2.0 );
    	    	tmpB = 32.0 * ( J_l + 1.0 ) * C_3( L_u, A_u, B_v_u, J_u + 1.0 ) * C_3( L_l, A_l, B_v_l, J_l );
    	    	tmpC = ( J_l - L_l + 1.0 ) * ( J_l + L_l ) * u_3_minus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_minus( L_l, A_l, B_v_l, J_l );
    	    	tmpD = ( J_l - L_l ) * ( J_l + L_l + 3.0 ) * u_3_plus( L_u, A_u, B_v_u, J_u + 1.0 ) * u_3_plus( L_l, A_l, B_v_l, J_l );
    	    	tmpE = 8.0 * ( J_l - L_l + 1.0 ) * ( J_l - L_l + 1.0 ) * ( J_l + L_l + 1.0 ) * ( J_l + L_l + 3.0 );
    	    }
    	    S = tmpA / tmpB * ( tmpC + tmpD + tmpE ) * ( tmpC + tmpD + tmpE );
    	}
    }
    
    return S;
}

double
HundABTripletLBLDiatomicBand::
u_1_plus( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) {
    	Y *= -1.0;
    }
    return sqrt( lambda*lambda*Y*( Y - 4.0 ) + 4.0*J*J ) + lambda * ( Y - 2.0 );
}

double
HundABTripletLBLDiatomicBand::
u_1_minus( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) {
    	Y *= -1.0;
    }
    return sqrt( lambda*lambda*Y*( Y - 4.0 ) + 4.0*J*J ) - lambda * ( Y - 2.0 );
}

double
HundABTripletLBLDiatomicBand::
u_3_plus( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) {
    	Y *= -1.0;
    }
    return sqrt( lambda*lambda*Y*( Y - 4.0 ) + 4.0*(J+1.0)*(J+1.0) ) + lambda * ( Y - 2.0 );
}

double
HundABTripletLBLDiatomicBand::
u_3_minus( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) {
    	Y *= -1.0;
    }
    return sqrt( lambda*lambda*Y*( Y - 4.0 ) + 4.0*(J+1.0)*(J+1.0) ) - lambda * ( Y - 2.0 );
}

double
HundABTripletLBLDiatomicBand::
C_1( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) {
    	Y *= -1.0;
    	lambda *= -1.0;
    }
    return lambda*lambda*Y*( Y - 4.0 ) * ( J - lambda + 1.0 ) * ( J + lambda ) + 2.0 * ( 2.0 * J + 1.0 ) * ( J - lambda ) * ( J + lambda );
}

double
HundABTripletLBLDiatomicBand::
C_2( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) {
    	Y *= -1.0;
    	// lambda *= -1.0;
    }
    return lambda*lambda*Y*( Y - 4.0 ) + 4.0 * J * ( J + 1.0 );
}

double
HundABTripletLBLDiatomicBand::
C_3( double lambda, double A, double B_v, double J )
{
    double Y = A / B_v;
    if ( A < 0.0 ) {
    	Y *= -1.0;
    	lambda *= -1.0;
    }
    return lambda*lambda*Y*( Y - 4.0 ) * ( J - lambda ) * ( J + lambda + 1.0 ) + 2.0 * ( 2.0 * J + 1.0 ) * ( J - lambda + 1.0 ) * ( J + 1.0 ) * ( J + lambda + 1.0 );
}

SRBDiatomicBand::
SRBDiatomicBand( int Vu, int Vl, DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, double Re_vib )
: DiatomicBand( Vu, Vl, elev_u, elev_l, Re_vib )
{}

SRBDiatomicBand::
~SRBDiatomicBand() {}


void
SRBDiatomicBand::
initialize(  double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    // Eventually calculate SB_constA, SB_constB and SB_constC
    
    return;
}

string
SRBDiatomicBand::
line_width_string(  double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    // return blank string
    
    return "\n";
}

void
SRBDiatomicBand::
write_lines_to_file( string fname )
{
    // Do nothing
    
    return;
}

void
SRBDiatomicBand::
calculate_spectrum( CoeffSpectra &X, double Tv, double Tr )
{
    
    return;
}

