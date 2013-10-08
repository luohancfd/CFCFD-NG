/** \file atomic_line.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 21-Aug-07
 *  \brief Definitions for Atomic_line class; simply a place-holder for transition details.
 *
 **/

#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#include "../../util/source/useful.h"

#include "atomic_line.hh"
#include "radiation_constants.hh"
#include "spectral_model.hh"
#include "atomic_radiator.hh"

using namespace std;

AtomicLine::
AtomicLine( vector<double> line_data, double m_w, double I, int npoints, int nwidths, double beta )
: m_w( m_w ), I( I ), npoints( npoints ), nwidths( nwidths )
{
    // 0. Check size of line data vector
    // -> now done by AtomicRadiator constructor
    
    // 1. Copy across data 
    // col1: E_l (cm-1), col2: E_u (cm-1), col3: g_l, col4: g_u, col5: Aul (1/s), col6: L_l, col7: L_u, col8: type
    E_l = line_data[0] * RC_c * RC_h_SI;
    E_u = line_data[1] * RC_c * RC_h_SI;
    nu_ul = ( line_data[1] - line_data[0] ) * RC_c;
    g_l   = int(line_data[2]);
    g_u   = int(line_data[3]);
    A_ul  = line_data[4];
    // additional data for CR modelling may be provided
    // NOTE: may be -1 indicating data is NA (see search below)
    ie_l = int(line_data[5]);
    ie_u = int(line_data[6]);
    type = int(line_data[7]);
    // Stark broadening and shift parameters
    double n = -1;
    double gamma_S0 = -1;
    if ( line_data.size()==10 ) {
        n = double(line_data[8]);
        gamma_S0 = double(line_data[9]) * nu_ul / ( 10.0 * nu2lambda(nu_ul) );       // Ang -> Hz
    }

    // 2. search for nearest energy levels for ie_l and ie_u if -1
    // -> now done by AtomicRadiator constructor
    
    // 3. Calculate and store B_lu and f_lu
    B_lu = A_ul * ( double(g_u) / double(g_l) ) * RC_c_SI * RC_c_SI / \
    		( 8.0 * M_PI * RC_h_SI ) / ( nu_ul * nu_ul * nu_ul );
    f_lu = RC_m * RC_c * RC_c * RC_c / ( 8.0 * M_PI * M_PI * RC_e * RC_e ) * \
    		( double(g_u) / double(g_l) ) / ( nu_ul * nu_ul ) * A_ul;

    // 4. Make the cluster function instance
    rcf = new RobertsClusterFunction(1,0,beta);

    // 5. Make the StarkWidthModel instance
    if( gamma_S0 < 0.0 ) {
        // Use the approximate model as no specific data has been given for this line
#       if ATOMIC_APPROX_STARK_WIDTH_METHOD==0
        double constA = 1.69e10; double constB = 2.623;
#       elif ATOMIC_APPROX_STARK_WIDTH_METHOD==1
        double constA = 9.27e07; double constB = 2.000;
#       elif ATOMIC_APPROX_STARK_WIDTH_METHOD==2
        double constA = 4.20e07; double constB = 2.000;
#       elif ATOMIC_APPROX_STARK_WIDTH_METHOD==3
        double constA = 5.00e07; double constB = 2.000;
#       endif
        swm = new ApproxStarkWidth( 0.33, constA, constB, I, E_u );
    }
    else {
        // Specific Stark parameters have been given for this line
        swm = new GriemStarkWidth( n, gamma_S0 );
    }
}

AtomicLine::
~AtomicLine()
{
    delete rcf;
    delete swm;
}

string
AtomicLine::
string()
{
    ostringstream ost;
    ost << setw(15) << E_l / ( RC_c * RC_h_SI )
        << setw(15) << E_u / ( RC_c * RC_h_SI )
        << setw(15) << nu_ul
        << setw(5) << g_l
        << setw(5) << g_u
        << setw(15) << A_ul
        << setw(10) << ie_l
        << setw(10) << ie_u
        << setw(10) << type;
	
    return ost.str();
}

void
AtomicLine::
initialise( double T, double T_e, double p, double N_hvy, double N_e, double mw_av  )
{
    double N_u = elev_u->calculate_boltz_N(T_e,g_u,E_u);
    double N_l = elev_u->calculate_boltz_N(T_e,g_l,E_l);
    
    // 1. Calculate spectral coefficients
    j_ul = N_u * RC_h_SI * nu_ul * A_ul / ( 4.0 * M_PI );	// W/m**3-sr
    double B_ul = B_lu * ( double(g_l) / double(g_u) );
    kappa_lu = ( N_l * B_lu - N_u * B_ul ) * RC_h_SI * nu_ul;
    
    // 2. Calculate line-widths
    gamma_L = calculate_lorentz_width( T, T_e, p, N_hvy, N_e, N_l*1.0e-6, mw_av );
    gamma_D = calculate_doppler_width( T );
    gamma_V = calculate_voigt_width();
    
    return;
}

void
AtomicLine::
spectral_distribution( vector<double> &nus )
{
    // One-sided line points and extent
    int    line_points = npoints;
    double line_extent = gamma_V * double(nwidths);

    // First the line-peak
    // NOTE: we sort the list later, so can do this in any order we like
    nus.push_back( nu_ul );

    // Now the points either side
    for ( int inu = 1; inu <= line_points; inu++ ) {
        double t = double(inu) / double(line_points);
        double delta_nu = line_extent * rcf->eval(t);
        nus.push_back( nu_ul + delta_nu );
        nus.push_back( nu_ul - delta_nu );
    }

    return;
}

void
AtomicLine::
calculate_spectrum( CoeffSpectra &X )
{
    // 1. Calculate desired spectral range
    int inu_start = 0;
    int inu_end = int ( X.nu.size() );
#   if LIMITED_ATOMIC_LINE_EXTENT
    double nu_lower = nu_ul - double(nwidths) * gamma_V;
    double nu_upper = nu_ul + double(nwidths) * gamma_V;
    inu_start = get_nu_index(X.nu,nu_lower,X.adaptive) + 1;
    inu_end = get_nu_index(X.nu,nu_upper,X.adaptive) + 1;
#   endif

    double nu, delta_nu, b_nu;

    // 2. Loop over spectral range
    for ( int inu=inu_start; inu<inu_end; inu++ ) {
	nu = X.nu[inu];
	delta_nu = fabs( nu - nu_ul );
	b_nu = get_voigt_point( delta_nu );
	X.j_nu[inu] += b_nu*j_ul;
	X.kappa_nu[inu] += b_nu*kappa_lu;
    }
    return;
}

double
AtomicLine::
get_lorentz_point( double delta_nu )
{
    /* Line shape as a function of delta_nu for a Lorentz profile */
    double bp_nu;
    
    bp_nu = gamma_L / M_PI / ( ( delta_nu * delta_nu ) + gamma_L * gamma_L );
    
    return bp_nu;
}

double
AtomicLine::
get_doppler_point( double delta_nu )
{
    /* Line shape as a function of delta_nu for a Doppler profile */
    double bd_nu;
    
    bd_nu = sqrt (log(2.0) / M_PI ) / gamma_D * exp ( -log(2.0) * pow( delta_nu/gamma_D, 2.0) );
    
    return bd_nu;
}

double
AtomicLine::
get_voigt_point( double delta_nu )
{
    // Ref: Whiting (1968) JQRST Vol. 8 pp 1379-1384
    double R_l = delta_nu / ( 2.0 * gamma_V );
    double R_d = gamma_L / gamma_V;
    
#   if ATOMIC_VOIGT_PROFILE_METHOD == 0
    // Accurate expression
    double tmpA = ( 1.0 - R_d ) * exp( -2.772 * R_l * R_l ) + R_d / ( 1.0 + 4.0 * R_l * R_l );
    double tmpB = 0.016 * ( 1.0 - R_d ) * R_d * ( exp( -0.4 * pow( R_l, 2.25 ) ) \
    			- 10.0 / ( 10.0 + pow( R_l, 2.25 ) ) );
    double tmpC = 2.0 * gamma_V * (1.065 + 0.447 * R_d + 0.058 * R_d * R_d);
    
    double b_nu = ( tmpA + tmpB ) / tmpC;
    
#   elif ATOMIC_VOIGT_PROFILE_METHOD == 1
    // Approximate expression
    double tmpA = ( 1.0 - R_d ) * exp( -2.772 * R_l * R_l ) + R_d / ( 1.0 + 4.0 * R_l * R_l );
    double tmpC = 2.0 * gamma_V * (1.065 + 0.447 * R_d + 0.058 * R_d * R_d);
    
    double b_nu = tmpA / tmpC;
    
#   endif
    
    return b_nu;
}

string
AtomicLine::
line_width_string( double T, double T_e, double p, double N_hvy, double N_e, double mw_av )
{
    double N_l = elev_u->calculate_boltz_N(T_e,g_l,E_l);
    
    ostringstream ost;
    ost << setprecision(12) << scientific << showpoint;
    
    double lambda_nm_ul = nu2lambda(nu_ul);
    double gamma_D_nm = lambda_nm_ul / nu_ul * calculate_doppler_width(T);
    double gamma_VW_nm = lambda_nm_ul / nu_ul * calculate_vanderwaals_width(T,p,mw_av,N_hvy);
    double gamma_R_nm = lambda_nm_ul / nu_ul * calculate_resonance_width( N_l*1.0e-6 );
    double gamma_S_nm = lambda_nm_ul / nu_ul * calculate_Stark_width(T_e,N_e);
    double gamma_N_nm = lambda_nm_ul / nu_ul * calculate_natural_width();
    ost << setw(20) << lambda_nm_ul << setw(20) << gamma_D_nm << setw(20) << gamma_VW_nm 
	<< setw(20) << gamma_R_nm   << setw(20) << gamma_S_nm << setw(20) << gamma_N_nm << endl;
    
    return ost.str();
}

double
AtomicLine::
calculate_lorentz_width( double T, double T_e, double p, double N_hvy, double N_e, double N_l, double mw_av )
{
    double gamma_R  = calculate_resonance_width(N_l);
    double gamma_VW = calculate_vanderwaals_width(T,p,mw_av,N_hvy);
    double gamma_S  = calculate_Stark_width(T_e,N_e);
    double gamma_N  = calculate_natural_width();
    double _gamma_L  = gamma_VW + gamma_R + gamma_S + gamma_N;
    
#   if DEBUG_RAD > 0
    if ( isnan(gamma_L) ) {
        cout << "  gamma_VW = " << gamma_VW
             << ", gamma_R = " << gamma_R
             << ", gamma_S = " << gamma_S
             << ", gamma_N = " << gamma_N
             << ", gamma_L = " << gamma_L << endl;
        exit( FAILURE );
    }
#   endif
    
    return _gamma_L;

}

double
AtomicLine::
calculate_resonance_width( double N_l )
{
    /* New formulation using Johnston's expression, Eq. 3.12 in PhD thesis.
       NOTE: Johnston's formula is wavelength based, this one is frequency based.
             This expression is actually the same as LORAN/Hartung/Nicolet, but 
             here the lower state population is used instead of 'n_perturbers' */
    double tmpA = 3.0 * M_PI * sqrt( double(g_l) / double(g_u) );
    double tmpB = RC_e * RC_e * f_lu / ( 2.0 * M_PI * RC_m * nu_ul );
    
    double gamma_R = tmpA * tmpB * N_l;
    
    return gamma_R;
    
}

double
AtomicLine::
calculate_vanderwaals_width( double T, double p, double mw_av, double N_hvy )
{
    /* AKA Pressure broadening by non-resonant collisions.
       From radipac6.f90 line 1626, referencing Lochte-Holtgreven (1968) pp 66-134 */
       
    double lambda_cl = nu2lambda( nu_ul ) * 1.0e1;			// Angstroms
    double gamma_VW_ang = 5.85e-30 * sqrt ( 2.0 * T / mw_av ) * N_hvy * pow( lambda_cl, 2);
    
    // Convert delta Ang to delta Hz
    double gamma_VW = gamma_VW_ang / lambda_cl * nu_ul;

    return gamma_VW;
}

double
AtomicLine::
calculate_natural_width()
{
    /* Method A: use constant value as used by Spradian07.
       Method B: use classical radiation theory expression, p.190 in Thorne et al "Spectrophysics..." */
    
#   if 0
    double gamma_N = 1.18e-5 * nu_ul / nu2lambda(nu_ul);
#   else
    double gamma_N = ( 2.0 * M_PI * RC_e_SI*RC_e_SI * nu_ul*nu_ul ) /
	          ( 3.0 * RC_eps0_SI * RC_m_SI * RC_c_SI*RC_c_SI*RC_c_SI );
#   endif
    
    return gamma_N;
}

double
AtomicLine::
calculate_doppler_width( double T )
{
    double m_s = m_w * 1000.0 / RC_Na;	// mass-per-particle in grams
    
    /* NOTE: Johnston uses T_e, but T should govern Maxwell distribution... */
    double _gamma_D = (nu_ul / RC_c) * sqrt ( (2.0 * RC_k * T * log(2.0)) / m_s);
    
    return _gamma_D;
}

double
AtomicLine::
calculate_voigt_width()
{
    double R_d, d, beta, alpha;

    d = ( gamma_L - gamma_D ) / ( gamma_L + gamma_D );
    beta = 0.023665 * exp ( 0.6 * d) + 0.0418 * exp ( -1.9 * d);
    alpha = 0.18121;
    R_d = 1.0 - alpha * ( 1.0 - d * d ) - beta * sin(M_PI * d);
    
    double _gamma_V = R_d * ( gamma_L + gamma_D);
    
    return _gamma_V;
}

