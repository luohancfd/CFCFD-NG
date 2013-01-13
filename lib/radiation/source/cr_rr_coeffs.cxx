/** \file cr_rr_coeffs.cxx
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 05-04-10: Port from lib/radiation
 *  \brief Definitions for collisional-radiative reaction rate coefficients
 *
 **/
 
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include <sstream>

#include "../../util/source/useful.h"
#include "../../nm/source/exponential_integrals.hh"
#include "../../util/source/lua_service.hh"

#include "cr_rr_coeffs.hh"
#include "radiation_constants.hh"
#include "spectral_model.hh"
#include "atomic_radiator.hh"

using namespace std;

// Gaussian 10 point parameters
const double gn10_x [] = { 0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285 };
const double gn10_w [] = { 0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443 };

ZeroRate::ZeroRate()
{
    type = "ZeroRate";
    equilibrium_flag = false;
}

double ZeroRate::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(T);
    UNUSED_VARIABLE(Q);
    
    return 0.0;
}

FromEquilibriumConstant::FromEquilibriumConstant()
{
    type = "FromEquilibriumConstant";
    equilibrium_flag = true;
}

double FromEquilibriumConstant::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(T);
    UNUSED_VARIABLE(Q);
        
    return 0.0;
}

RadGeneralisedArrhenius::RadGeneralisedArrhenius( lua_State *L )
{
    type = "RadGeneralisedArrhenius";
    equilibrium_flag = false;
    
    A = get_positive_number( L, -1, "A" );
    n = get_int( L, -1, "n" );
    T_a = get_number( L, -1, "T_a" );
}

RadGeneralisedArrhenius::
RadGeneralisedArrhenius( double A, double n, double T_a )
: A(A), n(n), T_a(T_a)
{
    type = "RadGeneralisedArrhenius";
    equilibrium_flag = false;
}

double RadGeneralisedArrhenius::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(Q);
    
    return A * pow(T,n) * exp( -T_a / T );
}

string RadGeneralisedArrhenius::get_latex_string()
{
    ostringstream tmpA, tmpB, A_oss, n_oss, T_a_oss;
    int len;
    
    // A
    tmpA << scientific << setprecision(3) << A / RC_Na;
    len = (int) tmpA.str().size();
    A_oss << "$" << tmpA.str().substr(0,4) << " \\times 10^{" << tmpA.str().substr(5,len-5) << "}$";
    
    // n
    n_oss << fixed << setprecision(2) << n;
    
    // T_a
    tmpB.clear();
    tmpB << fixed << setprecision(0) << T_a;
    len = tmpB.str().size();
    if ( len > 3 ) {
    	T_a_oss << tmpB.str().substr(0,len-3) << "," << tmpB.str().substr(len-3,3);
    }
    
    ostringstream oss;
    oss << A_oss.str() << "   &  " << n_oss.str() << "  &  " << T_a_oss.str();
    
    return oss.str();
}

// Park's version with (T/6000)**n

RadGeneralisedArrheniusPark::RadGeneralisedArrheniusPark( lua_State *L )
{
    type = "RadGeneralisedArrheniusPark";
    equilibrium_flag = false;
    
    A = get_positive_number( L, -1, "A" ) * RC_Na;
    n = get_int( L, -1, "n" );
    T_a = get_number( L, -1, "T_a" );
}

RadGeneralisedArrheniusPark::
RadGeneralisedArrheniusPark( double A, double n, double T_a )
: A(A*RC_Na), n(n), T_a(T_a)
{
    type = "RadGeneralisedArrheniusPark";
    equilibrium_flag = false;
}

double RadGeneralisedArrheniusPark::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(Q);
    
    return A * pow(T/6000.0,n) * exp( -T_a / T );
}

string RadGeneralisedArrheniusPark::get_latex_string()
{
    ostringstream tmpA, tmpB, A_oss, n_oss, T_a_oss;
    int len;
    
    // A
    tmpA << scientific << setprecision(3) << A * pow(1.0/6000.0,n);
    len = (int) tmpA.str().size();
    A_oss << "$" << tmpA.str().substr(0,4) << " \\times 10^{" << tmpA.str().substr(5,len-5) << "}$";
    
    // n
    n_oss << fixed << setprecision(2) << n;
    
    // T_a
    tmpB << fixed << setprecision(0) << T_a;
    len = tmpB.str().size();
    if ( len > 3 ) {
    	T_a_oss << tmpB.str().substr(0,len-3) << "," << tmpB.str().substr(len-3,3);
    }
    
    ostringstream oss;
    oss << A_oss.str() << "   &  " << n_oss.str() << "  &  " << T_a_oss.str();
    
    return oss.str();
}

DrawinOpticallyAllowedElectronImpactExcitation::
DrawinOpticallyAllowedElectronImpactExcitation( double E_l, double E_u )
: E_l( E_l ), E_u( E_u )
{
    type = "DrawinOpticallyAllowedElectronImpactExcitation";
}

double DrawinOpticallyAllowedElectronImpactExcitation::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(Q);
    
    // Using ABBA formulation as documented in Panesi (2007)
    double kTe = RC_k_SI * T;
    double a = ( E_u - E_l )/kTe;
    double I1a = 0.63255 * pow( a, -1.6454 ) * exp( -a );
    double alpha = 0.05;
    
    double tmpA = sqrt( ( 8.0 * kTe )/( M_PI * RC_m_SI ) );	// [m/s]
    double tmpB = 4.0 * M_PI * RC_a0_SI * RC_a0_SI * alpha;	// [m**2]
    double tmpC = RC_H_ionise_J/kTe * RC_H_ionise_J/kTe * I1a;	// [m**-6?]
    
    double K = tmpA * tmpB * tmpC;
    K *= RC_Na * 1.0e6;		// convert m**3/(particle)-s -> cm**3/mole-s
        
    return K;
}

DrawinOpticallyForbiddenElectronImpactExcitation::
DrawinOpticallyForbiddenElectronImpactExcitation( double E_l, double E_u )
: E_l( E_l ), E_u( E_u )
{
    type = "DrawinOpticallyForbiddenElectronImpactExcitation";
}

double DrawinOpticallyForbiddenElectronImpactExcitation::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(Q);
    
    // Using ABBA formulation as documented in Panesi (2007)
    double kTe = RC_k_SI * T;
    double a = ( E_u - E_l )/kTe;
    double I2a = 0.23933 * pow( a, -1.4933 ) * exp( -a );
    double alpha = 0.05;
    
    double tmpA = sqrt( ( 8.0 * kTe )/( M_PI * RC_m_SI ) );	// [m/s]
    double tmpB = 4.0 * M_PI * RC_a0_SI * RC_a0_SI * alpha;	// [m**2]
    double tmpC = a * a * I2a;	                                // [m**-6?]
    
    double K = tmpA * tmpB * tmpC;
    K *= RC_Na * 1.0e6;		// convert m**3/(particle)-s -> cm**3/mole-s
        
    return K;
}

GryzinskiElectronImpactExcitation::
GryzinskiElectronImpactExcitation( Radiator * rad, ElecLev * elev_l, ElecLev * elev_u )
: E_l( elev_l->E ), E_u( elev_u->E ),
  I( rad->I )
{
    type = "GryzinskiElectronImpactExcitation";
    
    // NOTE: storing all energy as cm**-1
    if ( elev_u->i<(rad->nlevs-1) ) {
	// The u+1 level exists
	E_up1 = rad->get_elev_pointer( elev_u->i + 1 )->E;
    }
    else {
	// The u+1 level is the ionised state - use I
	E_up1 = I;
    }
}

GryzinskiElectronImpactExcitation::
GryzinskiElectronImpactExcitation( double E_l, double E_u, double E_up1, double I )
: E_l( E_l ), E_u( E_u ), E_up1( E_up1 ), I( I ) {}

double GryzinskiElectronImpactExcitation::get_rate( double T, Gas_data &Q )
{
    double tmpA = 8.0 * M_PI / sqrt( RC_m ) * pow( (2.0*M_PI*RC_m*RC_k*T), -1.5 );
    // Try using Park derivation:
    // double tmpA = 5.45 * pow( T, -1.5 ) / ( RC_E0*RC_E0 * M_PI * RC_a0*RC_a0 );
    
    // Compute integral's
    double integral_lu = compute_integral( T, E_l, E_u );
    double integral_lup1 = compute_integral( T, E_l, E_up1 );
    
    double K = tmpA * ( integral_lu - integral_lup1 ) * pow( RC_h_SI * RC_c, 3) * 1.0e-4;
    K *= RC_Na;					// Convert cm**3/(particle)-s -> cm**3/mole-s
    
    // cout << "tmpA = " << tmpA << ", integral_lu = " << integral_lu << ", integral_lup1 = " << integral_lup1 << endl;
    
    return K;
}

double GryzinskiElectronImpactExcitation::sigma( double E, double E_i, double E_j )
{
    // NOTE: energy units need to be in cm-1 according to Johnston notation
    //       But only leading delta_E_ij**2 term affects dimensions
    double delta_E_ij = E_j - E_i;
    double tmpA = 4.2484e-6 / ( delta_E_ij * delta_E_ij ) * pow( E / ( I - E_i + E ), 1.5 );
    double tmpB = (2.0/3.0) * ( ( I - E_i ) / E + delta_E_ij / E * ( 1.0 - ( I - E_i ) / E ) - ( delta_E_ij / E ) * ( delta_E_ij / E ) );
    double tmpC = 1.0;
    if ( ( delta_E_ij + I - E_i ) > E ) {
        tmpC = sqrt( ( 1.0 + delta_E_ij / ( I - E_i ) ) * ( 1.0 - delta_E_ij / E ) );
    }
    
    double sigma_cm2 = tmpA * tmpB * tmpC;		// cm**2 / EU**2	[EU = energy unit]
    
    // cout << "E = " << E << ", E_i = " << E_i << ", E_j = " << E_j << ", E_ionise = " << I << ", tmpA = " << tmpA << ", tmpB = " << tmpB << ", tmpC = " << tmpC << ", sigma = " << sigma_cm2 << endl;
    
    return sigma_cm2;
}

/// \brief Integrand i,j with change-of-variables to map onto finite domain

double GryzinskiElectronImpactExcitation::compute_integral( double T, double E_i, double E_j )
{
    // Compute integral via mapping onto finite domain and applying Gaussian 10 point method
    
    double integral = 0.0;
    double x, dx, integrand_l, integrand_u, E;
    for ( int j=1; j<=5; ++j ) {
	dx = gn10_x[j];
	// Firstly lower point
	x = -dx;
	E = r( x, E_i, E_j );
	integrand_l = drdx( x, E_i, E_j ) * sigma(E,E_i,E_j) * exp( - E / ( RC_k_SI * T ) ) * E;
	// cout << "x = " << x << ", E = " << E << ", drdx = " << drdx( x, E_i, E_j ) << ", sigma = " << sigma(E,E_i_,E_j_) << ", integrand = " << integrand_l << endl;
	// Now upper point
	x = dx;
	E = r( x, E_i, E_j );
	integrand_u = drdx( x, E_i, E_j ) * sigma(E,E_i,E_j) * exp( - E / ( RC_k_SI * T ) ) * E;
	// combine and scale by weighting term
	integral += gn10_w[j] * ( integrand_l + integrand_u );
    }
	
    return integral;
}

/// \brief Change-of-variable mapping function

double GryzinskiElectronImpactExcitation::r( double x, double E_i, double E_j )
{
    // Map x from [-1,1] domain to energy [E_min,inf) domain
    double E_min = E_j - E_i;
    return ( 2.0 * E_min ) / ( 1.0 - x );
}

/// \brief Change-of-variable mapping function derivative

double GryzinskiElectronImpactExcitation::drdx( double x, double E_i, double E_j )
{
    // Derivative of r wrt x
    double E_min = E_j - E_i;
    return ( 2.0 * E_min ) / ( ( 1.0 - x ) * ( 1.0 - x ) );
}

DrawinElectronImpactIonization::
DrawinElectronImpactIonization( double E_l, double I )
: E_l( E_l ), I( I )
{
    type = "DrawinElectronImpactIonization";
}

double DrawinElectronImpactIonization::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(Q);
    
    // ABBA formulation as documented in Panesi (2007)
    double kTe = RC_k_SI * T;
    double a = ( I - E_l )/kTe;
    double I1a = 0.63255 * pow( a, -1.6454 ) * exp( -a );
    double alpha = 1.0;
    
    double tmpA = sqrt( ( 8.0 * kTe )/( M_PI * RC_m_SI ) );	// [m/s]
    double tmpB = 4.0 * M_PI * RC_a0_SI * RC_a0_SI * alpha;	// [m**2]
    double tmpC = RC_H_ionise_J/kTe * RC_H_ionise_J/kTe * I1a;	// [ND?]
    
    double K = tmpA * tmpB * tmpC;
    K *= RC_Na * 1.0e6;		// convert m**3/(particle)-s -> cm**3/mole-s
        
    return K;
}

CJDrawinElectronImpactIonization::
CJDrawinElectronImpactIonization( double E_l, double I )
: E_l( E_l ), I( I )
{
    type = "CJDrawinElectronImpactIonization";
}

double CJDrawinElectronImpactIonization::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(Q);
    
    // Johnston formulation from his PhD thesis
    double y = ( I - E_l ) / RC_k_SI / T;
    double psi = exp( -y )  / ( 1.0 + y ) * ( 1.0 / ( 20.0 + y ) + log( 1.25 * ( 1.0 + 1.0 / y ) ) );
    double zeta = 1.0;
    if ( E_l / RC_h_SI / RC_c < 100.0 ) zeta = 3.0;	// ~ test for ground state
    double K = 1.46e-10 * sqrt(T) *  RC_H_ionise_J * RC_H_ionise_J / ( I - E_l ) / ( I - E_l ) * zeta * psi;
    K *= RC_Na;		// convert cm**3/(particle)-s -> cm**3/mole-s
        
    return K;
}

OpticallyThinExponentialDecay::OpticallyThinExponentialDecay( lua_State * L ) 
{
    type = "OpticallyThinExponentialDecay";
    equilibrium_flag = false;
    
    tau = get_positive_number( L, -1, "tau" );
    lambda = 1.0;
}

OpticallyThinExponentialDecay::OpticallyThinExponentialDecay( double tau )
: tau( tau ), lambda( 1.0 ) 
{
    type = "OpticallyThinExponentialDecay";
    equilibrium_flag = false;
}

double OpticallyThinExponentialDecay::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(T);
    UNUSED_VARIABLE(Q);
    
    return lambda / tau;
}

string
OpticallyThinExponentialDecay::get_latex_string()
{
    ostringstream tmpA, A_oss;
    int len;
    
    // 1/tau
    tmpA << scientific << setprecision(3) << 1/tau;
    len = (int) tmpA.str().size();
    A_oss << "$" << tmpA.str().substr(0,4) << " \\times 10^{" << tmpA.str().substr(5,len-5) << "}$";

    ostringstream oss;
    oss << A_oss.str() << "   &  ";
    
    return oss.str();
}

OpticallyVariableExponentialDecay::
OpticallyVariableExponentialDecay( lua_State * L )
{
    cout << "OpticallyVariableExponentialDecay::OpticallyVariableExponentialDecay()" << endl
	 << "Not implemented yet." << endl;
    exit( BAD_INPUT_ERROR );
}

OpticallyVariableExponentialDecay::
OpticallyVariableExponentialDecay( double tau, double wavel, double wavel_switch, 
    				   double lambda_lower, double lambda_upper )
: tau( tau )
{
    type = "OpticallyVariableExponentialDecay";
    equilibrium_flag = false;
    
    // Determine which escape factor to use based on the transition wavelength
    if ( wavel < wavel_switch ) lambda = lambda_lower;
    else lambda = lambda_upper;
}

double OpticallyVariableExponentialDecay::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(T);
    UNUSED_VARIABLE(Q);
    
    return lambda / tau;
}

CurveFitExponentialDecay::CurveFitExponentialDecay( double tau )
: tau( tau )
{
    type = "CurveFitExponentialDecay";
    equilibrium_flag = false;
}

double CurveFitExponentialDecay::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(T);
    UNUSED_VARIABLE(Q);
    
    // Ref: Johnston (2006) PhD
    
    double lambda;
    if ( T < 4477.0 ) lambda = 0.0;
    else lambda = -4.6079e-8 * T * T + 8.052e-4 * T - 2.6814;
    
    return lambda / tau;
}

FrostNitrogenElectronImpactExcitation::
FrostNitrogenElectronImpactExcitation( lua_State * L, ElecLev * elev_l, ElecLev * elev_u )
: delta_E( elev_u->E - elev_l->E ), g_l( elev_l->g )
{
    type = "FrostNitrogenElectronImpactExcitation";
    
    lua_getfield(L, -1, "temperatures");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	oss << "Error locating 'temperatures' table" << endl;
	input_error(oss);
    }
    vector<double> T_vec;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	T_vec.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop temperatures table
    
    // get reaction field
    ostringstream oss;
    oss << "reac_EIE_ilev" << elev_l->i << "_to_ilev" << elev_u->i;
    lua_getfield(L,-1,oss.str().c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	ost << "Error locating " << oss.str() << " table" << endl;
	input_error(ost);
    }
    vector<double> data;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop reaction table
    
    if ( data.size() != T_vec.size( )) {
	ostringstream oss;
	oss << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	oss << "data.size() = " << data.size() << " is not equal to T_vec.size() = " << T_vec.size() + 2 << endl;
	input_error(oss);
    }
    
    // Create spline from vector of 'coords'
    // x: Electronic temperature in K
    // y: Dimensionless collision parameter gamma
    // NOTE: converting Te from eV to K
    vector<Vector3> points;
    for ( size_t i=0; i<T_vec.size(); ++i )
    	points.push_back( Vector3( T_vec[i] * RC_e_SI / RC_k_SI, data[i] ) );
    
    gamma_fit = new Spline( points );
}

FrostNitrogenElectronImpactExcitation::
FrostNitrogenElectronImpactExcitation( string fname, ElecLev * elev_l, ElecLev * elev_u )
: delta_E( elev_u->E - elev_l->E ), g_l( elev_l->g )
{
    type = "FrostNitrogenElectronImpactExcitation";
    
    lua_State *L = initialise_radiation_lua_State();
    if( do_gzfile(L, fname) != 0 ) {
	ostringstream ost;
	ost << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation():\n";
	ost << "Error in input file: " << fname << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "N");
    lua_getfield(L, -1, "QSS_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	ost << "Error locating 'QSS_model' table" << endl;
	input_error(ost);
    }
    
    // Search for this transition
    lua_getfield(L, -1, "frost_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	oss << "Error locating 'frost_data' table" << endl;
	input_error(oss);
    }
    
    lua_getfield(L, -1, "temperatures");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	oss << "Error locating 'temperatures' table" << endl;
	input_error(oss);
    }
    vector<double> T_vec;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	T_vec.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop temperatures table
    
    // get reaction field
    ostringstream oss;
    oss << "reac_EIE_ilev" << elev_l->i << "_to_ilev" << elev_u->i;
    lua_getfield(L,-1,oss.str().c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	ost << "Error locating " << oss.str() << " table" << endl;
	input_error(ost);
    }
    vector<double> data;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop reaction table
    
    if ( data.size() != T_vec.size() ) {
	ostringstream oss;
	oss << "FrostNitrogenElectronImpactExcitation::FrostNitrogenElectronImpactExcitation()\n";
	oss << "data.size() = " << data.size() << " is not equal to T_vec.size() = " << T_vec.size() << endl;
	input_error(oss);
    }
    
    // Create spline from vector of 'coords'
    // x: Electronic temperature in K
    // y: Dimensionless collision parameter gamma
    // NOTE: converting Te from eV to K
    vector<Vector3> points;
    for ( size_t i=0; i<T_vec.size(); ++i )
    	points.push_back( Vector3( T_vec[i] * RC_e_SI / RC_k_SI, data[i] ) );
    
    gamma_fit = new Spline( points );
    
    lua_pop(L,1); 	// pop frost_data
    lua_pop(L,1);	// pop QSS_model
    lua_pop(L,1);	// pop "N"
}

FrostNitrogenElectronImpactExcitation::
~FrostNitrogenElectronImpactExcitation()
{
    delete gamma_fit;
}

double FrostNitrogenElectronImpactExcitation::get_rate( double T, Gas_data &Q )
{
    double tmpA = 2.0 * sqrt( M_PI ) * RC_alpha * RC_c * RC_a0 * RC_a0;				// cm**3/s
    double tmpB = sqrt( RC_H_ionise_J / ( RC_k_SI * T ) ) * gamma_fit->eval_from_x( T ).y;	// ND
    double tmpC = exp( - delta_E / ( RC_k_SI * T ) ) / double(g_l);				// ND
    
    double K = tmpA * tmpB * tmpC * RC_Na;							// Convert cm**3/(particle)-s -> cm**3/mole-s
    
    return K;
}

ZatsarinnyTayalOxygenElectronImpactExcitation::
ZatsarinnyTayalOxygenElectronImpactExcitation( lua_State * L, ElecLev * elev_l, ElecLev * elev_u )
: delta_E( elev_u->E - elev_l->E ), g_l( elev_l->g )
{
    type = "ZatsarinnyTayalOxygenElectronImpactExcitation";
    
    // Get temperature vector
    lua_getfield(L, -1, "temperatures");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	oss << "Error locating 'temperatures' table" << endl;
	input_error(oss);
    }
    vector<double> T_vec;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	T_vec.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop temperatures table
    
    // get reaction field
    ostringstream oss;
    oss << "reac_EIE_ilev" << elev_l->i << "_to_ilev" << elev_u->i;
    lua_getfield(L,-1,oss.str().c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	ost << "Error locating " << oss.str() << " table" << endl;
	input_error(ost);
    }
    vector<double> data;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop reaction table
    
    if ( data.size() != T_vec.size( )) {
	ostringstream oss;
	oss << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	oss << "data.size() = " << data.size() << " is not equal to T_vec.size() = " << T_vec.size() << endl;
	input_error(oss);
    }
    
    // Create spline from vector of 'coords'
    // x: Electronic temperature in K
    // y: Dimensionless collision parameter gamma
    // NOTE: converting Te from eV to K
    vector<Vector3> points;
    for ( size_t i=0; i<T_vec.size(); ++i )
    	points.push_back( Vector3( T_vec[i], data[i] ) );
    
    gamma_fit = new Spline( points );
}

ZatsarinnyTayalOxygenElectronImpactExcitation::
ZatsarinnyTayalOxygenElectronImpactExcitation( string fname, ElecLev * elev_l, ElecLev * elev_u )
: delta_E( elev_u->E - elev_l->E ), g_l( elev_l->g )
{
    type = "ZatsarinnyTayalOxygenElectronImpactExcitation";
    
    lua_State *L = initialise_radiation_lua_State();
    if( do_gzfile(L, fname) != 0 ) {
	ostringstream ost;
	ost << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation():\n";
	ost << "Error in input file: " << fname << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "O");
    lua_getfield(L, -1, "QSS_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	ost << "Error locating 'QSS_model' table" << endl;
	input_error(ost);
    }
    
    // Search for this transition
    lua_getfield(L, -1, "zatsarinny_tayal_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	oss << "Error locating 'zatsarinny_tayal_data' table" << endl;
	input_error(oss);
    }

    // Get temperature vector
    lua_getfield(L, -1, "temperatures");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	oss << "Error locating 'temperatures' table" << endl;
	input_error(oss);
    }
    vector<double> T_vec;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	T_vec.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop temperatures table
    
    // get reaction field
    ostringstream oss;
    oss << "reac_EIE_ilev" << elev_l->i << "_to_ilev" << elev_u->i;
    lua_getfield(L,-1,oss.str().c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	ost << "Error locating " << oss.str() << " table" << endl;
	input_error(ost);
    }
    vector<double> data;
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L, 1);	// pop reaction table
    
    if ( data.size() != T_vec.size( )) {
	ostringstream oss;
	oss << "ZatsarinnyTayalOxygenElectronImpactExcitation::ZatsarinnyTayalOxygenElectronImpactExcitation()\n";
	oss << "data.size() = " << data.size() << " is not equal to T_vec.size() = " << T_vec.size() << endl;
	input_error(oss);
    }
    
    // Create spline from vector of 'coords'
    // x: Electronic temperature in K
    // y: Dimensionless collision parameter gamma
    // NOTE: converting Te from eV to K
    vector<Vector3> points;
    for ( size_t i=0; i<T_vec.size(); ++i )
    	points.push_back( Vector3( T_vec[i], data[i] ) );
    
    gamma_fit = new Spline( points );
    
    lua_pop(L,1); 	// pop zatsarinny_tayal_data
    lua_pop(L,1);	// pop QSS_model
    lua_pop(L,1);	// pop "O"
}

ZatsarinnyTayalOxygenElectronImpactExcitation::
~ZatsarinnyTayalOxygenElectronImpactExcitation()
{
    delete gamma_fit;
}

double ZatsarinnyTayalOxygenElectronImpactExcitation::get_rate( double T, Gas_data &Q )
{
    double tmpA = 8.629e-6 / double(g_l) / sqrt( T );		// cm**3 / s
    double gamma = gamma_fit->eval_from_x( T ).y;		// ND
    double tmpB = exp( - delta_E / RC_k_SI / T );		// ND
    
    
    double K = tmpA * gamma * tmpB * RC_Na;			// Convert cm**3/(particle)-s -> cm**3/mole-s
    
    return K;
}

SunoKatoCarbonElectronImpactExcitation::
SunoKatoCarbonElectronImpactExcitation( lua_State * L, ElecLev * elev_l, ElecLev * elev_u )
: delta_E( elev_u->E - elev_l->E ), g_l( elev_l->g )
{
    type = "SunoKatoCarbonElectronImpactExcitation";
    
    // Get parameters from the lua file
    V_if = get_number(L,-1,"V_if");
    eqn = get_int(L,-1,"eqn");
    A = get_number(L,-1,"A");
    B = get_number(L,-1,"B");
    C = get_number(L,-1,"C");
    D = get_number(L,-1,"D");
    E = get_number(L,-1,"E");
    F = get_number(L,-1,"F");
}

SunoKatoCarbonElectronImpactExcitation::
SunoKatoCarbonElectronImpactExcitation( string fname, ElecLev * elev_l, ElecLev * elev_u )
: delta_E( elev_u->E - elev_l->E ), g_l( elev_l->g )
{
    type = "SunoKatoCarbonElectronImpactExcitation";
    
    // Initialise the lua_State
    lua_State *L = initialise_radiation_lua_State();
    if( do_gzfile(L, fname) != 0 ) {
	ostringstream ost;
	ost << "SunoKatoCarbonElectronImpactExcitation::SunoKatoCarbonElectronImpactExcitation():\n";
	ost << "Error in input file: " << fname << endl;
	input_error(ost);
    }
    
    // Get to the suno and kato data section
    lua_getglobal(L, "C");
    lua_getfield(L, -1, "QSS_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "SunoKatoCarbonElectronImpactExcitation::SunoKatoCarbonElectronImpactExcitation()\n";
	ost << "Error locating 'QSS_model' table" << endl;
	input_error(ost);
    }
    lua_getfield(L, -1, "suno_kato_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "SunoKatoCarbonElectronImpactExcitation::SunoKatoCarbonElectronImpactExcitation()\n";
	oss << "Error locating 'suno_kato_data' table" << endl;
	input_error(oss);
    }
    
    // Search for this transition
    ostringstream oss;
    oss << "reac_EIE_ilev" << elev_l->i << "_to_ilev" << elev_u->i;
    lua_getfield(L,-1,oss.str().c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "SunoKatoCarbonElectronImpactExcitation::SunoKatoCarbonElectronImpactExcitation()\n";
	ost << "Error locating " << oss.str() << " table" << endl;
	input_error(ost);
    }
    
    // Get parameters from the lua file
    V_if = get_number(L,-1,"V_if");
    eqn = get_int(L,-1,"eqn");
    A = get_number(L,-1,"A");
    B = get_number(L,-1,"B");
    C = get_number(L,-1,"C");
    D = get_number(L,-1,"D");
    E = get_number(L,-1,"E");
    F = get_number(L,-1,"F");
    
    lua_pop(L,1);	// pop suno_and_kato_data
    lua_pop(L,1);	// pop QSS model
    lua_pop(L,1);	// pop 'C'
}

SunoKatoCarbonElectronImpactExcitation::
~SunoKatoCarbonElectronImpactExcitation()
{}

double SunoKatoCarbonElectronImpactExcitation::get_rate( double T, Gas_data &Q )
{
    // 0. Calculate reduced temperature y
    double T_eV = T * RC_k_SI / RC_e_SI;
    double y = V_if / T_eV;
    
    // 1. Calculate gamma
    double gamma = compute_gamma( y );
    
    // 2. Calculate the reaction rate coefficient
    double R = 8.010e-8 / double(g_l) / sqrt( T_eV ) * exp( -y ) * gamma;
    double K = fabs(R * RC_Na);	// Convert cm3/s -> cm3/s/mol and ensure positivity
    
    // 3. Return result in cm**3/mole-s
    return K;
}

double SunoKatoCarbonElectronImpactExcitation::compute_gamma( double y )
{
    // NOTE: the following function makes use of the Gaussian 10 point quadrature
    //       points and weights declared as gn10_x and gn10_w respectively
    if ( eqn==10 || eqn==11 ) {
    	// 1. Those transitions using eqn's 10 and 11 compute gamma directly
    	// 1a. Compute E_1(y)
    	double E_1 = 0.0;
    	double x, dx, integrand_l, integrand_u, t;
    	for ( int j=1; j<=5; ++j ) {
    	    dx = gn10_x[j];
    	    // Firstly the lower point
    	    x = -dx;
    	    t = eval_t_from_x( x, y );
    	    integrand_l = dtdx( t, y ) * exp( -t ) / t;
    	    // Now upper point
    	    x = dx;
    	    t = eval_t_from_x( x, y );
    	    integrand_u = dtdx( t, y ) * exp( -t ) / t;
    	    // combine and scale by weighting term
    	    E_1 += gn10_w[j] * ( integrand_l + integrand_u );
    	}
    	// 1b. Evaluate the gamma curve fits
    	if ( eqn==10 ) {
    	    return y * ( ( A/y + C ) + 0.5 * D * ( 1 - y ) + \
    	    	exp(y) * E_1 * ( B - C*y + 0.5 * D * y * y + E/y ) );
    	}
    	else {
    	    double exp_mF = exp( -F );
    	    return A * y * ( 1.0 - exp(y) * E_1 * y ) + \
    	        ( B * exp_mF / ( F + y ) + \
    	    	  C * pow( exp_mF, 2 ) / ( 2.0 * F + y ) + \
    	    	  D * pow( exp_mF, 3 ) / ( 3.0 * F + y ) + \
    	    	  E * pow( exp_mF, 4 ) / ( 4.0 * F + y ) ) * y;
    	}
    }
    else if ( eqn==6 || eqn==7 ) {
    	// 2. Those transitions using eqn's 6 and 7 compute gamma via omega
    	double integral = 0.0;
    	double x, dx, integrand_l, integrand_u, X, omega;
    	for ( int j=1; j<=5; ++j ) {
    	    dx = gn10_x[j];
    	    // Firstly the lower point
    	    x = -dx;
    	    X = eval_X_from_x( x, y );
    	    if ( eqn==6 ) {
    	    	omega = A + B / X + C / ( X * X ) + D / ( X * X * X ) + E * log( X );
    	    }
    	    else {
    	    	double exp_mF = exp( -F );
    	    	omega = A / ( X * X ) + B * exp_mF + C * pow( exp_mF, 2 ) + \
    	    		D * pow( exp_mF, 3 ) + E * pow( exp_mF, 4 );
    	    }
    	    integrand_l = dXdx( X, y ) * omega * exp( -y * X );
    	    // Now upper point
    	    x = dx;
    	    X = eval_X_from_x( x, y );
    	    if ( eqn==6 ) {
    	    	omega = A + B / X + C / ( X * X ) + D / ( X * X * X ) + E * log( X );
    	    }
    	    else {
    	    	double exp_mF = exp( -F );
    	    	omega = A / ( X * X ) + B * exp_mF + C * pow( exp_mF, 2 ) + \
    	    		D * pow( exp_mF, 3 ) + E * pow( exp_mF, 4 );
    	    }
    	    integrand_u = dXdx( X, y ) * omega * exp( -y * X );
    	    // combine and scale by weighting term
    	    integral += gn10_w[j] * ( integrand_l + integrand_u );
    	}
    	return integral * y * exp( y );
    }
    else {
	cout << "SunoKatoCarbonElectronImpactExcitation::compute_gamma()" << endl
	     << "The given eqn index: " << eqn << " is not known." << endl
	     << "Exiting program." << endl;
	exit( BAD_INPUT_ERROR );
    }
}

double SunoKatoCarbonElectronImpactExcitation::eval_t_from_x( double x, double y )
{
    // Map x from [-1,1] domain to [y,inf) domain of t
    return - 2.0 * y / ( x - 1.0 );
}

double SunoKatoCarbonElectronImpactExcitation::dtdx( double t, double y )
{
    // compute derivative of t wrt x
    return 1.0 / ( - 2.0 * y * log( t ) );
}

double SunoKatoCarbonElectronImpactExcitation::eval_X_from_x( double x, double y )
{
    // Map x from [-1,1] domain to [1,inf) domain of X
    return - 2.0 / ( x - 1.0 );
}

double SunoKatoCarbonElectronImpactExcitation::dXdx( double X, double y )
{
    // compute derivative of X wrt x
    return 1.0 / ( - 2.0 * log( X ) );
}

SunoKatoCarbonElectronImpactIonization::
SunoKatoCarbonElectronImpactIonization( lua_State * L, ElecLev * elev, double I )
: delta_E( I - elev->E ), g_l( elev->g )
{
    type = "SunoKatoCarbonElectronImpactIonization";
    
    // Get parameters from the lua file
    // double I_eV = get_number(L,-1,"I");
    // cout << "I_eV[given] = " << I_eV << ", I_ev[calc] = " << delta_E / RC_e_SI << endl;
    V_if = delta_E / RC_e_SI;
    A1 = get_number(L,-1,"A1");
    A2 = get_number(L,-1,"A2");
    A3 = get_number(L,-1,"A3");
    A4 = get_number(L,-1,"A4");
    A5 = get_number(L,-1,"A5");
}

SunoKatoCarbonElectronImpactIonization::
SunoKatoCarbonElectronImpactIonization( string fname, ElecLev * elev, double I )
: delta_E( I - elev->E ), g_l( elev->g )
{
    type = "SunoKatoCarbonElectronImpactIonization";
    
    // Initialise the lua_State
    lua_State *L = initialise_radiation_lua_State();
    if( do_gzfile(L, fname) != 0 ) {
	ostringstream ost;
	ost << "SunoKatoCarbonElectronImpactIonization::SunoKatoCarbonElectronImpactIonization():\n";
	ost << "Error in input file: " << fname << endl;
	input_error(ost);
    }
    
    // Get to the suno and kato data section
    lua_getglobal(L, "C");
    lua_getfield(L, -1, "QSS_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "SunoKatoCarbonElectronImpactIonization::SunoKatoCarbonElectronImpactIonization()\n";
	ost << "Error locating 'QSS_model' table" << endl;
	input_error(ost);
    }
    lua_getfield(L, -1, "suno_kato_data");
    if ( !lua_istable(L, -1) ) {
	ostringstream oss;
	oss << "SunoKatoCarbonElectronImpactIonization::SunoKatoCarbonElectronImpactIonization()\n";
	oss << "Error locating 'suno_kato_data' table" << endl;
	input_error(oss);
    }
    
    // Search for this transition
    ostringstream oss;
    oss << "reac_EII_ilev" << elev->i;
    lua_getfield(L,-1,oss.str().c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "SunoKatoCarbonElectronImpactIonization::SunoKatoCarbonElectronImpactIonization()\n";
	ost << "Error locating " << oss.str() << " table" << endl;
	input_error(ost);
    }
    
    // Get parameters from the lua file
    // double I_eV = get_number(L,-1,"I");
    // cout << "I_eV[given] = " << I_eV << ", I_ev[calc] = " << delta_E / RC_e_SI << endl;
    V_if = delta_E / RC_e_SI;
    A1 = get_number(L,-1,"A1");
    A2 = get_number(L,-1,"A2");
    A3 = get_number(L,-1,"A3");
    A4 = get_number(L,-1,"A4");
    A5 = get_number(L,-1,"A5");
    
    lua_pop(L,1);	// pop suno_and_kato_data
    lua_pop(L,1);	// pop QSS model
    lua_pop(L,1);	// pop 'C'
}

double SunoKatoCarbonElectronImpactIonization::get_rate( double T, Gas_data &Q )
{         
    // 0. Calculate reduced temperature y
    double T_eV = T * RC_k_SI / RC_e_SI;
    double y = V_if / T_eV;
    
    // 1. Calculate gamma
    double gamma = compute_gamma( y );
    
    // 2. Calculate the reaction rate coefficient
    double R = 8.010e-8 / double(g_l) / sqrt( T_eV ) * exp( -y ) * gamma;
    double K = fabs(R * RC_Na);	// Convert cm3/s -> cm3/s/mol and ensure positivity
    
    // 3. Return result in cm**3/mole-s
    return K;
}

double SunoKatoCarbonElectronImpactIonization::compute_omega( double X )
{
    double E_eV = X * delta_E / RC_e_SI;
    double I_eV = delta_E / RC_e_SI;
    
    double tmp = A1 * log( X );
    tmp += A2 * ( 1.0 - X );
    tmp += A3 * pow( 1.0 - X, 2 );
    tmp += A4 * pow( 1.0 - X, 3 );
    tmp += A5 * pow( 1.0 - X, 4 );

    // Convert the cross-section to the collision strength via eqn. 3
    double sigma = 1.0e-13 / ( I_eV * E_eV ) * tmp;
    return sigma * g_l * V_if * X / 1.1969e-15;
}

double SunoKatoCarbonElectronImpactIonization::compute_gamma( double y )
{
    // NOTE: the following function makes use of the Gaussian 10 point quadrature
    //       points and weights declared as gn10_x and gn10_w respectively
    
    // Use eqn. 9 to compute gamma from omega
    double integral = 0.0;
    double x, dx, integrand_l, integrand_u, X;
    for ( int j=1; j<=5; ++j ) {
	dx = gn10_x[j];
	// Firstly the lower point
	x = -dx;
	X = eval_X_from_x( x, y );
	integrand_l = dXdx( X, y ) * compute_omega(X) * exp( -y * X );
	// Now upper point
	x = dx;
	X = eval_X_from_x( x, y );
	integrand_u = dXdx( X, y ) * compute_omega(X) * exp( -y * X );
	// combine and scale by weighting term
	integral += gn10_w[j] * ( integrand_l + integrand_u );
    }
    
    return y * exp( y ) * integral;
}

double SunoKatoCarbonElectronImpactIonization::eval_t_from_x( double x, double y )
{
    // Map x from [-1,1] domain to [y,inf) domain of t
    return - 2.0 * y / ( x - 1.0 );
}

double SunoKatoCarbonElectronImpactIonization::dtdx( double t, double y )
{
    // compute derivative of t wrt x
    return 1.0 / ( - 2.0 * y * log( t ) );
}

double SunoKatoCarbonElectronImpactIonization::eval_X_from_x( double x, double y )
{
    // Map x from [-1,1] domain to [1,inf) domain of X
    return - 2.0 / ( x - 1.0 );
}

double SunoKatoCarbonElectronImpactIonization::dXdx( double X, double y )
{
    // compute derivative of X wrt x
    return 1.0 / ( - 2.0 * log( X ) );
}

KuncSoonElectronImpactIonization::
KuncSoonElectronImpactIonization( lua_State * L, ElecLev * elev, double I )
: E_l( elev->get_E() ), I( I )
{
    type = "KuncSoonElectronImpactIonization";
    
    // Get parameters from the lua file
    A = get_number(L,-1,"A");
    chi = get_number(L,-1,"chi");
    Q_val = get_number(L,-1,"Q");
    
    // NOTE: assuming the level belongs to an atomic radiator
    l = double( dynamic_cast<AtomicElecLev*>(elev)->get_l() );
}

KuncSoonElectronImpactIonization::
KuncSoonElectronImpactIonization( double A, double xi, double Q, ElecLev * elev, double I )
: E_l( elev->get_E() ), I( I ), A( A ), chi( chi ), Q_val( Q )
{
    type = "KuncSoonElectronImpactIonization";
    
    // NOTE: assuming the level belongs to an atomic radiator
    l = double( dynamic_cast<AtomicElecLev*>(elev)->get_l() );
}

double KuncSoonElectronImpactIonization::get_rate( double T, Gas_data &Q )
{
    UNUSED_VARIABLE(Q);

    // 1. Reduced temperature
    double beta = ( I - E_l ) / RC_k_SI / T;
    
    // 2. G parameter
    double G = sqrt( beta / ( beta + 1.0 ) ) * A / ( beta + chi );
    
    // 3. Rate coefficient
    double S = 1.0e-8 * pow( RC_H_ionise_J / ( I - E_l ), 1.5 ) * Q_val / ( 2.0 * l + 1.0 ) * exp( - beta ) * G;
    
    // 4. Convert cm**3 / s -> cm**3 / mole-s
    double K = S * RC_Na;
        
    return K;
}

CR_ReactionRateCoefficient * 
create_explicit_rate_coeff( lua_State * L, string direction )
{
    CR_ReactionRateCoefficient * rate_coeff = 0;
    
    string rate_coeff_model = "blank", parameter_field;
    if ( direction=="forward" ) {
    	rate_coeff_model = get_string( L, -1, "forward_rate_coeff_model" );
    	parameter_field = "k_f";
    }
    else if ( direction=="backward" ) {
    	rate_coeff_model = get_string( L, -1, "backward_rate_coeff_model" );
    	parameter_field = "k_b";
    }
    else {
    	ostringstream oss;
    	oss << "create_explicit_rate_coeff()" << endl
    	    << "direction: " << direction << " not recognised." << endl;
    	input_error( oss );
    }
    
    if ( rate_coeff_model=="ZeroRate" ) {
    	rate_coeff = new ZeroRate();
    }
    else if ( rate_coeff_model=="FromEquilibriumConstant" ) {
    	rate_coeff = new FromEquilibriumConstant();
    }
    else if ( rate_coeff_model=="GeneralisedArrhenius" ) {
    	lua_getfield( L, -1, parameter_field.c_str() );
    	rate_coeff = new RadGeneralisedArrhenius( L );
    	lua_pop( L, 1);		// pop parameter_field
    }
    else if ( rate_coeff_model=="GeneralisedArrheniusPark" ) { 
    	lua_getfield( L, -1, parameter_field.c_str() );
    	rate_coeff = new RadGeneralisedArrheniusPark( L );
    	lua_pop( L, 1);		// pop parameter_field
    }
    else if ( rate_coeff_model=="OpticallyThinExponentialDecay" ) {
    	lua_getfield( L, -1, parameter_field.c_str() );
    	rate_coeff = new OpticallyThinExponentialDecay( L );
    	lua_pop( L, 1);		// pop parameter_field
    }
    else if ( rate_coeff_model=="OpticallyVariableExponentialDecay" ) {
    	lua_getfield( L, -1, parameter_field.c_str() );
    	rate_coeff = new OpticallyVariableExponentialDecay( L );
    	lua_pop( L, 1);		// pop parameter_field
    }
    else {
    	ostringstream oss;
    	oss << "create_explicit_rate_coeff()" << endl
    	    << "rate_coeff_model: " << rate_coeff_model << " not recognised." << endl;
    	input_error( oss );
    }
    
    return rate_coeff;
}
