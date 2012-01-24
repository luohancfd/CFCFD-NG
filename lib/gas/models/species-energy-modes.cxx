// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../util/source/useful.h"

#include "species-energy-modes.hh"
#include "physical_constants.hh"

// Forward declaration
std::string get_library_species_name( int isp );

using namespace std;

Species_energy_mode::
Species_energy_mode( int isp, double R, double min_massf, string type, int iT )
 : isp_( isp ), R_( R ), min_massf_( min_massf ), type_( type ), iT_( iT ) {}
 
/* ------- One level electronic species energy mode ------- */
// NOTES: - This class uses the Multi_level_electronic expressions, simplified
//          for one level
//        - The primary use is for electrons where theta=0.0, which can actually
//          be reduced further

One_level_electronic::
One_level_electronic( int isp, double R, double min_massf, int g, double theta )
 : Species_energy_mode( isp, R, min_massf, "electronic" ), g_( g ), theta_( theta )
{}
 
double
One_level_electronic::
s_eval_energy( double T )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.31
    return R_ * theta_;
}

double
One_level_electronic::
s_eval_entropy( double T, double p )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.36
    return R_*log( g_*exp( - theta_ / T ) ) + (R_/T)*theta_;
}

double
One_level_electronic::
s_eval_Cv( double T )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.45
    return R_/(T*T) * ( theta_*theta_ - theta_ );
}

/* ------- Two level electronic species energy mode ------- */

Two_level_electronic::
Two_level_electronic( int isp, double R, double min_massf, int g0, 
                      double theta0, int g1, double theta1 )
 : Species_energy_mode( isp, R, min_massf, "electronic" ), g0_( g0 ), 
   g1_( g1 ), theta_( theta1 - theta0 )
{}
 
double
Two_level_electronic::
s_eval_energy( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq. 11.3 p 131
    double tmp = double(g1_)/double(g0_) * exp( - theta_ / T );
    return R_ * theta_ * tmp / ( 1.0 + tmp );
}

double
Two_level_electronic::
s_eval_entropy( double T, double p )
{
    // Ref: Vincenti and Kruger (1975) Ex. 12.3 p 139
    double tmp = double(g1_)/double(g0_) * exp( - theta_ / T );
    double n1 = tmp / ( 1.0 + tmp );
    
    return R_ * ( log(double(g0_)) + log( 1.0 + tmp ) + n1 );
}

double
Two_level_electronic::
s_eval_Cv( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq. 11.4 p 131
    double tmp = double(g1_)/double(g0_) * exp( - theta_ / T );
    return R_ * tmp / pow( T / theta_ * ( 1.0 + tmp ), 2 );
    // return R_ * pow( theta_ / T, 2 ) * tmp / pow( 1.0 + tmp, 2 );
}

/* ------- Multi level electronic species energy mode ------- */

Multi_level_electronic::
Multi_level_electronic( int isp, double R, double min_massf, 
    			vector<int> &g, vector<double> &theta )
 : Species_energy_mode( isp, R, min_massf, "electronic" ), g_( g ), theta_( theta )
{}
 
double
Multi_level_electronic::
s_eval_energy( double T )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.31
    double numerator = 0.0;		// sum of weighted energies
    double denominator = 0.0;		// total partition function
    double tmp;				// level partition function
    for ( size_t i=0; i<theta_.size(); ++i ) {
    	tmp = double(g_[i]) * exp( - theta_[i] / T );
    	numerator += theta_[i] * tmp;
    	denominator += tmp;
    }
    return R_ * numerator / denominator;
}

double
Multi_level_electronic::
s_eval_entropy( double T, double p )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.36
    double tmpA = 0.0;	// sum of weighted energies
    double tmpB = 0.0;	// total partition function
    double tmpC;	// level partition function
    for ( size_t i=0; i<theta_.size(); ++i ) {
    	tmpC = double(g_[i]) * exp( - theta_[i] / T );
    	tmpA += theta_[i] * tmpC;
    	tmpB += tmpC;
    }
    return R_*log( tmpB) + R_/T*log( tmpA/tmpB );
}

double
Multi_level_electronic::
s_eval_Cv( double T )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.45
    double tmpA = 0.0;	// sum of weighted energies
    double tmpB = 0.0;	// sum of weighted energies^2
    double tmpC = 0.0;	// total partition function
    double tmpD;	// level partition function
    for ( size_t i=0; i<theta_.size(); ++i ) {
    	tmpD = double(g_[i]) * exp( - theta_[i] / T );
    	tmpA += theta_[i] * tmpD;
    	tmpB += theta_[i] * theta_[i] * tmpD;
    	tmpC += tmpD;
    }
    return R_/(T*T)*( tmpB / tmpC - tmpA / tmpC );
}

/* ------- Fully excited translation ------- */

Fully_excited_translation::
Fully_excited_translation( int isp, double R, double min_massf )
 : Species_energy_mode( isp, R, min_massf, "translation" ), Cv_( 1.5 * R ), 
   Cp_( 2.5 * R )
{
    // calculate constant expression for entropy calculation
    double m = PC_R_u/R_/PC_Avogadro;
    double tmp = pow( 2.0 * M_PI * m * PC_k_SI / PC_h_SI / PC_h_SI, 1.5 ) * PC_k_SI;
    entropy_constant_ = R_ * ( log( tmp ) + 2.5 );
}

double
Fully_excited_translation::
s_eval_energy( double T )
{
    return Cv_*T;
}

double
Fully_excited_translation::
s_eval_enthalpy( double T )
{
    return Cp_*T;
}

double
Fully_excited_translation::
s_eval_entropy( double T, double p )
{
    // Ref: Vincenti and Kruger (1975) Eq. 9.5b p 123
    // NOTE: - check use of SI units in 'tmp'
    //       - if 'p' is zero then returned value will be '-inf'
    //       - constant part of expression has been pre-calculated in constructor
    // double m = PC_R_u/R_/PC_Avogadro;
    // double tmp = pow( 2.0 * M_PI * m * PC_k_SI / PC_h_SI / PC_h_SI, 1.5 ) * PC_k_SI;
    return Cp_*log(T) - R_*log(p) + entropy_constant_;
}

double
Fully_excited_translation::
s_eval_Cv( double T )
{
    return Cv_;
}

double
Fully_excited_translation::
s_eval_Cp( double T )
{
    return Cp_;
}

/* ------- Fully excited rotation ------- */

Fully_excited_rotation::
Fully_excited_rotation( int isp, double R, double min_massf, double theta, int sigma )
 : Species_energy_mode( isp, R, min_massf, "rotation" ), theta_( theta ), 
   sigma_( sigma )
{
    Cv_ = R_;
}

double
Fully_excited_rotation::
s_eval_energy( double T )
{
    return Cv_*T;
}

double
Fully_excited_rotation::
s_eval_entropy( double T, double p )
{
    // Ref: Vincenti and Kruger (1975) Ex 12.1 p 138
    return R_ * ( log( T / double(sigma_) / theta_ ) + 1.0 );
}

double
Fully_excited_rotation::
s_eval_Cv( double T )
{
    return Cv_;
}

/* ------- Fully excited nonlinear rotation ------- */

Fully_excited_nonlinear_rotation::
Fully_excited_nonlinear_rotation( int isp, double R, double min_massf, 
    				  double theta_A0, double theta_B0, double theta_C0,
    				  int sigma )
 : Species_energy_mode( isp, R, min_massf, "rotation" ), theta_A0_( theta_A0 ), 
   theta_B0_( theta_B0 ), theta_C0_( theta_C0 ), sigma_( sigma )
{
    Cv_ = 1.5 * R_;
}

double
Fully_excited_nonlinear_rotation::
s_eval_energy( double T )
{
    return Cv_*T;
}

double
Fully_excited_nonlinear_rotation::
s_eval_entropy( double T, double p )
{
    // Ref: My derivation from Capitelli ESA STR 246 p21 Q_rot expression
    double Q = sqrt( T*T*T * M_PI / ( theta_A0_ * theta_B0_ * theta_C0_ ) );
    
    return R_ * ( log( Q ) + 1.5 );
}

double
Fully_excited_nonlinear_rotation::
s_eval_Cv( double T )
{
    return Cv_;
}

/* ------- Generic vibration ------- */

Vibration::
Vibration(  int isp, double R, double min_massf, double theta )
: Species_energy_mode( isp, R, min_massf, "vibration" ), theta_( theta ) {}

double
Vibration::
s_eval_energy( double T )
{
    cout << "Vibration::s_eval_energy()" << endl
         << "This function is not meant for use." << endl
         << "Bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}

double
Vibration::
s_eval_entropy( double T, double p )
{
    cout << "Vibration::s_eval_entropy()" << endl
         << "This function is not meant for use." << endl
         << "Bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}

double
Vibration::
s_eval_Cv( double T )
{
    cout << "Vibration::s_eval_Cv()" << endl
         << "This function is not meant for use." << endl
         << "Bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}

/* ------- Un-truncated harmonic vibration ------- */

Harmonic_vibration::
Harmonic_vibration(  int isp, double R, double min_massf, double theta )
: Vibration( isp, R, min_massf, theta ) {}

double
Harmonic_vibration::
s_eval_energy( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq. 12.12 p 135
    return R_ * theta_ / ( exp( theta_ / T ) - 1.0 );
}

double
Harmonic_vibration::
s_eval_entropy( double T, double p )
{
    // Ref: Vincenti and Kruger (1975) Ex 12.2 p 138
    return R_ * ( - log ( 1 - exp( - theta_ / T ) ) + ( theta_ / T ) / ( exp( theta_ / T ) - 1.0 ) );
}

double
Harmonic_vibration::
s_eval_Cv( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq 12.13 p 135
    double theta_2T = theta_ / ( 2.0 * T );
    return R_ * pow( theta_2T / sinh( theta_2T ), 2 );
}

double
Harmonic_vibration::
s_eval_Q( double T, double A )
{
    // Ref: Vincenti and Kruger (1975) p 135
    return R_ * 1.0 / ( exp( theta_ / T ) - 1.0 );
}

/* ------- Truncated harmonic vibration ------- */

Truncated_harmonic_vibration::
Truncated_harmonic_vibration(  int isp, double R, double min_massf, double theta, double theta_D )
: Vibration( isp, R, min_massf, theta ), theta_D_( theta_D ) {}

double
Truncated_harmonic_vibration::
s_eval_energy( double T )
{
    // Ref: Gollan, R.G. PhD Thesis June 2008
    double exp_theta_D_T = exp( theta_D_ / T );
    if ( isnan(exp_theta_D_T) ) {
    	// use harmonic oscillator solution
    	return R_ * ( theta_ / ( exp( theta_ / T ) - 1.0 ) );
    }
    else {
    	return R_ * ( theta_ / ( exp( theta_ / T ) - 1.0 ) - theta_D_ / ( exp_theta_D_T - 1.0 ) );
    }
}

double
Truncated_harmonic_vibration::
s_eval_entropy( double T, double p )
{
    // Ref: My own derivation from Vincenti and Kruger's s(Q,T) entropy expression
    //      See myproject/thesis/LaTeX/thermal_nonequilibrium/version2/tneq.tex
    double exp_theta_v_T = exp( - theta_ / T );
    double exp_theta_D_T = exp( - theta_D_ / T );
    
    return R_ * ( - log ( ( 1.0 - exp_theta_v_T ) / ( 1.0 - exp_theta_D_T ) ) + \
    		 ( 1.0 / T ) * ( ( theta_ * exp_theta_v_T ) / ( 1.0 - exp_theta_v_T ) - \
    		     		   theta_D_ * exp_theta_D_T / ( 1.0 - exp_theta_D_T ) ) );
}

double
Truncated_harmonic_vibration::
s_eval_Cv( double T )
{
    // Ref: lib/gas_models2/source/molecule.cxx::Diatomic_molecule_THO
    //      (my own derivation)
    double tmpA = theta_ / T;
    double tmpB = theta_D_ / T;
    double exp_tmpA = exp( tmpA );
    double exp_tmpB = exp( tmpB );
    double tmpC = tmpA * tmpA * exp_tmpA /  ( (exp_tmpA - 1.0) * (exp_tmpA - 1.0) );
    double tmpD = tmpB * tmpB * exp_tmpB /  ( (exp_tmpB - 1.0) * (exp_tmpB - 1.0) );
    
    // if tmpD is nan then it is converging to zero
    if ( isnan(tmpD) ) tmpD = 0.0;
    
    return R_ * ( tmpC - tmpD );
}

double
Truncated_harmonic_vibration::
s_eval_Q( double T, double A )
{
    double theta_lim = A;
    if ( theta_lim < 0.0 ) theta_lim = theta_D_;
    // Ref: DFP thesis
    return ( 1.0 - exp( - theta_ / T ) ) / ( 1.0 - exp( - theta_lim / T ) );
}

double
Truncated_harmonic_vibration::
s_eval_HO_energy( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq. 12.12 p 135
    return R_ * theta_ / ( exp( theta_ / T ) - 1.0 );
}

/* ------- Fully coupled diatomic internal mode ------- */

Fully_coupled_diatom_internal::
Fully_coupled_diatom_internal( int isp, double R, double min_massf, int sigma_r,
    std::vector< std::vector<double> > &elev_data )
 : Species_energy_mode( isp, R, min_massf, "internal" ), sigma_r_( sigma_r ),
 m_( PC_R_u / R_ / PC_Avogadro )
{
    // Create the electronic levels to do the work
    for ( size_t ilev=0; ilev < elev_data.size(); ++ilev ) {
    	cout << "Fully_coupled_diatom_internal -> attempting to create ilev = " << ilev << endl;
    	elevs_.push_back( new Diatom_electronic_level( elev_data[ilev] ) );
    }
}

Fully_coupled_diatom_internal::
~Fully_coupled_diatom_internal()
{
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev )
    	delete elevs_[ilev];
}

double
Fully_coupled_diatom_internal::
s_eval_energy( double T )
{
    // Loop over all rovibronic levels, sum up contributions from definition
    // NOTE: referencing all energies from the respective ground states
    double E_weighted = 0.0;
    double Q_total = 0.0;
    double E_el_0 = elevs_[0]->E;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
    	double E_el = elevs_[ilev]->E - E_el_0;
    	double E_vib_0 = elevs_[ilev]->eval_E_vib( 0 );
    	for ( int iV=0; iV<elevs_[ilev]->V_max; ++iV ) {
    	    double E_vib = elevs_[ilev]->eval_E_vib( iV ) - E_vib_0;
    	    double E_rot_0 = elevs_[ilev]->eval_E_rot(iV,0);
    	    for ( int iJ=0; iJ<elevs_[ilev]->J_max[iV]; ++iJ ){
    	    	double E_rot = elevs_[ilev]->eval_E_rot(iV,iJ) - E_rot_0;
    	    	double E_rot_dash = E_el + E_vib + E_rot;
    	    	double Q_rot_dash = elevs_[ilev]->g * double(2*iJ+1) * exp( - E_rot_dash / PC_k_SI / T );
    	    	if ( isinf( Q_rot_dash ) || isnan( Q_rot_dash ) ) {
    	    	    cout << get_library_species_name(isp_) << ": ilev = " << ilev << ", iV = " << iV << ", iJ = " << iJ << ", E_rot_dash = " << E_rot_dash << ", Q_rot_dash = " << Q_rot_dash << endl;
    	    	    exit(0);
    	    	}
    	    	E_weighted += E_rot_dash * Q_rot_dash;
    	    	Q_total += Q_rot_dash;
    	    }
    	}
    }
    
    return E_weighted / Q_total / m_;	// convert J/particle -> J/kg
}

double
Fully_coupled_diatom_internal::
s_eval_entropy( double T, double p )
{
    // Loop over all rovibronic levels, sum up contributions from definition:
    // S = N k ln(Q) + E/T; s = R ln(Q) + e/T
    // NOTE: referencing all energies from the respective ground states
    double E_weighted = 0.0;
    double Q_total = 0.0;
    double E_el_0 = elevs_[0]->E;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
    	double E_el = elevs_[ilev]->E - E_el_0;
    	double E_vib_0 = elevs_[ilev]->eval_E_vib( 0 );
    	for ( int iV=0; iV<elevs_[ilev]->V_max; ++iV ) {
    	    double E_vib = elevs_[ilev]->eval_E_vib( iV ) - E_vib_0;
    	    double E_rot_0 = elevs_[ilev]->eval_E_rot(iV,0);
    	    for ( int iJ=0; iJ<elevs_[ilev]->J_max[iV]; ++iJ ){
    	    	double E_rot = elevs_[ilev]->eval_E_rot(iV,iJ) - E_rot_0;
    	    	double E_rot_dash = E_el + E_vib + E_rot;
    	    	double Q_rot_dash = elevs_[ilev]->g * double(2*iJ+1) * exp( - E_rot_dash / PC_k_SI / T );
    	    	if ( isinf( Q_rot_dash ) || isnan( Q_rot_dash ) ) {
    	    	    cout << get_library_species_name(isp_) << ": ilev = " << ilev << ", iV = " << iV << ", iJ = " << iJ << ", E_rot_dash = " << E_rot_dash << ", Q_rot_dash = " << Q_rot_dash << endl;
    	    	    exit(0);
    	    	}
    	    	E_weighted += E_rot_dash * Q_rot_dash;
    	    	Q_total += Q_rot_dash;
    	    }
    	}
    }
    
    double e = E_weighted / Q_total / m_;	// convert J/particle -> J/kg
    
    return R_*log(Q_total) + e/T;		// units of J/kg/K
}

double
Fully_coupled_diatom_internal::
s_eval_Cv( double T )
{
    // differentiation of energy wrt T using the quotient rule
    // NOTE: referencing all energies from the respective ground states
    double u = 0.0, u_dash = 0.0;
    double v = 0.0, v_dash = 0.0;
    double E_el_0 = elevs_[0]->E;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
    	double E_el = elevs_[ilev]->E - E_el_0;
    	double E_vib_0 = elevs_[ilev]->eval_E_vib( 0 );
    	for ( int iV=0; iV<elevs_[ilev]->V_max; ++iV ) {
    	    double E_vib = elevs_[ilev]->eval_E_vib( iV ) - E_vib_0;
    	    double E_rot_0 = elevs_[ilev]->eval_E_rot(iV,0);
    	    for ( int iJ=0; iJ<elevs_[ilev]->J_max[iV]; ++iJ ){
    	    	double E_rot = elevs_[ilev]->eval_E_rot(iV,iJ) - E_rot_0;
    	    	double E_rot_dash = E_el + E_vib + E_rot;
    	    	double g = elevs_[ilev]->g * double(2*iJ+1);
    	    	double Q_rot_dash = g * exp( - E_rot_dash / PC_k_SI / T );
    	    	u += E_rot_dash * Q_rot_dash;
    	    	u_dash += E_rot_dash / PC_k_SI / T / T * E_rot_dash * Q_rot_dash;
    	    	v += Q_rot_dash;
    	    	v_dash += E_rot_dash / PC_k_SI / T / T * Q_rot_dash;
    	    }
    	}
    }
    
    return ( u_dash * v - u * v_dash ) / ( v * v ) / m_;	// convert J/particle/K -> J/kg/K
}

