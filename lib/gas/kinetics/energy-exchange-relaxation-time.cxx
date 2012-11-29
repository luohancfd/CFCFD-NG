// Author: Daniel F. Potter
// Date: 18-Nov-2009

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "../models/gas_data.hh"
#include "../models/physical_constants.hh"
#include "../models/chemical-species.hh"
#include "../models/chemical-species-library.hh"
#include "../models/species-energy-modes.hh"

#include "energy-exchange-relaxation-time.hh"

#define THIVET_MW_FORMULATION          1
#define N_TOTAL_FOR_PARK_VT_CORRECTION 1

using namespace std;

Relaxation_time::
Relaxation_time() {}

Relaxation_time::
~Relaxation_time() {}

VT_MillikanWhite::
VT_MillikanWhite( lua_State * L )
    : Relaxation_time()
{
    // 1. Vibrating species data
    string p_name = get_string(L, -1, "p_name");
    Chemical_species * p = get_library_species_pointer_from_name( p_name );
    M_p_ = p->get_M();
    ip_ = p->get_isp();
    Vibration * p_vib =  dynamic_cast<Vibration*>(p->get_mode_pointer_from_type("vibration"));
    theta_ = p_vib->get_theta();
    iT_ = p->get_iT_trans();
    
    // 2. Colliding species data
    lua_getfield(L, -1, "q_names" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite::VT_MillikanWhite():\n";
	ost << "Error in the declaration of colliding species : a table 'q_names' is expected.\n";
	input_error(ost);
    }
    int nqs = lua_objlen(L, -1);
    for ( int i=0; i<nqs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	const char* q_name = luaL_checkstring(L, -1);
    	lua_pop(L, 1);
    	Chemical_species * X = get_library_species_pointer_from_name( q_name );
    	if ( p_name==q_name )
    	    homogenous_flags_.push_back( true );
    	else
    	    homogenous_flags_.push_back( false );
    	M_qs_.push_back( X->get_M() );
    	iqs_.push_back( X->get_isp() );
    	mus_.push_back( ((M_p_ * M_qs_.back()) / (M_p_ + M_qs_.back())) * 1000.0 );
    	a_.push_back( 1.16e-3*sqrt(mus_[i])*pow(theta_, 4.0/3.0) );
	b_.push_back( 0.015*pow(mus_[i], 1.0/4.0) );
    }
    lua_pop(L,1);	// pop q_names
    
    // 3. Check for manually set a values
    lua_getfield(L, -1, "a_values" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite::VT_MillikanWhite():\n";
	ost << "Error in the declaration of the 'a' parameter values : a table 'a_values' is expected.\n";
	input_error(ost);
    }
    int nas = lua_objlen(L, -1);
    for ( int i=0; i< nas; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	double a = luaL_checknumber(L,-1);
    	lua_pop(L, 1);
    	if ( a==-1.0 ) continue;
    	else {
    	    cout << "p = " << p_name << ", i = " << i << ", MW a = " << a_[i] << ", a given = " << a << endl;
    	    a_[i] = a;
    	}
    }
    lua_pop(L,1);	// pop a_values
    
    // 4. Check for manually set b values
    lua_getfield(L, -1, "b_values" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite::VT_MillikanWhite():\n";
	ost << "Error in the declaration of the 'b' parameter values : a table 'b_values' is expected.\n";
	input_error(ost);
    }
    int nbs = lua_objlen(L, -1);
    for ( int i=0; i<nbs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	double b = luaL_checknumber(L,-1);
    	lua_pop(L, 1);
    	if ( b==-1.0 ) continue;
    	else {
    	    cout << "p = " << p_name << ", i = " << i << ", MW b = " << b_[i] << ", b given = " << b << endl;
    	    b_[i] = b;
    	}
    }
    lua_pop(L,1);	// pop b_values
}

VT_MillikanWhite::
~VT_MillikanWhite() {}

double
VT_MillikanWhite::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double p_atm;
    double tau_q;
    double tau_inv = 0.0;
    
    // NOTE: assuming thermally perfect gases
    for ( size_t i=0; i<iqs_.size(); ++i ) {
#       if THIVET_MW_FORMULATION==1
        // Use the mixture partial pressure
	if ( homogenous_flags_[i] ) {
	    p_atm = Q.rho * Q.massf[ip_] * PC_R_u / M_p_ * Q.T[iT_] / PC_P_atm;
	}
	else {
	    p_atm = Q.rho * ( Q.massf[ip_] * PC_R_u / M_p_ +  Q.massf[iqs_[i]] * PC_R_u / M_qs_[i] ) * Q.T[iT_] / PC_P_atm;
	}
#       else
        // Use the collider partial pressure
	p_atm = Q.rho * Q.massf[iqs_[i]] * PC_R_u / M_qs_[i] * Q.T[iT_] / PC_P_atm;
#       endif
	
	tau_q = (1.0/p_atm)*exp(a_[i]*(pow(Q.T[iT_], -1.0/3.0) - b_[i]) - 18.42);

#       if THIVET_MW_FORMULATION==1
	// Scale by collider mole-fraction
	tau_inv += molef[iqs_[i]] / tau_q;
#       else
	// Omit mole-fraction scaling
	tau_inv += 1.0 / tau_q;
#       endif
    }

    return 1.0/tau_inv;
}

VT_high_temperature_cross_section *
create_VT_high_temperature_cross_section_model( lua_State *L )
{
    string type = get_string(L, -1, "type");
    if ( type=="Park" ) return new Park_VT_HTCS(L);
    else if ( type=="Fujita" ) return new Fujita_VT_HTCS();
    else {
	cout << "create_VT_high_temperature_correction_model():\n";
	cout << "HTC model type: " << type << " was not recognised.\n";
	cout << "Bailing out!\n";
	exit( BAD_INPUT_ERROR );
    }
}

Park_VT_HTCS::Park_VT_HTCS( lua_State *L )
{
    sigma_dash_ = get_positive_number(L, -1, "sigma_dash" );
}

Park_VT_HTCS::~Park_VT_HTCS() {}

double
Park_VT_HTCS::s_eval_sigma( double T )
{
    // if ( T > 5.0e4 ) T = 5.0e4;
    return sigma_dash_ * ( 50000.0 / T ) * ( 50000.0 / T );		// cm**2
}

Fujita_VT_HTCS::Fujita_VT_HTCS() {}

Fujita_VT_HTCS::~Fujita_VT_HTCS() {}

double
Fujita_VT_HTCS::s_eval_sigma( double T )
{
    return ( 1.8e-14 ) / sqrt( T ) + ( 1.0e-4 ) / pow( T, 3.0 );	// cm**2
}

VT_MillikanWhite_HTC::
VT_MillikanWhite_HTC( lua_State * L )
    : Relaxation_time()
{
    // 1. Vibrating species data
    string p_name = get_string(L, -1, "p_name");
    Chemical_species * p = get_library_species_pointer_from_name( p_name );
    M_p_ = p->get_M();
    ip_ = p->get_isp();
    Vibration * p_vib =  dynamic_cast<Vibration*>(p->get_mode_pointer_from_type("vibration"));
    theta_ = p_vib->get_theta();
    iT_ = p->get_iT_trans();
    
    // 2. Colliding species data
    lua_getfield(L, -1, "q_names" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite::VT_MillikanWhite():\n";
	ost << "Error in the declaration of colliding species : a table 'q_names' is expected.\n";
	input_error(ost);
    }
    int nqs = lua_objlen(L, -1);
    for ( int i=0; i<nqs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	const char* q_name = luaL_checkstring(L, -1);
    	lua_pop(L, 1);
    	Chemical_species * X = get_library_species_pointer_from_name( q_name );
    	if ( p_name==q_name )
    	    homogenous_flags_.push_back( true );
    	else
    	    homogenous_flags_.push_back( false );
    	M_qs_.push_back( X->get_M() );
    	iqs_.push_back( X->get_isp() );
    	mus_.push_back( ((M_p_ * M_qs_.back()) / (M_p_ + M_qs_.back())) * 1000.0 );
    	a_.push_back( 1.16e-3*sqrt(mus_[i])*pow(theta_, 4.0/3.0) );
	b_.push_back( 0.015*pow(mus_[i], 1.0/4.0) );
    }
    
    lua_pop(L,1);	// pop q_names
    
    // 3. Check for manually set a values
    lua_getfield(L, -1, "a_values" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite::VT_MillikanWhite():\n";
	ost << "Error in the declaration of the 'a' parameter values : a table 'a_values' is expected.\n";
	input_error(ost);
    }
    int nas = lua_objlen(L, -1);
    for ( int i=0; i< nas; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	double a = luaL_checknumber(L,-1);
    	lua_pop(L, 1);
    	if ( a==-1.0 ) continue;
    	else {
    	    cout << "p = " << p_name << ", i = " << i << ", MW a = " << a_[i] << ", a given = " << a << endl;
    	    a_[i] = a;
    	}
    }
    lua_pop(L,1);	// pop a_values
    
    // 4. Check for manually set b values
    lua_getfield(L, -1, "b_values" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite::VT_MillikanWhite():\n";
	ost << "Error in the declaration of the 'b' parameter values : a table 'b_values' is expected.\n";
	input_error(ost);
    }
    int nbs = lua_objlen(L, -1);
    for ( int i=0; i<nbs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	double b = luaL_checknumber(L,-1);
    	lua_pop(L, 1);
    	if ( b==-1.0 ) continue;
    	else {
    	    cout << "p = " << p_name << ", i = " << i << ", MW b = " << b_[i] << ", b given = " << b << endl;
    	    b_[i] = b;
    	}
    }
    lua_pop(L,1);	// pop b_values
    
    // 5. High-temperature correction cross-section model
    lua_getfield(L,-1,"HTCS_model");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite_HTC::VT_MillikanWhite_HTCS():\n";
	ost << "Error in the declaration of HTCS model: a table is expected.\n";
	input_error(ost);
    }
    HTC_model_ = create_VT_high_temperature_cross_section_model( L );
    lua_pop(L,1);	// pop HTCS_model
}

VT_MillikanWhite_HTC::
~VT_MillikanWhite_HTC()
{
    delete HTC_model_;
}

double
VT_MillikanWhite_HTC::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double p_atm;
    double m_av;
    double tau_q;
    double tau_inv = 0.0;
    double sigma;
#   if N_TOTAL_FOR_PARK_VT_CORRECTION==1
    double n_col = Q.p / PC_k_SI / Q.T[0] * 1.0e-6;
#   else
    double n_col;
#   endif

    
    // NOTE: assuming thermally perfect gases for the p_atm calc
    for ( size_t i=0; i<iqs_.size(); ++i ) {
#       if THIVET_MW_FORMULATION==1
        // Use binary mixture partial pressure
	if ( homogenous_flags_[i] ) {
	    p_atm = Q.rho * Q.massf[ip_] * PC_R_u / M_p_ * Q.T[iT_] / PC_P_atm;		// atmospheres
	    m_av = M_p_ / PC_Avogadro * 1.0e3;						// grams
	}
	else {
	    p_atm = Q.rho * ( Q.massf[ip_] * PC_R_u / M_p_ +  Q.massf[iqs_[i]] * PC_R_u / M_qs_[i] ) * Q.T[iT_] / PC_P_atm;
	    m_av = 0.5 * ( M_p_ + M_qs_[i] ) / PC_Avogadro * 1.0e3;
	}
#       if N_TOTAL_FOR_PARK_VT_CORRECTION==0
	n_col = p_atm * PC_P_atm / PC_k_SI / Q.T[0] * 1.0e-6;                           // colliding particles / cm**3
#       endif
#       else
	// Use collider partial pressure
	p_atm = Q.rho * Q.massf[iqs_[i]] * PC_R_u / M_qs_[i] * Q.T[iT_] / PC_P_atm;
        if ( homogenous_flags_[i] ) {
#           if N_TOTAL_FOR_PARK_VT_CORRECTION==0
            n_col = Q.rho * Q.massf[ip_] / M_p_ * PC_Avogadro * 1.0e-6;                 // colliding particles / cm**3
#           endif
            m_av = M_p_ / PC_Avogadro * 1.0e3;                                          // grams
        }
        else {
#           if N_TOTAL_FOR_PARK_VT_CORRECTION==0
            n_col = Q.rho * ( Q.massf[ip_] * PC_R_u / M_p_ +  Q.massf[iqs_[i]] * PC_R_u / M_qs_[i] ) * PC_Avogadro * 1.0e-6;
#           endif
            m_av = 0.5 * ( M_p_ + M_qs_[i] ) / PC_Avogadro * 1.0e3;
        }
#       endif
	
	tau_q = (1.0/p_atm)*exp(a_[i]*(pow(Q.T[iT_], -1.0/3.0) - b_[i]) - 18.42);
	// high-temperature correction
	sigma = HTC_model_->eval_sigma(Q.T[iT_]);				// cm**-2
	tau_q += 1.0 / ( n_col * sqrt( 8.0 * PC_k_CGS * Q.T[iT_] / ( M_PI * m_av ) ) * sigma );
#       if THIVET_MW_FORMULATION==1
	// mole-fraction scaling (see Thivet paper)
	tau_inv += molef[iqs_[i]] / tau_q;
#       else
	// Omit mole-fraction scaling
	tau_inv += 1.0 / tau_q;
#       endif
    }

    return 1.0/tau_inv;
}

ET_AppletonBray::
ET_AppletonBray( lua_State * L )
    : Relaxation_time()
{
    // Initialise ions
    lua_getfield(L,-1,"ions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ET_AppletonBray::ET_AppletonBray():\n";
	ost << "Error in declaration of ions: a table is expected.\n";
	input_error(ost);
    }
    int n_ions = lua_objlen(L, -1);
    for ( int i=0; i<n_ions; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	tau_ETs_.push_back( new ET_Ion(L) );
    	lua_pop(L, 1);
    }
    lua_pop(L, 1);	// pop ions
    
    // Initialise neutrals
    lua_getfield(L,-1,"neutrals");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ET_AppletonBray::ET_AppletonBray():\n";
	ost << "Error in declaration of neutrals: a table is expected.\n";
	input_error(ost);
    }
    int n_neutrals = lua_objlen(L, -1);
    for ( int i=0; i<n_neutrals; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	tau_ETs_.push_back( new ET_Neutral(L) );
    	lua_pop(L, 1);
    }
    lua_pop(L, 1);	// pop neutrals
}

double
ET_AppletonBray::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double tau_inv = 0.0;
    for ( size_t i=0; i<tau_ETs_.size(); ++i ) {
    	double tau = tau_ETs_[i]->compute_relaxation_time(Q,molef);
    	if ( tau < 0.0 ) continue;
    	tau_inv += 1.0/tau;
    }
    
    return 1.0/tau_inv;
}

ET_AppletonBray::
~ET_AppletonBray()
{
    for ( size_t i=0; i<tau_ETs_.size(); ++i )
    	delete tau_ETs_[i];
}

ET_Ion::
ET_Ion( lua_State * L )
    : Relaxation_time()
{
    // 1. Colliding species data
    string c_name = get_string(L, -1, "c_name");
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    if ( X->get_Z()<1 ) {
	ostringstream ost;
	ost << "ET_Ion::ET_Ion():\n";
	ost << "Error in the declaration of colliding species: " << c_name << " is not an ion.\n";
	input_error(ost);
    }    
    ic_ = X->get_isp();
    M_c_ = X->get_M();
    
    // 2. Electron data
    X = get_library_species_pointer_from_name( "e_minus" );
    ie_ = X->get_isp();
    iTe_ = X->get_iT_trans();
    M_e_ = X->get_M();
}

ET_Ion::
~ET_Ion() {}

double
ET_Ion::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;	// Number density of colliders (in cm**-3)
    double n_e = Q.massf[ie_] * Q.rho / M_e_ * PC_Avogadro * 1.0e-6;	// Number density of electrons (in cm**-3)
    
    // If there are no participating species set a negative relaxation time
    if ( n_e==0.0 || n_c==0.0 ) return -1.0;

    double tmpA = 8.0 / 3.0 * sqrt( M_PI / PC_m_CGS );
    double tmpB = n_c * pow ( PC_e_CGS, 4.0 ) / pow ( 2.0 * PC_k_CGS * Q.T[iTe_], 1.5 );
    double tmpC = log( pow( PC_k_CGS * Q.T[iTe_], 3.0 ) / ( M_PI * n_e * pow ( PC_e_CGS, 6.0 ) ) );
    double nu_ec = tmpA * tmpB * tmpC;
        
    double tau_ec = M_c_ / nu_ec;
    
    return tau_ec;
    
}

ET_Neutral::
ET_Neutral( lua_State * L )
    : Relaxation_time()
{
    int nCs;

    // 1. Colliding species data
    string c_name = get_string(L, -1, "c_name");
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    if ( X->get_Z()!=0 ) {
	ostringstream ost;
	ost << "ET_Neutral::ET_Neutral():\n";
	ost << "Error in the declaration of colliding species: " << c_name << " is not a neutral.\n";
	input_error(ost);
    }    
    ic_ = X->get_isp();
    M_c_ = X->get_M();
    
    // 2. Electron data
    X = get_library_species_pointer_from_name( "e_minus" );
    iTe_ = X->get_iT_trans();
    
    // 3.  Sigma quadratic curve fit coefficients
    // 3a. Temperature switch
    T_switch_ = get_positive_number(L, -1, "T_switch");
    // 3b. Low temperature range
    lua_getfield(L, -1, "LT_sigma_coefficients" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ET_Neutral::ET_Neutral():\n";
	ost << "Error in the declaration of low temperature sigma coefficients: a table is expected.\n";
	input_error(ost);
    }
    nCs = lua_objlen(L, -1);
    if ( nCs!=3 ) {
	ostringstream ost;
	ost << "ET_Neutral::ET_Neutral():\n";
	ost << "Error in the declaration of low temperature sigma coefficients: 3 coefficients are expected.\n";
	input_error(ost);
    }
    for ( int i=0; i<nCs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	LT_C_.push_back( luaL_checknumber(L,-1) );
    	lua_pop(L,1);
    }
    lua_pop(L,1);	// pop 'LT_sigma_coefficients'
    // 3c. High temperature range
    lua_getfield(L, -1, "HT_sigma_coefficients" );
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "ET_Neutral::ET_Neutral():\n";
        ost << "Error in the declaration of high temperature sigma coefficients: a table is expected.\n";
        input_error(ost);
    }
    nCs = lua_objlen(L, -1);
    if ( nCs!=3 ) {
        ostringstream ost;
        ost << "ET_Neutral::ET_Neutral():\n";
        ost << "Error in the declaration of high temperature sigma coefficients: 3 coefficients are expected.\n";
        input_error(ost);
    }
    for ( int i=0; i<nCs; ++i ) {
        lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
        HT_C_.push_back( luaL_checknumber(L,-1) );
        lua_pop(L,1);
    }
    lua_pop(L,1);       // pop 'HT_sigma_coefficients'
}

ET_Neutral::
~ET_Neutral()
{
    LT_C_.resize(0);
    HT_C_.resize(0);
}

double
ET_Neutral::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;	// Number density of colliders (in cm**-3)
    
    // If there are no colliding species set a negative relaxation time
    if ( n_c==0.0 ) return -1.0;

    // calculate sigma using the appropriate curve fit for the electron temperature
    double T = Q.T[iTe_];
    double sigma_ec = 0.0;
    if ( T < T_switch_ )
        sigma_ec = LT_C_[0] + LT_C_[1] * T + LT_C_[2] * T * T;
    else
        sigma_ec = HT_C_[0] + HT_C_[1] * T + HT_C_[2] * T * T;
    /* Convert sigma_ec from m**2 -> cm**2 */
    sigma_ec *= 1.0e4;
    double nu_ec = n_c * sigma_ec * sqrt( 8.0 * PC_k_CGS * T / ( M_PI * PC_m_CGS ) );
    double tau_ec = M_c_ / nu_ec;
    
    return tau_ec;
}

/**
 * \brief Calculate the collider distance for like colliders.
 *
 * This function (A) is applicable for like colliders, ie.
 * polar-polar or nonpolar-nonpolar.
 * See Hirschfelder. p. 222
 *
 * \[ r^0_{pq} = \frac{1}{2} \left( r^0_{pp} + r^{0}_{qq} \right) \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double collider_distance_A(Diatomic_species &p, Diatomic_species &q)
{
    return (0.5 * (p.get_r0() + q.get_r0()));
}

/**
 * \brief Calculate the collider distance for unlike colliders.
 *
 * This function (B) is applicable when the colliders are unalike,
 * ie. polar-nonpolar
 *
 * \[ r^{0}_{pn} = \frac{1}{2} \left( r^{0}_{pp} + r^{0}_{nn} \right) \xi^{-1/6} \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double collider_distance_B(Diatomic_species &p, Diatomic_species &n)
{
    double r_pn, xi;

    xi = xi_correction(p, n);
    r_pn = 0.5 * (p.get_r0() + n.get_r0()) * pow(xi, -1.0/6.0);

    return (r_pn);
}

/**
 * \brief Calculate xi (a correction factor)
 *
 * This function calculates a correction factor needed for the 
 * the empirical combining laws when there is a collision
 * between polar and nonpolar molecules.
 *
 * \[ \xi = \left\[ 1 + \frac{\alpha_n \mu^{*2}_p}{4 \r^3_n}
 *                      \sqrt{\frac{\epsilon_p}{\epsilon_n}} \right\] \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double xi_correction(Diatomic_species &p, Diatomic_species &n)
{
    double mu_star, tmp_a, tmp_b, xi;

    mu_star = p.get_mu_B() / sqrt( p.get_eps0() * pow(p.get_r0(), 3.0));
    tmp_a = (n.get_alpha() * mu_star) / ( 4.0 * pow(n.get_r0(), 3.0));
    tmp_b = sqrt(p.get_eps0() / n.get_eps0());

    xi = 1.0 + tmp_a * tmp_b;

    return (xi);
}

/**
 * \brief Calculate eps0_pq: potential well for like colliders.
 *
 * This function (A) calculates the potential well in a Lennard-Jones
 * type collision for two alike molecules, ie. polar-polar or
 * nonpolar-nonpolar.
 *
 * \[ \epsilon^0_{pq} = \sqrt{ \epsilon^0_{pp} \epsilon^0_{qq} } \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double potential_well_A(Diatomic_species &p, Diatomic_species &q)
{
    return ( sqrt(p.get_eps0() * q.get_eps0()) );
}

/**
 * \brief Calculate eps0_pq: potential well for unlike colliders.
 *
 * This function (B) calculates the potential well for collisions
 * between unalike molecules, ie. polar-nonpolar.
 *
 * \[ \epsilon^0_{pn} = \sqrt{\epsilon^0_{pp} \epsilon^0_{qq}} \xi^2 \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double potential_well_B(Diatomic_species &p, Diatomic_species &n)
{
    double eps_pq, xi;
    
    xi = xi_correction(p, n);
    eps_pq = sqrt(p.get_eps0() * n.get_eps0()) * pow(xi, 2.0);

    return (eps_pq);
}

/**
 * \brief Calculate the collision frequency based on the hard-sphere model.
 *
 * \[ Z^{pq}_{\mbox{coll}} = 2\sqrt{2} n \sigma^2_{pq} \sqrt{\pi k_b T / \mu_{pq}} \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double collision_frequency(double sigma, double mu, double T, double nd)
{
    double tmp;
    double Z;

    tmp = sqrt(2.0 * M_PI * PC_k_SI * T / mu);
    Z = 2.0 * nd * sigma * sigma * tmp;

    //    cout << "Z= " << Z << " nd= " << nd << " sigma= " << sigma << " mu= " << mu << endl;

    return (Z);
}

/**
 * \brief Calculate expansion parameter beta.
 *
 * \[ \beta = \{ \frac{1}{2} (2 \epsilon^0_{pq} / \mu_{pq})^9 
 *              (3h\mu_{pq}/\pi^2 r^0_{pq} k T \Delta E)^6 \}^{1/19} \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_beta(double eps0, double mu, double r0,
		double del_E, double T)
{

    double tmp_a, tmp_b;
    double beta;

    tmp_a = 2.0 * eps0 / mu;
    tmp_b = 3.0 * PC_h_SI * mu /
	(M_PI * M_PI * r0 * PC_k_SI * T * del_E);

    beta = pow( 0.5 * pow(tmp_a, 9) * pow(tmp_b, 6), 1.0/19.0);

    return (beta);
}

/**
 * \brief Calculate D_star using a first-order expression.
 *
 * Implements:
 * \[ D^* = (1/\beta^2)( 1 - \frac{14}{19} \beta ) \]
 *
 * \author Rowan J Gollan
 * \date 01-March-2005
 **/

double SSH_D_star(double beta)
{
    return ((1.0 / (beta * beta)) * (1.0 - (14.0/19.0) * beta));
}

/**
 * \brief Calculate r_c_star.
 *
 * r_c_star is the classical distance of closest
 * approach based on kinetic theory.
 * Implements:
 * \[ r^{c*}_{pq} = r^0_{pq}/ ( (1/(2\beta)^{1/6})
 *                              (1 + 2/19 \beta) ) \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_r_c_star(double beta, double r0)
{
    double tmp_a, tmp_b;

    tmp_a = 1.0 / pow(2.0 * beta, 1.0/6.0);
    tmp_b = 1.0 + 2.0/19.0 * beta;

    return ( r0 / (tmp_a * tmp_b) );
}

double SSH_r_c_star2(double D_star, double r0)
{

    return (r0 * pow( (0.5 + 0.5*sqrt(D_star)), 6.0));

}

/**
 * \brief Calculate del_star.
 *
 * del_star is a parameter of the exponential potential.
 * Implements:
 * \[ \delta^*_{pq} = \frac{1}{r^0_{pq}} ( 12/(2\beta)^{1/6}) (1 + 2/19 \beta) \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_del_star(double beta, double r0)
{
    double tmp_a, tmp_b;

    tmp_a = 12.0 / pow(2.0 * beta, 1.0/6.0);
    tmp_b = 1.0 + 21.0 / 19.0 * beta;

    return ( (tmp_a * tmp_b) / r0 );
}

double SSH_del_star2(double D_star, double r0)
{
    double tmp_a, tmp_b;

    tmp_a = sqrt(D_star);
    tmp_b = pow((0.5 + 0.5*tmp_a), 1.0/6.0);
    return ( (12.0/r0) * tmp_b * (1.0 + 1.0/tmp_a) );
}

/**
 * \brief Calculate alpha_pq.
 *
 * alpha_pq is the main parameter which describes
 * the relaxation tims'e dependence on temperature.
 * Implements:
 * \[ \alpha_{pq} = \frac{16 \pi^2 \mu_{pq} \Delta E^2}
 *                       {\delta^{*2}_{pq} h^2 k}  \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_alpha_pq(double mu, double del_E, double del_star)
{
    double tmp_a, tmp_b;

    tmp_a = 16.0 * pow(M_PI, 4) * mu * del_E * del_E;
    tmp_b = del_star * del_star * PC_h_SI * PC_h_SI * PC_k_SI;

    return ( tmp_a/tmp_b );
}

/**
 * \brief Calculate chi_pq.
 *
 * chi_pq is an interaction parameter of molecules
 * p and q.
 * Implements:
 * \[ \chi^*_{pq} = \frac{1}{2} \left( \frac{alpha_{pq}}{T} \right)^{1/3} \]
 * 
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_chi_pq(double alpha_pq, double T)
{
    return ( 0.5 * pow(alpha_pq/T, 1.0/3.0) );
}

/**
 * \brief Calculate Z_0.
 *
 * Z_0 is a steric factor reqpresenting the orientation
 * of collisions.
 * Implements:
 * \[ Z^p_0 = \delta^*_{pq} * r^p_{eq} +
 *            \frac{5}{2} \left[1/(\delta^*_{pq} r^p_{eq})^2 \right] \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_0(double del_star, double r_eq)
{
    double tmp_a;
    
    tmp_a = del_star * r_eq;

    return ( tmp_a + 5.0/2.0 * ( 1.0 / (tmp_a * tmp_a)) );
}

/**
 * \brief Calculate Z_V.
 *
 * Z_V is the vibrational factor.
 * Implements:
 * \[ Z^p_V = \frac{f^p_m}{\pi^2} \frac{\mu_{pp}}{\mu_{pq}}
 *            \frac{\alpha_{pq}}{\theta^p_v}
 *            \left( \frac{k \theta^p_v}{\Delta E} \right)^2 \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_V(double f_m, double mu_pp, double mu_pq, double alpha_pq,
	       double theta_v, double del_E)
{
    double tmp_a, tmp_b, tmp_c, tmp_d;
    tmp_a = f_m / (1.0*M_PI * M_PI);
    tmp_b = mu_pp / mu_pq;
    tmp_c = alpha_pq / theta_v;
    tmp_d = PC_k_SI * theta_v / del_E;
    return ( tmp_a * tmp_b * tmp_c * pow(tmp_d, 2) );
}

/**
 * \brief Calculate Z_T
 *
 * Z_T is the translational factor.
 * Implements:
 * \[ Z^{pq}_T = \pi^2 \left( \frac{3}{2\pi} \right)^{1/2}
 *               \left( \frac{\Delta E}{k \alpha_{pq}} \right)^2
 *               \left( \frac{T}{\alpha_{pq}} \right)^{1/6}
 *               \times \exp \left[ \frac{3}{2}
 *                                  \left( \frac{\alpha_{pq}}{T} \right)^{1/3}
 *                                   - \frac{\Delta E}{2kT} \right] \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_T(double del_E, double alpha_pq, double T)
{
    double tmp_a, tmp_b, tmp_c;

    tmp_a = del_E / (PC_k_SI * alpha_pq);
    tmp_b = 3.0/2.0 * pow(alpha_pq/T, 1.0/3.0) - (del_E / (2.0*PC_k_SI*T));
    tmp_c = sqrt(3.0/(2.0*M_PI)) * tmp_a * tmp_a * pow(T/alpha_pq, 1.0/6.0);

    return (M_PI*M_PI * tmp_c * exp(tmp_b));
}

/**
 * \brief Calculate Z_plus.
 *
 * Z_plus represents the contribution from attractive forces.
 * This function implements:
 * \[ Z^{pq}_{+} = \exp \left[ - \frac{4}{\pi} \left(\frac{\epsilon^0_{pq} \chi^{*}_{pq}}{kT} \right)^{1/2}
 *                             - \frac{16}{3\pi^2} \frac{ \epsilon^0_{pq}}{kT} \right] \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_plus(double eps0, double chi, double T)
{
    double tmp_a, tmp_b;

    tmp_a = (-4.0 / M_PI) * sqrt(eps0*chi / (PC_k_SI*T));
    tmp_b = (-16.0*eps0) / (3.0*M_PI*M_PI*PC_k_SI*T);

    return ( exp(tmp_a + tmp_b) );
}

/**
 * \brief Calculate A.
 *
 * The scripted A represents a collision cross-reference
 * factor.
 * This function implements:
 * \[ \cal{A} = \left( \frac{r_{c*}_{pq}}{\sigma_{pq}} \right)^2 \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_A_factor(double r_c, double sigma)
{
    return ( (r_c*r_c) / (sigma*sigma) );
}

VV_SSH::
VV_SSH( lua_State * L )
    : Relaxation_time()
{
    // 1. Vibrating species 'p'
    string p_name = get_string(L, -1, "p_name");
    Diatomic_species * P = get_library_diatom_pointer_from_name( p_name );
    ip_ = P->get_isp();
    iT_ = P->get_mode_pointer_from_type("translation")->get_iT();
    iTvp_ = P->get_mode_pointer_from_type("vibration")->get_iT();
    // FIXME: Need to test here that p_vib is truncated harmonic vibrational mode
    Truncated_harmonic_vibration * p_vib = dynamic_cast<Truncated_harmonic_vibration*>(P->get_mode_pointer_from_type("vibration"));
    theta_v_p_ = p_vib->get_theta();
    M_p_ = P->get_M();
    r_eq_p_ = P->get_r_eq();
    f_m_p_ = P->get_f_m();
    mu_p_ = P->get_mu();
    
    // 2. Vibrating species 'q'
    string q_name = get_string(L, -1, "q_name");
    Diatomic_species * Q = get_library_diatom_pointer_from_name( q_name );
    iq_ = Q->get_isp();
    iTvq_ = Q->get_mode_pointer_from_type("vibration")->get_iT();
    // FIXME: Need to test here that q_vib is truncated harmonic vibrational mode
    Truncated_harmonic_vibration * q_vib = dynamic_cast<Truncated_harmonic_vibration*>(Q->get_mode_pointer_from_type("vibration"));
    theta_v_q_ = q_vib->get_theta();
    M_q_ = Q->get_M();
    r_eq_q_ = Q->get_r_eq();
    f_m_q_ = Q->get_f_m();
    mu_q_ = Q->get_mu();
    
    // 3. Derived parameters
    if( P->get_polar_flag() ) {
	if( Q->get_polar_flag() ) {
	    r_ = collider_distance_A(*P, *Q);
	    eps_ = potential_well_A(*P, *Q);
	}
	else {
	    r_ = collider_distance_B(*P, *Q);
	    eps_ = potential_well_B(*P, *Q);
	}
    }
    else {
	if( Q->get_polar_flag() ) {
	    r_ = collider_distance_B(*Q, *P);
	    eps_ = potential_well_B(*Q, *P);
	}
	else {
	    r_ = collider_distance_A(*Q, *P);
	    eps_ = potential_well_A(*Q, *P);
	}
    }

    if( p_name == "N2" && q_name == "N2" ) {
    	cout << "VV_SSH::VV_SSH()" << endl
	     << "Applying the correction to the Lennard-Jones potential parameter, r,\n"
	     << "for N2-N2 interactions as suggested by Thivet et al. (1990).\n";
	r_ = 4.2e-10;
	cout << "r_0(N2,N2)= " << r_ << endl;
    }

    mu_ = ((M_p_ * M_q_) / (M_p_ + M_q_)) / PC_Avogadro;
    delta_E_ = PC_k_SI * (theta_v_p_ - theta_v_q_);
    sigma_ = r_;
}

VV_SSH::
~VV_SSH()
{}

double
VV_SSH::
VV_SSH_transition_probability( double T )
{
    /* 1. Values that change with T */    
    double beta = SSH_beta(eps_, mu_, r_, delta_E_, T);
    double r_c_star = SSH_r_c_star(beta, r_);
    double del_star = SSH_del_star(beta, r_);
    double alpha_pq = SSH_alpha_pq(mu_, delta_E_, del_star);
    double chi_pq = SSH_chi_pq(alpha_pq, T);

    /* 2, The steric factors */
    double Z_0_p = SSH_Z_0(del_star, r_eq_p_);
    double Z_0_q = SSH_Z_0(del_star, r_eq_q_);
    double Z_V_p = SSH_Z_V(f_m_p_, mu_p_, mu_, alpha_pq, theta_v_p_, delta_E_);
    double Z_V_q = SSH_Z_V(f_m_q_, mu_q_, mu_, alpha_pq, theta_v_q_, delta_E_);
    double Z_T = SSH_Z_T(delta_E_, alpha_pq, T);
    double Z_plus = SSH_Z_plus(eps_, chi_pq, T);

    /* collision cross-reference factor */

    double A = SSH_A_factor(r_c_star, sigma_);
    double P = A / (Z_0_p * Z_0_q * Z_V_p * Z_V_q * Z_T * Z_plus);

    return (P);
}

double
VV_SSH::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double T = Q.T[iT_];

    // If either collider has zero concentration
    // then a relaxation time has no meaning
    if( molef[ip_] <= 0.0 || molef[iq_] <= 0.0 ) {
	return -1.0;
    }

    double nd = 0.0;	// colliders / m**3
    double Z = 0.0;	// collision frequency
    
    if( ip_ == iq_ ) {
	nd = Q.massf[ip_] * Q.rho / M_p_ * PC_Avogadro;
	Z = collision_frequency(sigma_, mu_, T, nd);
    }
    else {
	nd = ( Q.massf[ip_] / M_p_ + Q.massf[iq_] / M_q_ ) * Q.rho * PC_Avogadro;
	Z = collision_frequency(sigma_, mu_, T, nd);
    }

    // 3. Compute the probability of transition
    double P = VV_SSH_transition_probability(T);

    double P_rev = P * exp(1.0*(theta_v_p_ - theta_v_q_)/T);

    // 4. Compute the relaxation time
    double tau = 1.0 / (Z * P_rev);
    //cout << "VV_SSH_time:: p= " << p_->ivib << " q= " << q_->ivib << " tau= " << tau << endl;

    return tau;
}

VV_Candler::
VV_Candler( lua_State * L )
    : Relaxation_time()
{
    // 1. Vibrating species 'p'
    string p_name = get_string(L, -1, "p_name");
    Diatomic_species * P = get_library_diatom_pointer_from_name( p_name );
    ip_ = P->get_isp();
    iT_ = P->get_mode_pointer_from_type("translation")->get_iT();
    double M_p = P->get_M();
    double r0_p = P->get_r0();
    
    // 2. Vibrating species 'q'
    string q_name = get_string(L, -1, "q_name");
    Diatomic_species * Q = get_library_diatom_pointer_from_name( q_name );
    iq_ = Q->get_isp();
    M_q_ = Q->get_M();
    double r0_q = Q->get_r0();
    
    // 3. derived parameters
    mu_ = ((M_p * M_q_) / (M_p + M_q_)) / PC_Avogadro;
    // CHECKME: Candler uses d_p*dq here. is there a difference between r0 and d? 
    sigma_ = r0_p * r0_q;	// collision cross-section
    
    // 4. Transition probability (use fixed value as recommended in Candler thesis)
    P_ = 1.0e-2;
}

VV_Candler::
~VV_Candler()
{}

double
VV_Candler::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    // If either collider has zero concentration
    // then a relaxation time has no meaning
    if( molef[ip_] <= 0.0 || molef[iq_] <= 0.0 ) {
	return -1.0;
    }
    
    // Calculate collision frequency x average vibrational energy
    double Z_eps = PC_Avogadro * sigma_ * sqrt( 8.0 * PC_R_u * Q.T[iT_] / ( M_PI * mu_ ) ) 
    			* Q.rho * Q.massf[iq_] / M_q_;
    
    // Calculate relaxation time
    double tau = 1.0 / ( P_ * Z_eps );
    
    return tau;
}

VE_Lee::
VE_Lee( lua_State * L )
    : Relaxation_time()
{
    // 1. Free electron and vibrational species index (and iTe)
    Chemical_species * E = get_library_species_pointer_from_name( "e_minus" );
    ie_ = E->get_isp();
    iTe_ = E->get_mode_pointer_from_type("translation")->get_iT();
    string v_name = get_string(L,-1,"v_name");
    iv_ = get_library_species_pointer_from_name( v_name )->get_isp();
    
    // 2. Temperature switches
    lua_getfield(L, -1, "T_switches" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VE_Lee::VE_Lee():\n";
	ost << "Error in the declaration of T switches: a table is expected.\n";
	input_error(ost);
    }
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	T_switches_.push_back( luaL_checknumber(L,-1) );
    	lua_pop(L,1);
    }
    lua_pop(L,1);	// pop 'T_switches'  
    
    // 3. p.tau coefficients
    ptau_coeffs_.resize( T_switches_.size() + 1 );
    for ( size_t ic=0; ic<T_switches_.size()+1; ++ic ) {
    	ostringstream ptau_label;
    	ptau_label << "ptau_coefficients_" << ic;
	lua_getfield(L, -1, ptau_label.str().c_str() );
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "VE_Lee::VE_Lee():\n";
	    ost << "Error in the declaration of " << ptau_label.str() << ": a table is expected.\n";
	    input_error(ost);
	}
	int nCs = lua_objlen(L, -1);
	if ( nCs!=3 ) {
	    ostringstream ost;
	    ost << "VE_Lee::VE_Lee():\n";
	    ost << "Error in the declaration of p.tau coefficients: 3 coefficients are expected.\n";
	    input_error(ost);
	}
	for ( int i=0; i<nCs; ++i ) {
	    lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
	    ptau_coeffs_[ic].push_back( luaL_checknumber(L,-1) );
	    lua_pop(L,1);
	}
	lua_pop(L,1);	// pop 'ptau_coefficients'
    }
}

VE_Lee::
~VE_Lee()
{
    T_switches_.resize(0);
    ptau_coeffs_.resize(0);
}

double
VE_Lee::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    // If either collider has zero concentration
    // then a relaxation time has no meaning
    if( molef[ie_] <= 0.0 || molef[iv_] <= 0.0 ) {
	return -1.0;
    }
    
    // Select the ptau coefficients to use from Te
    double Te = Q.T[iTe_];
    size_t iptc;
    for ( iptc=0; iptc<T_switches_.size(); ++iptc ) {
    	if ( Te < T_switches_[iptc] ) break;
    }
    
    // Calculate tau
    // NOTE: pe_atm should be > 0 as molef[ie_] > 0
    double pe_atm = Q.p_e / PC_P_atm;
    double log_ptau = ptau_coeffs_[iptc][0] * log10(Te) * log10(Te) + ptau_coeffs_[iptc][1] * log10(Te) + ptau_coeffs_[iptc][2];
    double tau = pow( 10.0, log_ptau ) / pe_atm;
    
    return tau;
}

RT_Parker::
RT_Parker( lua_State * L )
    : Relaxation_time()
{
    // 1. Rotating species data
    string p_name = get_string(L, -1, "p_name");
    Chemical_species * p = get_library_species_pointer_from_name( p_name );
    double M_p = p->get_M();
    ip_ = p->get_isp();
    iT_ = p->get_iT_trans();
    
    // 2. Colliding species data
    lua_getfield(L, -1, "q_names" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "RT_Parker::RT_Parker():\n";
	ost << "Error in the declaration of colliding species : a table 'q_names' is expected.\n";
	input_error(ost);
    }
    int nqs = lua_objlen(L, -1);
    for ( int i=0; i<nqs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	const char* q_name = luaL_checkstring(L, -1);
    	lua_pop(L, 1);
    	Chemical_species * X = get_library_species_pointer_from_name( q_name );
    	M_qs_.push_back( X->get_M() );
    	iqs_.push_back( X->get_isp() );
    	// CHECKME: SI units?
    	mus_.push_back( ((M_p * M_qs_.back()) / (M_p + M_qs_.back())) );
    	// NOTE: using constant value from Abe et al (2002) for collision cross-section
    	sigmas_.push_back( 1.0e-19 );
    }
    
    lua_pop(L,1);	// pop q_names
}

RT_Parker::
~RT_Parker() {}

double
RT_Parker::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double n_q;		// colliding species number density
    double tau_q_inv;
    double tau_inv = 0.0;
    
    // NOTE: assuming thermally perfect gases
    for ( size_t i=0; i<iqs_.size(); ++i ) {
	n_q = Q.rho * Q.massf[iqs_[i]] / M_qs_[i] * PC_Avogadro;
	tau_q_inv = sigmas_[i] * n_q * sqrt( 8.0 * PC_k_SI * Q.T[iT_] / ( M_PI * mus_[i] ) );
	tau_inv += tau_q_inv;
	// cout << "n_q = " << n_q << ", tau_q_inv = " << tau_q_inv << endl;
    }

    return 1.0/tau_inv;
}

RE_Abe::
RE_Abe( lua_State * L )
    : Relaxation_time()
{
    // Initialise ions
    lua_getfield(L,-1,"ions");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "RE_Abe::RE_Abe():\n";
	ost << "Error in declaration of ions: a table is expected.\n";
	input_error(ost);
    }
    int n_ions = lua_objlen(L, -1);
    for ( int i=0; i<n_ions; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	tau_REs_.push_back( new RE_Ion(L) );
    	lua_pop(L, 1);
    }
    lua_pop(L, 1);	// pop ions
    
    // Initialise neutrals
    lua_getfield(L,-1,"neutrals");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "RE_Abe::RE_Abe():\n";
	ost << "Error in declaration of neutrals: a table is expected.\n";
	input_error(ost);
    }
    int n_neutrals = lua_objlen(L, -1);
    for ( int i=0; i<n_neutrals; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	tau_REs_.push_back( new RE_Neutral(L) );
    	lua_pop(L, 1);
    }
    lua_pop(L, 1);	// pop neutrals
}

RE_Abe::
~RE_Abe()
{
    for ( size_t i=0; i<tau_REs_.size(); ++i )
    	delete tau_REs_[i];
}

double
RE_Abe::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double tau_inv = 0.0;
    for ( size_t i=0; i<tau_REs_.size(); ++i ) {
    	double tau = tau_REs_[i]->compute_relaxation_time(Q,molef);
    	if ( tau < 0.0 ) continue;
    	tau_inv += 1.0/tau;
    }
    
    return 1.0/tau_inv;
}

RE_Ion::
RE_Ion( lua_State * L )
    : Relaxation_time()
{
    // 1. Colliding species data
    string c_name = get_string(L, -1, "c_name");
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    if ( X->get_Z()<1 ) {
	ostringstream ost;
	ost << "RE_Ion::RE_Ion():\n";
	ost << "Error in the declaration of colliding species: " << c_name << " is not an ion.\n";
	input_error(ost);
    }    
    ic_ = X->get_isp();
    M_c_ = X->get_M();
    // Initialise g_rot factor (see Abe et al 2002 paper)
    if ( c_name=="N2_plus" )
    	g_rot_ = 10.0;
    else if ( c_name=="O2_plus" )
    	g_rot_ = 10.0;
    else if ( c_name=="NO_plus" )
    	g_rot_ = 100.0;
    else {
	ostringstream ost;
	ost << "RE_Ion::RE_Ion():\n";
	ost << "Species: " << c_name << " does not have a known g_rot value.\n";
	ost << "Setting g_rot to 10 (N2_plus value)\n";
	// input_error(ost);
	cout << ost.str();
	g_rot_ = 10.0;
    }
    
    // 2. Electron data
    X = get_library_species_pointer_from_name( "e_minus" );
    ie_ = X->get_isp();
    iTe_ = X->get_iT_trans();
    M_e_ = X->get_M();
}

RE_Ion::
~RE_Ion() {}

double
RE_Ion::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;	// Number density of colliders (in cm**-3)
    double n_e = Q.massf[ie_] * Q.rho / M_e_ * PC_Avogadro * 1.0e-6;	// Number density of electrons (in cm**-3)
    
    // If there are no participating species set a negative relaxation time
    if ( n_e==0.0 || n_c==0.0 ) return -1.0;

    double tmpA = 8.0 / 3.0 * sqrt( M_PI / PC_m_CGS );
    double tmpB = n_c * pow ( PC_e_CGS, 4.0 ) / pow ( 2.0 * PC_k_CGS * Q.T[iTe_], 1.5 );
    double tmpC = log( pow( PC_k_CGS * Q.T[iTe_], 3.0 ) / ( M_PI * n_e * pow ( PC_e_CGS, 6.0 ) ) );
    double nu_ec = tmpA * tmpB * tmpC;
    
    // NOTE: applying g_rot factor to denominator
    double tau_ec = M_c_ / nu_ec / g_rot_;
    
    return tau_ec;
    
}

RE_Neutral::
RE_Neutral( lua_State * L )
    : Relaxation_time()
{
    // 1. Colliding species data
    string c_name = get_string(L, -1, "c_name");
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    if ( X->get_Z()!=0 ) {
	ostringstream ost;
	ost << "ET_Neutral::ET_Neutral():\n";
	ost << "Error in the declaration of colliding species: " << c_name << " is not a neutral.\n";
	input_error(ost);
    }    
    ic_ = X->get_isp();
    M_c_ = X->get_M();
    // Initialise g_rot factor (see Abe et al 2002 paper)
    if ( c_name=="N2" )
    	g_rot_ = 10.0;
    else if ( c_name=="O2" )
    	g_rot_ = 10.0;
    else if ( c_name=="NO" )
    	g_rot_ = 100.0;
    else {
	ostringstream ost;
	ost << "RE_Ion::RE_Ion():\n";
	ost << "Species: " << c_name << " does not have a known g_rot value.\n";
	ost << "Setting g_rot to 10 (N2 value)\n";
	// input_error(ost);
	cout << ost.str();
	g_rot_ = 10.0;
    }
    
    // 2. Electron data
    X = get_library_species_pointer_from_name( "e_minus" );
    iTe_ = X->get_iT_trans();
    
    // 3. Sigma quadratic curve fit coefficients
    lua_getfield(L, -1, "sigma_coefficients" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "RE_Neutral::RE_Neutral():\n";
	ost << "Error in the declaration of sigma coefficients: a table is expected.\n";
	input_error(ost);
    }
    int nCs = lua_objlen(L, -1);
    if ( nCs!=3 ) {
	ostringstream ost;
	ost << "RE_Neutral::RE_Neutral():\n";
	ost << "Error in the declaration of sigma coefficients: 3 coefficients are expected.\n";
	input_error(ost);
    }
    for ( int i=0; i<nCs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	C_.push_back( luaL_checknumber(L,-1) );
    	lua_pop(L,1);
    }
    lua_pop(L,1);	// pop 'sigma_coefficients'    	
}

RE_Neutral::
~RE_Neutral()
{
    C_.resize(0);
}

double
RE_Neutral::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;	// Number density of colliders (in cm**-3)
    
    // If there are no colliding species set a very large relaxation time
    if ( n_c==0.0 ) return 1.0e99;

    double T = Q.T[iTe_];
    double simga_ec = C_[0] + C_[1] * T + C_[2] * T * T;
    double nu_ec = n_c * simga_ec * sqrt( 8.0 * PC_k_CGS * T / ( M_PI * PC_m_CGS ) );
    // NOTE: applying g_rot factor to denominator
    double tau_ec = M_c_ / nu_ec / g_rot_;
    
    return tau_ec;
}

Relaxation_time* create_new_relaxation_time( lua_State * L )
{
    // 1. get name of relaxation time type to be created
    string relaxation_time = get_string(L, -1, "type");

    // 2. create the relaxation time instance
    if( relaxation_time == "VT_MillikanWhite" ) {
	return new VT_MillikanWhite(L);
    }
    else if( relaxation_time == "VT_MillikanWhite_HTC" ) {
	return new VT_MillikanWhite_HTC(L);
    }
    else if( relaxation_time == "ET_AppletonBray" ) {
	return new ET_AppletonBray(L);
    }
    else if( relaxation_time == "ET_Ion" ) {
	return new ET_Ion(L);
    }
    else if( relaxation_time == "ET_Neutral" ) {
	return new ET_Neutral(L);
    }
    else if( relaxation_time == "VV_SSH" ) {
	return new VV_SSH(L);
    }
    else if( relaxation_time == "VV_Candler" ) {
	return new VV_Candler(L);
    }
    else if( relaxation_time == "VE_Lee" ) {
	return new VE_Lee(L);
    }
    else if( relaxation_time == "RT_Parker" ) {
	return new RT_Parker(L);
    }
    else if( relaxation_time == "RE_Abe" ) {
	return new RE_Abe(L);
    }
    else if( relaxation_time == "RE_Ion" ) {
	return new RE_Ion(L);
    }
    else if( relaxation_time == "RE_Neutral" ) {
	return new RE_Neutral(L);
    }
    else {
    	cout << "create_new_relaxation_time()\n";
    	cout << "The relaxation_time: " << relaxation_time << " is not known.\n";
	exit(BAD_INPUT_ERROR);
    }
}

