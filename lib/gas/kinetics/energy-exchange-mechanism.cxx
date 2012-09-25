// Author: Daniel F. Potter
// Date: 18-Nov-2009

#include <iostream>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "energy-exchange-mechanism.hh"
#include "../models/chemical-species-library.hh"

using namespace std;

Energy_exchange_mechanism::
Energy_exchange_mechanism() {}

Energy_exchange_mechanism::
~Energy_exchange_mechanism() {}

double
Energy_exchange_mechanism::
compute_rate(const valarray<double> &y, Gas_data &Q, std::vector<double> &molef)
{
    return specific_compute_rate(y, Q, molef);
}

double
Energy_exchange_mechanism::
py_compute_rate(const vector<double> &y, Gas_data &Q, vector<double> &molef)
{
  valarray<double> y_(y.size());
  for ( size_t i=0; i<y.size(); ++i)
    y_[i] = y[i];

  return this->compute_rate(y_,Q,molef);
}

VT_exchange::
VT_exchange(lua_State *L, int ip, int imode)
    : Energy_exchange_mechanism(), ip_(ip), iTv_(imode)
{
    Chemical_species *X = get_library_species_pointer(ip_);
    p_vib_ = X->get_mode_pointer_from_type("vibration");

    int iq = get_int(L, -1, "iq");
    iT_ = get_int(L, -1, "itrans");

    lua_getfield(L,-1,"relaxation_time");
    tau_VT_ = create_new_relaxation_time(L, ip, iq, iT_);
    lua_pop(L, 1 );
}

VT_exchange::
~VT_exchange()
{
    delete tau_VT_;
}

double
VT_exchange::
specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, vector<double> &molef)
{
    // tau_ will be present and correct before beginning this
    // function ie. a call to compute_tau is expected earlier.
    double e_vib_eq = p_vib_->eval_energy_from_T(Q.T[iT_]);
    double e_vib = p_vib_->eval_energy_from_T(Q.T[iTv_]);
    // NOTE: - tau_ needs to be (already) weighted by colliding mole fractions
    //       - massf scaling is to convert J/s/kg-of-species-ip to J/s/kg-of-mixture
    double rate = Q.massf[ip_] * (e_vib_eq - e_vib) / tau_;
    //    cout << "e_vib_eq= " << e_vib_eq << endl;
    //    cout << "e_vib= " << e_vib << endl;
    //    cout << "massf= " << Q.massf[ip_] << endl;
    //    cout << "rate= " << rate << endl;

    return rate;
}

// Polymodal_VT_exchange::
// Polymodal_VT_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Get vibrating species pointer
//     string p_name = get_string(L, -1, "p_name");
//     Chemical_species * X = get_library_species_pointer_from_name( p_name );
    
//     // 2. Vibrating species isp index
//     ip_ = X->get_isp();
    
//     // 3. Translational energy mode index
//     iT_ = X->get_mode_pointer_from_type("translation")->get_iT();
    
//     // 4. Pointers to vibrational modes
//     for ( int imode=0; imode<X->get_n_modes(); ++imode ) {
//     	if ( X->get_mode_pointer(imode)->get_type()=="vibration" ) {
//     	    p_vib_.push_back( X->get_mode_pointer(imode) );
//     	}
//     }
    
//     // 5. Vibrational energy mode index
//     iTv_ = p_vib_[0]->get_iT();
    
//     // 6. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_VT_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// Polymodal_VT_exchange::
// ~Polymodal_VT_exchange()
// {
//     delete tau_VT_;
// }

// double
// Polymodal_VT_exchange::
// specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // tau_ will be present and correct before beginning this
//     // functionm ie. a call to compute_tau is expected earlier.

//     double e_vib_eq = 0.0;
//     double e_vib = 0.0;
//     for ( size_t ivib=0; ivib<p_vib_.size(); ++ivib ) {
//     	e_vib_eq += p_vib_[ivib]->eval_energy_from_T(Q.T[iT_]);
//         e_vib += p_vib_[ivib]->eval_energy_from_T(Q.T[iTv_]);
//     }
    
//     // NOTE: - tau_ needs to be (already) weighted by colliding mole fractions
//     //       - massf scaling is to convert J/s/kg-of-species-ip to J/s/kg-of-mixture
//     double rate = Q.massf[ip_] * (e_vib_eq - e_vib) / tau_;

//     return rate;
// }

// ET_exchange::
// ET_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Assume translational index is 0
//     iT_ = 0;
    
//     // 2. Initialise electron from name
//     Chemical_species * e = get_library_species_pointer_from_name( "e_minus" );
//     ie_ = e->get_isp();
//     iTe_ = e->get_mode_pointer_from_type("translation")->get_iT();;

//     // 3. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_ET_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// ET_exchange::
// ~ET_exchange()
// {
//     delete tau_ET_;
// }

// double
// ET_exchange::
// specific_compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // tau_ET_ will be present and correct before beginning this
//     // functionm ie. a call to compute_tau is expected earlier.
//     // CHECKME: - from Abe and Panesi, maybe we should be scaling by molef[ie_]
//     //            and NOT massf[ie_]?
//     // NOTE: - scaling by molef gives a way to fast E-T equilibriation
//     double rate = Q.massf[ie_] * 3.0 * PC_R_u * ( Q.T[iT_] - Q.T[iTe_] ) / tau_;
    
//     // cout << "ET_exchange::specific_compute_rate()" << endl
//     //      << "rate = " << rate << endl;
    
//     return rate;
// }

// VV_THO_exchange::
// VV_THO_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Vibrating species 'p'
//     // 1a. Species pointer
//     string p_name = get_string(L, -1, "p_name");
//     Chemical_species * P = get_library_species_pointer_from_name( p_name );
//     // 1b. Thermal mode indices
//     iT_ = P->get_mode_pointer_from_type("translation")->get_iT();
//     iTvp_ = P->get_mode_pointer_from_type("vibration")->get_iT();
//     // 1c. Pointer to vibrational mode cast as THO
//     // FIXME: Need to test here that p_vib is truncated harmonic vibrational mode
//     p_vib_ = dynamic_cast<Truncated_harmonic_vibration*>(P->get_mode_pointer_from_type("vibration"));
//     // 1d. Characteristic vibrational temperature for p
//     theta_v_p_ = p_vib_->get_theta();
    
//     // 2. Vibrating species 'q'
//     // 2a. Species pointer
//     string q_name = get_string(L, -1, "q_name");
//     Chemical_species * Q = get_library_species_pointer_from_name( q_name );
//     // 2b. Species Q species index
//     iq_ = Q->get_isp();
//     // 2c. Thermal mode indices
//     iTvq_ = Q->get_mode_pointer_from_type("vibration")->get_iT();
//     // 2d. Pointer to vibrational mode cast as THO
//     // FIXME: Need to test here that q_vib is truncated harmonic vibrational mode
//     q_vib_ = dynamic_cast<Truncated_harmonic_vibration*>(Q->get_mode_pointer_from_type("vibration"));
//     // 1d. Characteristic vibrational temperature for q
//     theta_v_q_ = q_vib_->get_theta();
    
//     // 3. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_VV_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// VV_THO_exchange::
// ~VV_THO_exchange()
// {
//     delete tau_VV_;
// }

// double
// VV_THO_exchange::
// specific_compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // See VV_exchange::specific_compute_rate() in mTg_mechanisms.cxx for
//     // Rowan's initial implementation of Thivet's equations
    
//     // 0. Set the translational temperature
//     double T = Q.T[iT_];
    
//     // 1. Calculate vibrational energies
//     double e_p = p_vib_->eval_energy_from_T( Q.T[iTvp_] );	// p THO energy at Tv_p
//     double e_q = p_vib_->eval_energy_from_T( Q.T[iTvp_] );	// q THO energy at Tv_q
//     double e_p_bar = p_vib_->eval_energy_from_T( T );		// p THO energy at T
//     double e_q_bar = q_vib_->eval_energy_from_T( T );		// q THO energy at T
//     double e_q_hat = q_vib_->eval_HO_energy_from_T( T );	// q HO energy at T
    
//     // 2. Calculate the rate
//     double tmp_a = molef[iq_] / tau_;
//     double tmp_b = (1.0 - exp(-1.0 * theta_v_p_/T)) / (1.0 - exp(-1.0 * theta_v_q_/T))
// 	* ((e_q/e_q_hat)*(e_p_bar - e_p));
//     double tmp_c = (e_p/e_q_hat)*(e_q_bar - e_q);
//     double rate = tmp_a * (tmp_b - tmp_c);
    
//     return rate;
// }

// VV_HO_exchange::
// VV_HO_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Vibrating species 'p'
//     // 1a. Species pointer
//     string p_name = get_string(L, -1, "p_name");
//     Chemical_species * P = get_library_species_pointer_from_name( p_name );
//     // 1b. Species Q species index
//     ip_ = P->get_isp();
//     // 1c. Thermal mode indices
//     iT_ = P->get_mode_pointer_from_type("translation")->get_iT();
//     iTvp_ = P->get_mode_pointer_from_type("vibration")->get_iT();
//     // 1d. Pointer to vibrational mode cast as THO
//     // FIXME: Need to test here that p_vib is a harmonic vibrational mode
//     p_vib_ = dynamic_cast<Harmonic_vibration*>(P->get_mode_pointer_from_type("vibration"));
    
//     // 2. Vibrating species 'q'
//     // 2a. Species pointer
//     string q_name = get_string(L, -1, "q_name");
//     Chemical_species * Q = get_library_species_pointer_from_name( q_name );
//     // 2b. Thermal mode indices
//     iTvq_ = Q->get_mode_pointer_from_type("vibration")->get_iT();
//     // 2c. Pointer to vibrational mode cast as THO
//     // FIXME: Need to test here that q_vib is a harmonic vibrational mode
//     q_vib_ = dynamic_cast<Harmonic_vibration*>(Q->get_mode_pointer_from_type("vibration"));
    
//     // 3. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_VV_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// VV_HO_exchange::
// ~VV_HO_exchange()
// {
//     delete tau_VV_;
// }

// double
// VV_HO_exchange::
// specific_compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // Refs: Panesi et al (2008) AIAA 2008-1205
//     //       Candler et al (1991) JTHT Vol. 5 No. 11 pp. 266
//     //       Knab et al (1994) 2nd ESASV pp 129
    
//     // 0. Set the translational temperature
//     double T = Q.T[iT_];
    
//     // 1. Calculate vibrational energies
//     double e_p = p_vib_->eval_energy_from_T( Q.T[iTvp_] );	// p THO energy at Tv_p
//     double e_q = p_vib_->eval_energy_from_T( Q.T[iTvp_] );	// q THO energy at Tv_q
//     double e_p_bar = p_vib_->eval_energy_from_T( T );		// p THO energy at T
//     double e_q_bar = q_vib_->eval_energy_from_T( T );		// q THO energy at T
    
//     // 2. Compute the rate
//     double rate = Q.massf[ip_] * ( e_p_bar * e_q / e_q_bar - e_p ) / tau_;
    
//     return rate;
// }

// VE_exchange::
// VE_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Free electron
//     Chemical_species * E = get_library_species_pointer_from_name( "e_minus" );
//     ie_ = E->get_isp();
//     iTe_ = E->get_mode_pointer_from_type("translation")->get_iT();
    
//     // 2. Vibrating species
//     string v_name = get_string(L, -1, "v_name");
//     Chemical_species * V = get_library_species_pointer_from_name( v_name );
//     iv_ = V->get_isp();
//     iTv_ = V->get_mode_pointer_from_type("vibration")->get_iT();
//     v_vib_ = V->get_mode_pointer_from_type("vibration");

//     // 3. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_VE_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// VE_exchange::
// ~VE_exchange()
// {
//     delete tau_VE_;
// }

// double
// VE_exchange::
// specific_compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // Refs: Gnoffo (1989)
//     //       Lee (?)
    
//     // 1. Calculate vibrational energies
//     double e_vib = v_vib_->eval_energy_from_T( Q.T[iTv_] );
//     double e_vib_bar = v_vib_->eval_energy_from_T( Q.T[iTe_] );
    
//     // 2. Compute the rate
//     double rate = Q.massf[iv_] * ( e_vib_bar - e_vib ) / tau_;
    
//     // cout << "VE_exchange::specific_compute_rate()" << endl
//     //      << "rate = " << rate << endl;
//     // Q.print(false);
    
//     return rate;
// }

// EV_exchange::
// EV_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Free electron
//     Chemical_species * E = get_library_species_pointer_from_name( "e_minus" );
//     ie_ = E->get_isp();
//     iTe_ = E->get_mode_pointer_from_type("translation")->get_iT();
    
//     // 2. Vibrating species
//     string v_name = get_string(L, -1, "v_name");
//     Chemical_species * V = get_library_species_pointer_from_name( v_name );
//     iv_ = V->get_isp();
//     iTv_ = V->get_mode_pointer_from_type("vibration")->get_iT();
//     v_vib_ = V->get_mode_pointer_from_type("vibration");

//     // 3. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_EV_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// EV_exchange::
// ~EV_exchange()
// {
//     delete tau_EV_;
// }

// double
// EV_exchange::
// specific_compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // Refs: Gnoffo (1989)
//     //       Lee (?)
    
//     // 1. Calculate vibrational energies
//     double e_vib = v_vib_->eval_energy_from_T( Q.T[iTv_] );
//     double e_vib_bar = v_vib_->eval_energy_from_T( Q.T[iTe_] );
    
//     // 2. Compute the rate (note e_vib - e_vib_bar )
//     double rate = Q.massf[iv_] * ( e_vib - e_vib_bar ) / tau_;
    
//     // cout << "EV_exchange::specific_compute_rate()" << endl
//     //      << "rate = " << rate << endl;
//     // Q.print(false);
    
//     return rate;
// }

// RT_exchange::
// RT_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Get rotating species pointer
//     string p_name = get_string(L, -1, "p_name");
//     Chemical_species * X = get_library_species_pointer_from_name( p_name );
    
//     // 2. Value for Z_inf_ from paper of Abe et al (2002)
//     if ( p_name == "N2" ) {
//     	Z_R_inf_ = 15.7;
//     }
//     else if ( p_name == "O2" || p_name == "NO" ) {
//     	Z_R_inf_ = 14.4;
//     }
//     else {
// 	ostringstream ost;
// 	ost << "RT_exchange::RT_exchange()\n";
// 	ost << "Species: " << p_name << " does not have a known Z_R_inf value." << endl;
// 	ost << "Setting Z_R_inf to 14.4 (O2 value)" << endl;
// 	// input_error(ost);
// 	cout << ost.str();
// 	Z_R_inf_ = 14.4;
//     }
    
//     // 3. Rotating species isp index
//     ip_ = X->get_isp();
    
//     // 4. Translational energy mode index
//     iT_ = X->get_mode_pointer_from_type("translation")->get_iT();
    
//     // 5. Pointer to rotational mode
//     p_rot_ = X->get_mode_pointer_from_type("rotation");
    
//     // 6. Rotational energy mode index
//     iTr_ = p_rot_->get_iT();
    
//     // 7. Characteristic rotational temperature
//     theta_r_ = dynamic_cast<Fully_excited_rotation*>(p_rot_)->get_theta();
    
//     // 8. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_RT_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );

// }

// RT_exchange::
// ~RT_exchange()
// {
//     delete tau_RT_;
// }

// static const double half_Pi_pow_3_2 = 2.7841639984158539;
// static const double quarter_Pi_pow_2_plus_Pi = 4.5336746527977203;

// double 
// RT_exchange::
// calculate_collision_number( double T )
// {
//     double denom = 1.0 + half_Pi_pow_3_2 * sqrt( theta_r_ / T ) + 
//     			quarter_Pi_pow_2_plus_Pi * ( theta_r_ / T );
    
//     return Z_R_inf_ / denom;
// }

// double
// RT_exchange::
// specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // tau_ will be present and correct before beginning this
//     // functionm ie. a call to compute_tau is expected earlier.

//     double e_rot_eq = p_rot_->eval_energy_from_T(Q.T[iT_]);
//     double e_rot = p_rot_->eval_energy_from_T(Q.T[iTr_]);
//     double Z_R  = calculate_collision_number(Q.T[iT_]);
//     // cout << "Z_R = " << Z_R << endl;
//     // cout << "tau_ = " << tau_ << endl;
//     // NOTE: - tau_ needs to be (already) weighted by colliding mole fractions
//     //       - massf scaling is to convert J/s/kg-of-species-ip to J/s/kg-of-mixture
//     double rate = Q.massf[ip_] * (e_rot_eq - e_rot) / ( Z_R * tau_ );
    
    

//     return rate;
// }

// RE_exchange::
// RE_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Initialise rotational species data
//     string q_name = get_string(L, -1, "q_name");
//     Chemical_species * X = get_library_species_pointer_from_name( q_name );
//     iTr_ = X->get_mode_pointer_from_type("rotation")->get_iT();;
    
//     // 2. Initialise electron species data
//     Chemical_species * e = get_library_species_pointer_from_name( "e_minus" );
//     ie_ = e->get_isp();
//     iTe_ = e->get_mode_pointer_from_type("translation")->get_iT();;

//     // 3. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_RE_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// RE_exchange::
// ~RE_exchange()
// {
//     delete tau_RE_;
// }

// double
// RE_exchange::
// specific_compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // tau_ should be present and correct before beginning this
//     // function ie. a call to compute_tau is expected earlier.
//     // CHECKME: - using Gnoffo (1989) ET expression as Abe et al (2002) expression looks wrong
//     //          - from Abe and Panesi, maybe we should be scaling by molef[ie_]
//     //            and NOT massf[ie_]?
//     double rate = Q.massf[ie_] * 3.0 * PC_R_u * ( Q.T[iTe_] - Q.T[iTr_] ) / tau_;
    
//     return rate;
// }

// ER_exchange::
// ER_exchange( lua_State *L )
//     : Energy_exchange_mechanism()
// {
//     // 1. Initialise rotational species data
//     string q_name = get_string(L, -1, "q_name");
//     Chemical_species * X = get_library_species_pointer_from_name( q_name );
//     iTr_ = X->get_mode_pointer_from_type("rotation")->get_iT();;
    
//     // 2. Initialise electron species data
//     Chemical_species * e = get_library_species_pointer_from_name( "e_minus" );
//     ie_ = e->get_isp();
//     iTe_ = e->get_mode_pointer_from_type("translation")->get_iT();;

//     // 3. Initialise relaxation time
//     lua_getfield(L,-1,"relaxation_time");
//     tau_RE_ = create_new_relaxation_time( L );
//     lua_pop(L, 1 );
// }

// ER_exchange::
// ~ER_exchange()
// {
//     delete tau_RE_;
// }

// double
// ER_exchange::
// specific_compute_rate(const valarray<double> &y, Gas_data &Q, vector<double> &molef)
// {
//     // tau_ will be present and correct before beginning this
//     // functionm ie. a call to compute_tau is expected earlier.
//     // CHECKME: - using Gnoffo (1989) ET expression as Abe et al (2002) expression looks wrong
//     //          - from Abe and Panesi, maybe we should be scaling by molef[ie_]
//     //            and NOT massf[ie_]?
//     double rate = Q.massf[ie_] * 3.0 * PC_R_u * ( Q.T[iTr_] - Q.T[iTe_] ) / tau_;
    
//     return rate;
// }

Energy_exchange_mechanism* create_energy_exhange_mechanism(lua_State *L, int ip, int imode)
{
    string type = get_string(L, -1, "type");
    
    if( type == "VT" ) {
	return new VT_exchange(L, ip, imode);
    }

//     // if( type == "Polymodal_VT_exchange" ) {
//     // 	return new Polymodal_VT_exchange(L);
//     // }
//     // else if( type == "ET_exchange" ) {
//     // 	return new ET_exchange(L);
//     // }
//     // else if( type == "VV_THO_exchange" ) {
//     // 	return new VV_THO_exchange(L);
//     // }
//     // else if( type == "VV_HO_exchange" ) {
//     // 	return new VV_HO_exchange(L);
//     // }
//     // else if( type == "VE_exchange" ) {
//     // 	return new VE_exchange(L);
//     // }
//     // else if( type == "EV_exchange" ) {
//     // 	return new EV_exchange(L);
//     // }
//     // else if( type == "RT_exchange" ) {
//     // 	return new RT_exchange(L);
//     // }
//     // else if( type == "RE_exchange" ) {
//     // 	return new RE_exchange(L);
//     // }
//     // else if( type == "ER_exchange" ) {
//     // 	return new ER_exchange(L);
//     // }
     else {
     	cout << "create_energy_exhange_mechanism()" << endl
 	     << "The type given as an energy exchange mechanism: " << type << endl
 	     << "is not known or not yet implemented.\n";
 	exit(BAD_INPUT_ERROR);
     }
}

// Energy_exchange_mechanism* create_energy_exhange_mechanism_from_file( string input_file )
// {
//     lua_State *L = luaL_newstate();
//     luaL_openlibs(L);
    
//     // Parse the input file...
//     if( luaL_dofile(L, input_file.c_str()) != 0 ) {
// 	ostringstream ost;
// 	ost << "create_energy_exhange_mechanism_from_file():\n";
// 	ost << "Error in input file: " << input_file << endl;
// 	input_error(ost);
//     }
    
//     lua_getglobal(L, "mechanism");
//     if ( !lua_istable(L, -1) ) {
// 	ostringstream ost;
// 	ost << "create_energy_exhange_mechanism_from_file()\n";
// 	ost << "Error interpreting 'mechanism'; a table is expected.\n";
// 	input_error(ost);
//     }
    
//     Energy_exchange_mechanism* eem = create_energy_exhange_mechanism(L);
    
//     lua_pop(L,1);
    
//     lua_close(L);
    
//     return eem;
// }

