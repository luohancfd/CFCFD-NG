/** \file polyatomic_radiator.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 20-Feb-2013: Initial version
 *
 *  \brief Definitions for polyatomic radiator classes
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
#include "polyatomic_radiator.hh"
// #include "diatomic_system.hh"
#include "radiation_constants.hh"
#include "photaura.hh"

using namespace std;

PolyatomicElecLev::PolyatomicElecLev( int ilev, double E, int g, lua_State * L )
 : ElecLev( ilev, E, g )
{
    // get the rest of the data from the lua_State
}

string
PolyatomicElecLev::
string()
{
    ostringstream ost;
    ost << "PolyatomicElecLev::string() -> Implement me";
	
    return ost.str();
}

double
PolyatomicElecLev::
calculate_equilibrium_Q_total( double T )
{
    return calculate_Q_el( T ) * calculate_Q_vib( T, T );
}

double
PolyatomicElecLev::
calculate_Q_vib( double T_vib, double T_rot )
{
    /* Full summation expression for vibrational partition function */
    double QvibQrot = 0.0;
    
    cout << "PolyatomicElecLev::calculate_Q_vib()" << endl
         << "Fix me!" << endl;
    exit( NOT_IMPLEMENTED_ERROR );
    
    return QvibQrot;
}

double
PolyatomicElecLev::
calculate_and_store_Q_vib( double T_vib, double T_rot )
{
    /* Full summation expression for vibrational partition function */
    QvQr = 0.0;

    cout << "calculate_and_store_Q_vib::calculate_and_store_Q_vib()" << endl
         << "Fix me!" << endl;
    exit( NOT_IMPLEMENTED_ERROR );

    return QvQr;
}

PolyatomicRadiator::
PolyatomicRadiator( lua_State * L, string name )
 : Radiator(L, name)
{
    iTv = get_int( L, -1, "iTv" );
    
    iTr = get_int( L, -1, "iTr" );
    
    D = get_number( L, -1, "eta_D" );
    D *= RC_c * RC_h_SI;		// Convert cm**-1 -> J
    
    read_photoionization_data( L );

    // read_elev_data( L ); ???
}

PolyatomicRadiator::
~PolyatomicRadiator()
{
    for ( int ilev=0; ilev<nlevs; ++ilev )
    	delete elevs[ilev];
    
    // for ( size_t isys=0; isys<systems.size(); ++isys )
    //	delete systems[isys];
}

void
PolyatomicRadiator::
set_e_index( int iel )
{
    Radiator::set_e_index( iel );
    
    // for ( size_t isys=0; isys<systems.size(); ++isys ) {
    // 	systems[isys]->e_index = iel;
    // }
}


ElecLev *
PolyatomicRadiator::
get_elev_pointer( int ie )
{
    return elevs[ie];
}

double
PolyatomicRadiator::
get_D()
{
    return D;
}

void
PolyatomicRadiator::
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
PolyatomicRadiator::
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
PolyatomicRadiator::
initialise_mechanisms( Gas_data &Q )
{
    // 0. Pre-calculate n_hvy, n_elecs, mw_av
    // NOTE: we are assuming Q.p includes the electron contribution, and that 
    //       Q.p_e is present and correct
    //       [if a 1T model is being used p_e will be zero!]
//    double N_elecs = 0.0;
//    if ( e_index >= 0 )
//    	N_elecs = Q.p_e * RC_Na/(RC_R_u*Q.T[iTe]) * 1.0e-6;
//    double N_hvy = (Q.p - Q.p_e ) * RC_Na/(RC_R_u*Q.T[iT]) * 1.0e-6;
//    double mw_av = ( Q.rho * RC_Na ) / ( ( N_hvy + N_elecs ) * 1.0e6 ) * 1.0e3;
    
    // for (size_t isys=0; isys<systems.size(); ++isys) {
	// systems[isys]->initialize( Q.T[iT], Q.T[iTe], Q.p, N_hvy, N_elecs, mw_av );
    // }
    
    return;
}

double
PolyatomicRadiator::
calculate_unresolved_emission_coefficient( Gas_data &Q )
{
    double j_ul = 0.0;
    
    /* Pre-initialised, just sum the contributions */

    // for (size_t isys=0; isys<systems.size(); ++isys)
	// j_ul += systems[isys]->calculate_j_ul(Q.T[iTv],Q.T[iTr]);

    return j_ul;
}

double
PolyatomicRadiator::
calculate_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u )
{
    double j_ul = 0.0;
    
    /* Pre-initialised, just sum the contributions */

    // for (size_t isys=0; isys<systems.size(); ++isys)
	// j_ul += systems[isys]->calculate_OV_j_ul(Q.T[iTv],Q.T[iTr],wavel_switch,Lambda_l,Lambda_u);

    return j_ul;
}

void
PolyatomicRadiator::
spectral_distribution( vector<double> &nus )
{
    return;
}

void
PolyatomicRadiator::
calculate_spectrum( Gas_data &Q, CoeffSpectra &X )
{
    // Calculate the unresolved total emission to measure the importance of each band
    // double j_av = 0.0;
    // for (size_t isys=0; isys<systems.size(); ++isys) {
    	// j_av += systems[isys]->calculate_j_ul(Q.T[iTv],Q.T[iTr]) /
    	// 		double(systems[isys]->lRe_dim * systems[isys]->uRe_dim);
    // }
    
    // Loop over the systems and add the contributions
    // for (size_t isys=0; isys<systems.size(); ++isys) {
   	// systems[isys]->calculate_spectrum(X,Q.T[iTv],Q.T[iTr],j_av);
    // }
    
    return;
}

string
PolyatomicRadiator::
line_width_string( Gas_data &Q )
{
    // 0. Pre-calculate n_hvy, n_elecs, mw_av
    // NOTE: we are assuming Q.p includes the electron contribution, and that 
    //       Q.p_e is present and correct
    //       [if a 1T model is being used p_e will be zero!]
//    double N_elecs = 0.0;
//    if ( e_index >= 0 )
//    	N_elecs = Q.p_e * RC_Na/(RC_R_u*Q.T[iTe]) * 1.0e-6;
//    double N_hvy = (Q.p - Q.p_e ) * RC_Na/(RC_R_u*Q.T[iT]) * 1.0e-6;
//    double mw_av = ( Q.rho * RC_Na ) / ( ( N_hvy + N_elecs ) * 1.0e6 ) * 1.0e3;
    
    // 1. Loop through systems, call line_width_string function
    string lws = "";
    // for (size_t isys=0; isys<systems.size(); ++isys) {
	// lws += systems[isys]->line_width_string( Q.T[iT], Q.T[iTe], Q.p, N_hvy, N_elecs, mw_av );
    // }
    
    return lws;
}

/************************** BoltzLinearPolyatomicRadiator **************************/

BoltzLinearPolyatomicRadiator::
BoltzLinearPolyatomicRadiator( lua_State * L, std::string name )
 : PolyatomicRadiator(L, name) {}
 
BoltzLinearPolyatomicRadiator::
~BoltzLinearPolyatomicRadiator() {}

void
BoltzLinearPolyatomicRadiator::
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
