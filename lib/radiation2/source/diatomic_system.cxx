/** \file diatomic_system.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 21-Aug-07
 *  \brief Deifinitions for Diatomic_system class; simply a place-holder for transition details.
 *
 **/

#include <stdlib.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <math.h>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"
#include "../../gas/models/gas_data.hh"
#include "diatomic_system.hh"
#include "radiation_constants.hh"
#include "photaura.hh"
#include "diatomic_radiator.hh"

using namespace std;

DiatomicSystem::
DiatomicSystem( lua_State * L, string name, string band_method, 
    		DiatomicElecLev * elev_u, DiatomicElecLev * elev_l, 
    		int iT, int iTe, int iTv, int iTr, double m_w, double I  )
 : name( name ), elev_u( elev_u ), elev_l( elev_l), iT( iT ), iTe( iTe ), iTv( iTv ), iTr( iTr )
{
    // Set the transition type flag
    // NOTE: perhaps use to create band types
    transition_type = get_diatomic_transition_type(elev_u,elev_l);
    if ( ECHO_RAD_INPUT > 0 )
    	cout << name << " has transition type: " << transition_type << endl;
    
    // For the collision width
    double sigma_nm = get_positive_number( L, -1, "sigma_nm" );
    
    // Level indices (pointers set above)
    ie_l = get_int( L, -1, "ie_l" );
    ie_u = get_int( L, -1, "ie_u" );
    
    // The band data table
    lua_getfield( L, -1, "band_data" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "DiatomicSystem::DiatomicSystem():\n";
	ost << "Error getting field: band_data - a table is expected.\n";
	input_error(ost);
    }
    
    string format = get_string( L, -1, "format" );
    
    uRe_dim = get_int( L, -1, "uRe_dim" );
    lRe_dim = get_int( L, -1, "lRe_dim" );
    
    if ( ECHO_RAD_INPUT > 1 ) {
	cout << "ie_l = " << ie_l << ", ie_u = " << ie_u << endl
	     << "lRe_dim = " << lRe_dim << ", uRe_dim = " << uRe_dim << endl;
    }
    
    // Electronic transition moments 
    if ( uRe_dim<=0 || lRe_dim <=0 ) {
	cout << "No electronic transition moments present!" << endl;
    }
    else {
	/* Create valarray with space for all bands with data */
	for ( int iVu=0; iVu < uRe_dim; ++iVu ) {
	    for ( int iVl=0; iVl < lRe_dim; ++iVl ) {
		double Re_vib = get_Re_vib_from_file( L, iVu, iVl, format );
		if ( band_method=="linebyline" ) {
		    if ( transition_type == HUND_A )
		    	bands.push_back( new HundALBLDiatomicBand( iVu, iVl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm ) );
		    else if ( transition_type == SIGMA_TRIPLET )
		    	bands.push_back( new SigmaTripletLBLDiatomicBand( iVu, iVl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm ) );
		    else if ( transition_type == HUND_B_DOUBLET ) 
		    	bands.push_back( new HundBDoubletLBLDiatomicBand( iVu, iVl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm ) );
		    else if ( transition_type == HUND_AB_DOUBLET ) 
		    	bands.push_back( new HundABDoubletLBLDiatomicBand( iVu, iVl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm ) );
		    else if ( transition_type == HUND_AB_TRIPLET ) 
		    	bands.push_back( new HundABTripletLBLDiatomicBand( iVu, iVl, elev_u, elev_l, Re_vib, m_w, I, sigma_nm ) );
		}
		else if ( band_method=="smeared_band" )
		    bands.push_back( new SRBDiatomicBand( iVu, iVl, elev_u, elev_l, Re_vib ) );
#		if WRITE_LINES_TO_FILE
		// option to write lines to file if this is the 0-0 transition
		if ( iVu==0 && iVl==0 ) {
		    ostringstream oss;
		    oss << name << "_Vu" << iVu << "-Vl" << iVl << "-bands.txt";
		    bands.back()->write_lines_to_file( oss.str() );
		}
#               endif
	    }
	}
#       if DEBUG_RAD > 0
	// System dimensions check
        if ( uRe_dim > ( elev_u->get_V_max() + 1 ) || lRe_dim > ( elev_l->get_V_max() + 1 ) ) {
            cout << "DiatomicSystem::DiatomicSystem()" << endl
                 << "WARNING: System: " << name << " has an ETM matrix of size " << uRe_dim << " x " << lRe_dim << endl
                 << "but elev_u->V_max = " << elev_u->get_V_max() << " and elev_l->V_max = " << elev_l->get_V_max() << endl
                 << "The calculation is possible but the partition functions will not be consistent with the number bands in the spectra." << endl;
        }
#       endif

	// Check vector size
	if ( bands.size() != size_t(uRe_dim*lRe_dim) ) {
	    ostringstream ost;
	    ost << "DiatomicSystem::DiatomicSystem():\n";
	    ost << "Number of bands does not match given dimensions.\n";
	    input_error(ost);
	}
    }
    
    double nu_00 = ( elev_u->E + elev_u->calculate_E_vib(0) - elev_l->E - elev_l->calculate_E_vib(0) ) / RC_h_SI;
    cout << name << " lambda_00 = " << nu2lambda(nu_00) << endl;
    
    lua_pop(L,1);	// pop band_data
}

DiatomicSystem::
~DiatomicSystem()
{
    for ( size_t ib = 0; ib < bands.size(); ib++ )
	delete bands[ib];
}

double
DiatomicSystem::
get_Re_vib_from_file( lua_State * L, int iVu, int iVl, string format )
{
    ostringstream Vu_oss;
    Vu_oss << "Vu_" << iVu;
    lua_getfield(L, -1, Vu_oss.str().c_str());
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "DiatomicSystem::get_Re_vib_from_file()\n";
	ost << "Error locating " << Vu_oss.str() << " table" << endl;
	input_error(ost);
    }
    // Check number of Re_vib's
    if ( lua_objlen(L,-1) != size_t(lRe_dim) ) {
	ostringstream ost;
	ost << "DiatomicSystem::get_Re_vib_from_file()\n";
	ost << "Dimension mismatch (Vl's) for Re_vib table." << endl;
	input_error(ost);
    }
    // Get requested Re_vib element
    lua_rawgeti(L, -1, iVl+1);
    
    // Apply format specific conversions
    double Re_vib = 0.0;
    if ( format=="transition moment" ) 
    	Re_vib = luaL_checknumber(L, -1);
    else if ( format=="Einstein coefficient" ) {
    	double A = luaL_checknumber(L, -1);
    	Re_vib = convert_VEC_to_ETM( A, iVu, iVl );
    }
    else if ( format=="absorption oscillator strength" ) {
    	double f = luaL_checknumber(L, -1);
    	Re_vib = convert_AOS_to_ETM( f, iVu, iVl );
    }
    
    lua_pop(L,1);	// pop Re_vib element
    
    lua_pop(L,1);	// pop Vu
    
    return Re_vib;
}

DiatomicBand *
DiatomicSystem::
band_pointer( int Vu, int Vl )
{
    /* Return pointer to specified band (unravels 1D valarray) */
    if ( Vu >= uRe_dim || Vl >= lRe_dim ) {
	cout << "Requested " << name << " band [" << Vu << ", " << Vl
	     << "] is outside valarray of dimensions [" << uRe_dim << " x "
	     << lRe_dim  << "].  Bailing out!" 
	     << endl;
	exit( BAD_INPUT_ERROR );
    }

    int index = Vu * lRe_dim + Vl;

    return bands[index];
}

void
DiatomicSystem::
initialize( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    // Loop over bands and call initialize function
    for ( int Vu=0; Vu<uRe_dim; ++Vu ) {
    	for ( int Vl=0; Vl<lRe_dim; ++Vl ) {
    	    if ( band_pointer(Vu,Vl)->reverse_band ) continue;
    	    band_pointer(Vu,Vl)->initialize(T,Te,p,N_hvy,N_elecs,mw_av);
    	}
    }
    
    return;
}

double
DiatomicSystem::
calculate_j_ul( double T_vib, double T_rot )
{
    double j_ul = 0.0;
    
    // Loop over bands, add contributions
    for ( int Vu=0; Vu<uRe_dim; ++Vu ) {
    	double N_u = elev_u->calculate_N_vib( T_vib, T_rot, Vu );
    	for ( int Vl=0; Vl<lRe_dim; ++Vl ) {
    	    if ( band_pointer(Vu,Vl)->reverse_band ) continue;
    	    j_ul += band_pointer( Vu, Vl )->calculate_j_ul( N_u );
    	}
    }
    
    return j_ul;
}

double
DiatomicSystem::
calculate_OV_j_ul( double T_vib, double T_rot, double wavel_switch, double Lambda_l, double Lambda_u )
{
    double j_ul = 0.0;
    
    // Loop over bands, add contributions
    for ( int Vu=0; Vu<uRe_dim; ++Vu ) {
    	double N_u = elev_u->calculate_N_vib( T_vib, T_rot, Vu );
    	for ( int Vl=0; Vl<lRe_dim; ++Vl ) {
    	    if ( band_pointer(Vu,Vl)->reverse_band ) continue;
    	    if ( nu2lambda(band_pointer(Vu,Vl)->nu_00) <= wavel_switch ) 
    	    	j_ul += Lambda_l * band_pointer( Vu, Vl )->calculate_j_ul( N_u );
    	    else
    	    	j_ul += Lambda_u * band_pointer( Vu, Vl )->calculate_j_ul( N_u );
    	}
    }
    
    return j_ul;
}

void
DiatomicSystem::
calculate_spectrum( CoeffSpectra &X, double T_vib, double T_rot, double j_av )
{
    // Loop over bands, add contributions
    for ( int Vu=0; Vu<uRe_dim; ++Vu ) {
    	if ( Vu > elev_u->get_V_max() ) break;
    	double N_u = elev_u->calculate_N_vib( T_vib, T_rot, Vu );
    	for ( int Vl=0; Vl<lRe_dim; ++Vl ) {
    	    if ( Vl > elev_l->get_V_max() ) break;
    	    if ( band_pointer(Vu,Vl)->reverse_band ) continue;
    	    // NOTE: - This gives an ~ 3% reduction in computation time
    	    //         but will miss bands on the edges of spectra
#           if TEST_BAND_NU_AV
    	    if ( band_pointer(Vu,Vl)->nu_av < X.nu[0] || band_pointer(Vu,Vl)->nu_av > X.nu.back() ) continue;
#           endif
	    if ( band_pointer( Vu, Vl )->calculate_j_ul( N_u ) / j_av < F_BAND_LIMIT ) continue;
    	    band_pointer(Vu,Vl)->calculate_spectrum(X,T_vib,T_rot);
    	}
    }
    
    return;
}

std::string
DiatomicSystem::
line_width_string( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av )
{
    // Loop over bands and add contributions to the line-width string
    string lws = "";
    for ( int Vu=0; Vu<uRe_dim; ++Vu ) {
    	for ( int Vl=0; Vl<lRe_dim; ++Vl ) {
    	    if ( band_pointer(Vu,Vl)->reverse_band ) continue;
    	    lws += band_pointer(Vu,Vl)->line_width_string(T,Te,p,N_hvy,N_elecs,mw_av);
    	}
    }
    
    return lws;
}

double
DiatomicSystem::
convert_VEC_to_ETM( double A, int iVu, int iVl )
{
    double nu_av = calculate_band_average_frequency( iVu, iVl );
    double tmpA = 64.0 * pow( M_PI, 4 ) * pow( nu_av, 3 ) / ( 3.0 * RC_h * pow( RC_c, 3 ) );
    double tmpB = pow( RC_e*RC_a0, 2 );
    double tmpC = double( 2 - elev_u->eval_kronecker_delta() * elev_l->eval_kronecker_delta() ) \
                / double( 2 - elev_u->eval_kronecker_delta() );
    
    return sqrt( A / tmpA / tmpB / tmpC );
}

double
DiatomicSystem::
convert_AOS_to_ETM( double f, int iVu, int iVl )
{
    double nu_av = calculate_band_average_frequency( iVu, iVl );
    double tmpA = 8.0 * pow( M_PI, 2 ) * RC_m * nu_av / ( 3.0 * RC_h * pow( RC_e, 2 ) );
    double tmpB = pow( RC_e*RC_a0, 2 );
    double tmpC = double( 2 - elev_u->eval_kronecker_delta() * elev_l->eval_kronecker_delta() ) \
                / double( 2 - elev_u->eval_kronecker_delta() );
    
    return sqrt( f / tmpA / tmpB / tmpC );
}

double
DiatomicSystem::
calculate_band_average_frequency( int iVu, int iVl )
{
    // NOTE: this function should replicate the behaviour of 
    //       DiatomicBand::calculate_average_frequency()
    double E_u = elev_u->calculate_E_vib(iVu) + elev_u->E; 
    double E_l = elev_l->calculate_E_vib(iVl) + elev_l->E;
    
    double nu_av = 0.0;
#   if BAND_AVERAGE_FREQUENCY_METHOD==0    
    // 3. Pull out lambda and J_max values from the electronic levels
    int tJu_max = 2*elev_u->get_J_max(iVu);
    int tJl_max = 2*elev_l->get_J_max(iVl);
    
    int lambda_u = elev_u->get_lambda();
    int lambda_l = elev_l->get_lambda();

    // 4. Calculate characteristic frequencies
    // 4a. 00 transition frequency
    double nu_00 = ( E_u - E_l ) / RC_h_SI;
    
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
	    double nu_cl = nu_00 + ( elev_u->calculate_E_rot_singlet(iVu,tJu) - elev_l->calculate_E_rot_singlet(iVl,tJl) ) / RC_h_SI;
	    nu_acc += fabs(nu_cl);
	    ++count;
	}
    }
    nu_av = nu_acc / double( count );
#   elif BAND_AVERAGE_FREQUENCY_METHOD==1
    double nu_00 = ( E_u - E_l ) / RC_h_SI;
    double Bv_u  = elev_u->calculate_B_v( iVu );
    double Bv_l  = elev_l->calculate_B_v( iVl );
    int J_max_l = (int) elev_l->get_J_max(iVl);
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

double
DiatomicSystem::
calculate_transition_probability( double Tv )
{
    double numer = 0.0;
    double denom = 0.0;
    for ( int Vu=0; Vu<uRe_dim; ++Vu ) {
    	double A_v = 0.0;
    	for ( int Vl=0; Vl<lRe_dim; ++Vl ) {
    	    if ( band_pointer(Vu,Vl)->reverse_band ) continue;
    	    A_v += band_pointer(Vu,Vl)->A_ul_av;
    	}
	double tmp = exp( - elev_u->calculate_E_vib( Vu ) / RC_k_SI / Tv );
	numer += A_v * tmp;
	denom += tmp;
    }
    
    return numer / denom;
}

DiatomicSystem * create_new_diatomic_system( lua_State * L, std::string name, 
    DiatomicElecLev * elev_u, DiatomicElecLev * elev_l,
    int iT, int iTe, int iTv, int iTr, double m_w, double I )
{
    // the system should now be at the top of the stack
    
    string band_method = get_string(L, -1, "band_method");
    
    if ( band_method == "linebyline" )  {
	cout << "  * Creating '" << name << "' as a new DiatomicSystem with LBL bands" << endl;
    }
    else if ( band_method == "smeared_band" )  {
	cout << "  * Creating '" << name << "' as a new DiatomicSystem with SRB bands" << endl;
    }
    else {
	ostringstream ost;
	ost << "The specified band method: " << band_method << endl
	    << "is not available or no yet implemented.\n" << endl
	    << "Bailing Out!\n";
	input_error(ost);
    }
    
    return new DiatomicSystem(L, name, band_method, elev_u, elev_l, iT, iTe, iTv, iTr, m_w, I);
}

