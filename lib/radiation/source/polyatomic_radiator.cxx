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
// #include "polyatomic_system.hh"
#include "radiation_constants.hh"
#include "photaura.hh"

using namespace std;

PolyatomicVibState::PolyatomicVibState( string type, int iV, lua_State * L )
  : type( type ), iV( iV )
{
    // The first entry is the state label, a string
    lua_rawgeti(L, -1, 1);
    label = luaL_checkstring(L, -1);
    lua_pop(L, 1 );
    // Now the rest which are numbers
    vector<double> state_data;
    for ( size_t i=1; i<lua_objlen(L, -1); ++i ) {
        lua_rawgeti(L, -1, i+1);
        state_data.push_back( luaL_checknumber(L, -1) );
        lua_pop(L, 1 );
    }

    g = (int) state_data[0];
    p_v = (int) state_data[1];
    G_v = state_data[2] * RC_c * RC_h_SI;                // cm-1 -> J
    A_v = state_data[3] * RC_c * RC_h_SI;
    B_v = state_data[4] * RC_c * RC_h_SI;
    C_v = state_data[5] * RC_c * RC_h_SI;
    D_v = state_data[6] * 1.0e-7 * RC_c * RC_h_SI;
    H_v = state_data[7] * 1.0e-13 * RC_c * RC_h_SI;
    J_max = (int) state_data[8];
}

double
PolyatomicVibState::
calculate_Q_rot( double T )
{
    double Q_rot = 0.0;    
  
    for ( int iJ=0; iJ<J_max; ++iJ ) {
        for ( int iK=-iJ; iK<=iJ; ++iK ) {
	    double E_rot = this->calculate_E_rot(iJ,iK);
            Q_rot += double(2*iJ+1) * exp( - E_rot / RC_k_SI / T );
        }
    }

    return Q_rot;
}

SphericalTopPolyatomicVibState::SphericalTopPolyatomicVibState( int iV, lua_State * L )
  : PolyatomicVibState( "SphericalTop", iV, L )
{}

double
SphericalTopPolyatomicVibState::calculate_E_rot( int iJ, int iK )
{
    /*
     * Rotational state energy for a spherical top (no spin splitting)
     * Ref: Equation 2 in Rothman et al (1992) JQSRT v48 n5/6 pp 527-566
     */

    UNUSED_VARIABLE( iK );

    // Quantum numbers
    double J = double(iJ);

    double E_r = B_v*J*(J+1.0) - D_v*J*J*(J+1.0)*(J+1.0) \
                + H_v*J*J*J*(J+1.0)*(J+1.0)*(J+1.0);

    return E_r;
}

SymmetricalTopPolyatomicVibState::SymmetricalTopPolyatomicVibState( int iV, lua_State * L )
  : PolyatomicVibState( "SymmetricalTop", iV, L )
{}

double
SymmetricalTopPolyatomicVibState::calculate_E_rot( int iJ, int iK )
{
    /*
     * Rotational state energy for a symmetric top (no spin splitting)
     * Ref: Equation 43 in Capitelli (2005)
     */

    // Quantum numbers
    double J = double(iJ);
    double K = double(iK);

    double E_r = B_v*J*(J+1.0) + ( A_v - B_v )*K*K;

    // CHECKME: perhaps we should add these anharmonicity terms?
    // E_r += D_v*J*J*(J+1.0)*(J+1.0) + H_v*J*J*J*(J+1.0)*(J+1.0)*(J+1.0);

    return E_r;
}

AsymmetricTopPolyatomicVibState::AsymmetricTopPolyatomicVibState( int iV, lua_State * L )
  : PolyatomicVibState( "AsymmetricTop", iV, L )
{}

double
AsymmetricTopPolyatomicVibState::calculate_E_rot( int iJ, int iK )
{
    /*
     * Rotational state energy for an Asymmetric top (no spin splitting)
     * Ref: Equation 44 in Capitelli (2005)
     */

    // Quantum numbers
    double J = double(iJ);
    double K = double(iK);

    double E_r = 0.5*( B_v + C_v )*J*(J+1.0) + ( A_v - B_v )*K*K;

    // CHECKME: perhaps we should add these anharmonicity terms?
    // E_r += D_v*J*J*(J+1.0)*(J+1.0) + H_v*J*J*J*(J+1.0)*(J+1.0)*(J+1.0);

    return E_r;
}

PolyatomicVibState * create_new_polyatomic_vib_state( string type, int iV, lua_State * L )
{
    PolyatomicVibState * PVS = 0;

    if ( type.find("SphericalTop")!=string::npos )
        PVS = new SphericalTopPolyatomicVibState( iV, L );
    else if ( type.find("SymmetricalTop")!=string::npos )
        PVS = new SymmetricalTopPolyatomicVibState( iV, L );
    else if ( type.find("AsymmetricTop")!=string::npos )
        PVS = new AsymmetricTopPolyatomicVibState( iV, L );
    else {
        ostringstream ost;
        ost << "create_new_polyatomic_vib_state()\n";
        ost << "Symmetry type not known: " << type << endl;
        input_error(ost);
    }

    return PVS;
}

PolyatomicElecLev::PolyatomicElecLev( int ilev, double E, int g, lua_State * L )
 : ElecLev( ilev, E, g )
{
    // get the rest of the data from the lua_State
    int Nvib_states = get_int(L,-1,"Nvib_states");
    type = get_string(L,-1,"type");

    for ( int iV=0; iV<Nvib_states; ++iV) {
        ostringstream state_oss;
        state_oss << "vibstate_" << iV;
        lua_getfield(L, -1, state_oss.str().c_str());
        if ( !lua_istable(L, -1) ) {
            ostringstream ost;
            ost << "PolyatomicElecLev::PolyatomicElecLev()\n";
            ost << "Error locating " << state_oss.str()<< " table" << endl;
            input_error(ost);
        }
        vstates.push_back( create_new_polyatomic_vib_state(type,iV,L) );
        lua_pop(L,1);   // pop vibstate
    }
    UNUSED_VARIABLE(D);
}

PolyatomicElecLev::~PolyatomicElecLev()
{
    for ( size_t iV=0; iV<vstates.size(); ++iV )
        delete vstates[iV];
}

string
PolyatomicElecLev::
string()
{
    ostringstream ost;
    ost << setw(15) << E / ( RC_c * RC_h_SI )
        << setw(5)  << g;

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
    
    for (size_t iV=0; iV<vstates.size(); iV++) {
        // Vibrational energy
        double E_v = vstates[iV]->get_E_vib();
        // Vibrational statistical weight
        int p_v = vstates[iV]->get_p_vib();
        // Rotational partition function
        double Q_rot = vstates[iV]->calculate_Q_rot(T_rot);
        QvibQrot += Q_rot * p_v * exp( - E_v / ( RC_k_SI * T_vib ) );
    }
    
    return QvibQrot;
}

double
PolyatomicElecLev::
calculate_and_store_Q_vib( double T_vib, double T_rot )
{
    QvQr = calculate_Q_vib(T_vib,T_rot);

    return QvQr;
}

PolyatomicRadiator::
PolyatomicRadiator( lua_State * L, string name )
 : Radiator(L, name)
{

#   if DEBUG_RAD > 0
    cout << "PolyatomicRadiator::PolyatomicRadiator()" << endl
         << "Species name: " << name << endl;
#   endif

    iTv = get_int( L, -1, "iTv" );
    
    iTr = get_int( L, -1, "iTr" );
    
    D = get_number( L, -1, "eta_D" );
    D *= RC_c * RC_h_SI;		// Convert cm**-1 -> J
    
    read_elevel_data( L );

    read_photoionization_data( L );

#   if DEBUG_RAD > 0
    cout << "PolyatomicRadiator::PolyatomicRadiator()" << endl
         << "Done." << endl;
#   endif


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
read_elevel_data( lua_State * L )
{
    lua_getfield(L, -1, "level_data");
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "PolyatomicRadiator::read_elevel_data()\n";
        ost << "Error locating 'elec_levels' table" << endl;
        input_error(ost);
    }

    nlevs = get_positive_int(L, -1, "n_levels");

    for ( int ilev=0; ilev<nlevs; ++ilev ) {
        ostringstream lev_oss;
        lev_oss << "ilev_" << ilev;
        lua_getfield(L, -1, lev_oss.str().c_str());
        if ( !lua_istable(L, -1) ) {
            ostringstream ost;
            ost << "PolyatomicRadiator::read_elevel_data()\n";
            ost << "Error locating " << lev_oss.str() << " table" << endl;
            input_error(ost);
        }
        double E = get_number(L, -1, "theta" ) * RC_c * RC_h_SI;
        int g = get_int(L, -1, "g");
        elevs.push_back( new PolyatomicElecLev( ilev, E, g, L ) );
        lua_pop(L,1);   // pop ilev
    }

    lua_pop(L,1);       // pop elec_levels

    if ( ECHO_RAD_INPUT > 1 ) {
        cout << setw(10) << "[ilev]"
             << setw(15) << "[E_el (1/cm)]"
             << setw(5)  << "[g_el]" << endl;
        for ( int ilev=0; ilev<nlevs; ++ilev ) {
            cout << "ilev_" << ilev << " = " << elevs[ilev]->string() << endl;
        }
    }
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
 : PolyatomicRadiator(L, name)
{
#   if DEBUG_RAD > 0
    cout << "BoltzLinearPolyatomicRadiator::BoltzLinearPolyatomicRadiator()" << endl
         << "Species name: " << name << endl;
#   endif
}
 
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
