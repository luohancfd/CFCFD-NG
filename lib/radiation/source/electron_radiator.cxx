/** \file electron_radiator.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: Improved port from old lib/radiation
 *
 *  \brief Definitions for electron radiator classes
 *
 **/

#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>

#include "electron_radiator.hh"
#include "radiation_constants.hh"
#include "photaura.hh"

#include "../../gas/models/gas_data.hh"
#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

using namespace std;

ContinuumMechanism::ContinuumMechanism()
{}

ContinuumMechanism::~ContinuumMechanism() {}

FreeFreeHydrogenic::FreeFreeHydrogenic( Radiator * R_ion )
: Ri( R_ion )
{}

FreeFreeHydrogenic::~FreeFreeHydrogenic() {}

// ff_constA = ( 8.0 / 3.0 ) * sqrt( 2.0 * M_PI / ( 3.0 * RC_m * RC_k ) ) *
//                pow( RC_e, 6 ) / ( RC_m * pow( RC_c, 3 ) );
const double ff_constA = 5.4429513110835207e-39;
// ff_constB = ( 4.0 / 3.0 ) * sqrt( 2.0 * M_PI / ( 3.0 * RC_m * RC_k ) ) *
//                  pow( RC_e, 6 ) / ( RC_h * RC_c * RC_m );
const double ff_constB = 3.6913858809e+08;

void
FreeFreeHydrogenic::calculate_spectrum( Gas_data &Q, CoeffSpectra &X, double n_elecs )
{
    /* 1. Initialise frequency independent parameters */
    int nnu = int ( X.nu.size() ); 
    
    double T = Q.T[Ri->iTe];
    double n_ions = Q.massf[Ri->isp] * Q.rho / Ri->m_w * RC_Na * 1.0e-6;
    
    double j_ff_tmpA = n_elecs * ff_constA;
    double j_ff_tmpB = double(Ri->Z * Ri->Z) / sqrt(T) * n_ions;
    double kappa_ff_tmpA = n_elecs * ff_constB;
    double kappa_ff_tmpB = j_ff_tmpB;
    
    /* 2. Loop over frequency */
    for ( int inu=0; inu<nnu; ++inu ) {
	double nu = X.nu[inu];
	double j_ff_nu = j_ff_tmpA * j_ff_tmpB * exp ( - RC_h_SI * nu / (RC_k_SI * T) ) * 1.0e-1;     // erg / cm**-3 - Hz -> J / m**-3 - Hz
	double kappa_ff_nu = kappa_ff_tmpA * kappa_ff_tmpB / ( nu * nu * nu ) * 1.0e2;		  	// 1 / cm -> 1 / m

	if ( !std::isnan(j_ff_nu) && !std::isinf(j_ff_nu) && 
	     !std::isnan(kappa_ff_nu) && !std::isinf(kappa_ff_nu) ) {
	    X.j_nu[inu] += j_ff_nu;
	    X.kappa_nu[inu] += kappa_ff_nu;
	}
	else {
	    cout << "FreeFreeHydrogenic::calculate_spectrum()" << endl
	         << "j_ff_nu = " << j_ff_nu << ", kappa_ff_nu = " << kappa_ff_nu << endl;
	    exit( FAILURE );
	}
    }
    
    return;
}

void
FreeFreeHydrogenic::initialise_mechanisms( Gas_data &Q )
{
    return;
}

BoundFreeHydrogenic::BoundFreeHydrogenic( Radiator * R_ion, Radiator * R_neutral )
: Ri( R_ion ), Rn( R_neutral )
{}

BoundFreeHydrogenic::~BoundFreeHydrogenic() {}

// bf_constA = 2.0 * RC_h / pow( RC_c, 2 )
const double bf_constA = 1.4745007665631073e-47;

// bf_constB = pow( RC_h * RC_h / ( 2.0 * M_PI * RC_m * RC_k ), 1.5 )
const double bf_constB = 4.1412964647957627e-16;

void
BoundFreeHydrogenic::calculate_spectrum( Gas_data &Q, CoeffSpectra &X, double n_elecs )
{
    /* 1. Initialise frequency independent parameters */
    int nnu = int ( X.nu.size() );

    double T = Q.T[Ri->iTe];
    double n_ions = Q.massf[Ri->isp] * Q.rho / Ri->m_w * RC_Na * 1.0e-6;
    double Q_ion = Ri->Q_el;

    // combine all frequency independent params
    double bf_constC = n_ions * n_elecs * bf_constA / ( 2.0 * Q_ion ) * bf_constB * pow( T, -1.5 );

    /* 2. Loop over neutral levels */
    for ( int ilev=0; ilev<Rn->nlevs; ++ilev ) {
        double N_i = Rn->get_elev_pointer(ilev)->N * 1.0e-6;
        double E_i = Rn->get_elev_pointer(ilev)->E;
        // Determine minimum frequency
        double nu_min = ( Rn->I - E_i ) / RC_h_SI;
        int inu_min = get_nu_index( X.nu, nu_min, X.adaptive ) + 1;
        /* 2a. Loop over frequency and add contributions */
        for ( int inu=inu_min; inu<nnu; ++inu ) {
            double nu = X.nu[inu];
            double E_min = Rn->I - RC_h_SI * nu;
            double sigma_bf = Rn->get_elev_pointer(ilev)->calculate_sigma_bf( nu );
            // cout << "ilev = " << ilev << ", inu = " << inu << ", sigma_bf = " << sigma_bf << endl;
            X.j_nu[inu] += bf_constC * pow( nu, 3 ) * Rn->get_elev_pointer(ilev)->g * sigma_bf * \
                        exp( ( E_min - E_i ) / ( RC_k_SI * T ) ) * 1.0e-1;              // convert erg - s / cm**3 - sr -> W / m**3 - sr
            X.kappa_nu[inu] += sigma_bf * N_i * 1.0e2;  // convert 1 / cm -> 1 / m
        }
    }

    return;
}

void
BoundFreeHydrogenic::initialise_mechanisms( Gas_data &Q )
{
    return;
}

ElectronRadiator::
ElectronRadiator( lua_State * L, string name )
 : Radiator(L, name)
{
    // Set internal partition function as it is a constant
    Q_int = 2.0;
    
    // Create elev
    elev = new ElecLev(0,0.0,2);
    elev->PICS_model = new NoPICSModel();
    
    // Get the systems list
    lua_getfield(L, -1, "systems_list" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ElectronRadiator::ElectronRadiator():\n";
	ost << "Error in the declaration of systems_list: a table is expected.\n";
	input_error(ost);
    }
    
    int nsys = lua_objlen(L, -1);
        
    for ( int isys = 0; isys < nsys; ++isys ) {
	lua_rawgeti(L, -1, isys+1); // A Lua list is offset one from the C++ vector index
	const char* sys = luaL_checkstring(L, -1);
	systems_list.push_back(string(sys));
	lua_pop(L, 1);
    }
    
    lua_pop(L,1);	// pop systems list
}

ElectronRadiator::
~ElectronRadiator()
{
    delete elev;
    
    for ( size_t icm=0; icm<continuum_mechanisms.size(); ++icm )
    	delete continuum_mechanisms[icm];
}

void
ElectronRadiator::
create_continuum_mechanisms( vector<Radiator*> &radiators )
{
    // Find neutral-ion pairs and create continuum mechanisms
    for ( size_t irad = 0; irad < radiators.size(); ++irad ) {
    	Radiator * iR = radiators[irad];
    	ostringstream ion_name;
    	ion_name << iR->name << "_plus";
    	for ( size_t jrad = 0; jrad < radiators.size(); ++jrad ) {
    	    if ( jrad==irad ) continue;
    	    Radiator * jR = radiators[jrad];
    	    if ( jR->name == ion_name.str() ) {
    	    	// jR is the ion of iR
    	    	for ( int isys = 0; isys < int(systems_list.size()); ++isys ) {
    	    	    // NOTE: assuming no detailed photoionization cross section data is available
    	    	    if ( systems_list[isys]=="free-free" ) {
    	    	    	 cout << " - Created FreeFreeHydrogenic mechanism: " << jR->name << endl;
    	    	    	continuum_mechanisms.push_back( new FreeFreeHydrogenic( jR ) );
    	    	    }
    	    	    else if ( systems_list[isys]=="bound-free" ) {
    	    	    	cout << " - Creating BoundFreeHydrogenic mechanism: " << jR->name << " + e_minus <-> " << iR->name << endl;
    	    	    	continuum_mechanisms.push_back( new BoundFreeHydrogenic( jR, iR ) );
    	    	    	// set the bf_ion_flag to true for jR
    	    	    	jR->bf_ion_flag = true;
    	    	    }
    	        }
    	    }
    	}
    }
    
    return;
}

ElecLev *
ElectronRadiator::
get_elev_pointer( int ie )
{
    if ( ie != 0 ) {
	cout << "ElectronRadiator::get_elev_pointer()" << endl
	     << "Electrons only have a single level with index 0!" << endl
	     << "Bailing out!" << endl;
	exit( FAILURE );
    }
    
    return elev;
}

void
ElectronRadiator::
calculate_Q_int( Gas_data &Q )
{
    // Do nothing (constant, set in constructor)
}

void
ElectronRadiator::
calculate_n_e( Gas_data &Q )
{
    elev->set_N( Q.massf[isp] * Q.rho / RC_m_SI );	// convert kg/m**3 -> electrons/m**3
    return;
}

double
ElectronRadiator::
calculate_total_equil_partition_function( double T )
{
    // NOTE: Formation energy is zero and makes no contribution
    
    // 1. Translational contribution
    double Q_tr = pow( 2.0 * M_PI * m_w / RC_Na * RC_k_SI * T / RC_h_SI / RC_h_SI, 1.5 );
    
    // 3. Electronic contribution (set as a constant)
    double Q_elec = Q_int;
    
    // 4. Return the product of the modal contributions
    return Q_tr * Q_elec;
}

void
ElectronRadiator::
spectral_distribution( std::vector<double> &nus )
{
    return;
}

void
ElectronRadiator::
calculate_spectrum( Gas_data &Q, CoeffSpectra &X )
{
    double n_elecs = elev->get_N() * 1.0e-6;
    
    // loop over continuum mechanisms
    for ( size_t icm=0; icm<continuum_mechanisms.size(); ++icm )
    	continuum_mechanisms[icm]->calculate_spectrum(Q,X,n_elecs);
    	
    return;
}

double
ElectronRadiator::
calculate_unresolved_emission_coefficient( Gas_data &Q )
{
    return 0.0;
}

void
ElectronRadiator::
initialise_mechanisms( Gas_data &Q )
{
    // loop over continuum mechanisms
    // for ( size_t icm=0; icm<continuum_mechanisms.size(); ++icm )
    	// continuum_mechanisms[icm]->initialise_mechanisms(Q);
    
    return;
}

string
ElectronRadiator::
line_width_string( Gas_data &Q )
{
    return "";
}

