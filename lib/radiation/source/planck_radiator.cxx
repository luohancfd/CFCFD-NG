/** \file planck_radiator.cxx
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: Improved port from old lib/radiation
 *
 *  \brief Definitions for planck radiator classes
 *
 **/

#include <stdlib.h>
#include <iostream>
 
#include "../../gas/models/gas_data.hh"
#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "planck_radiator.hh"
#include "spectral_model.hh"

using namespace std;

PlanckRadiator::
PlanckRadiator( lua_State * L, string name )
 : Radiator(L, name)
{
    /* Constant Absorption  */
    // Default to kappa = 1.0 1/m
    kappa_const = get_positive_number( L, -1, "kappa_const" );
    if ( ECHO_RAD_INPUT > 0 )
	cout << "kappa_const = " << kappa_const << " 1/m" << endl;
}

PlanckRadiator::
~PlanckRadiator() {}

void
PlanckRadiator::
set_e_index( int iel )
{
    Radiator::set_e_index( iel );
}

ElecLev *
PlanckRadiator::
get_elev_pointer( int ie )
{
    cout << "PlanckRadiator::get_elev_pointer()" << endl
	 << "PlanckRadiator's don't have any levels!" << endl
	 << "Bailing out!" << endl;
    exit( FAILURE );
    
    return 0;
}

void
PlanckRadiator::
calculate_Q_int( Gas_data &Q )
{}

void
PlanckRadiator::
calculate_n_e( Gas_data &Q )
{}

double
PlanckRadiator::
calculate_total_equil_partition_function( double T )
{
    return 0.0;
}

void
PlanckRadiator::
initialise_mechanisms( Gas_data &Q )
{
    return;
}

double
PlanckRadiator::
calculate_unresolved_emission_coefficient( Gas_data &Q )
{
    return 0.0;
}

void
PlanckRadiator::
spectral_distribution( std::vector<double> &nus )
{}

void
PlanckRadiator::
calculate_spectrum( Gas_data &Q, CoeffSpectra &X )
{
    double T = Q.T[iTe];
    
    for ( size_t inu=0; inu<X.nu.size(); ++inu ) {
	double nu = X.nu[inu];
	X.j_nu[inu] += kappa_const * planck_intensity( nu, T );
	X.kappa_nu[inu] += kappa_const;
    }
    return;
}

string
PlanckRadiator::
line_width_string( Gas_data &Q )
{
    return "";
}
