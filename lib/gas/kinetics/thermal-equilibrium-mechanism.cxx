// Author: Daniel F. Potter
// Date: 10-May-2010

#include <iostream>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "thermal-equilibrium-mechanism.hh"

using namespace std;

Thermal_equilibrium_mechanism::
Thermal_equilibrium_mechanism( lua_State * L )
{
    string species_name = get_string(L, -1, "species");
    Chemical_species * X = get_library_species_pointer_from_name( species_name );
    isp_ = X->get_isp();
    
    string mode_type = get_string(L, -1, "mode_type");
    
    mode_ = 0;
    for ( int imode=0; imode<X->get_n_modes(); ++imode ) {
        if ( X->get_mode_pointer(imode)->get_type()==mode_type ) {
            mode_ = X->get_mode_pointer(imode);
            break;
        }
    }
    
    if ( mode_==0 ) {
    	ostringstream oss;
    	oss << "Thermal_equilibrium_mechanism::Thermal_equilibrium_mechanism" << endl
    	    << "mode_type: " << mode_type << " was not found for species: " << species_name << endl;
    	input_error( oss );
    }
    
    iT_ = mode_->get_iT();
    iT_eq_ = get_int( L, -1, "iT_eq" );
}


Thermal_equilibrium_mechanism::
~Thermal_equilibrium_mechanism()
{}

void
Thermal_equilibrium_mechanism::
s_apply( Gas_data &Q )
{
    // 1. Eval the old energy
    double e_old = Q.massf[isp_] * mode_->eval_energy_from_T( Q.T[iT_] );
    
    // 2. Eval the new energy
    double e_new = Q.massf[isp_] * mode_->eval_energy_from_T( Q.T[iT_eq_] );
    
    // 3. Update the bulk energy modes
    if ( iT_!=0 ) Q.e[iT_] += e_new - e_old;
    if ( iT_eq_!=0 ) Q.e[iT_eq_] += e_new - e_old;
    
    return;
}

