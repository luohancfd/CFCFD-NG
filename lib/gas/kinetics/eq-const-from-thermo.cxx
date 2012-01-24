// Author: Rowan J. Gollan
// Date: 21-Oct-2008
// Place: Poquoson, Virginia, USA
//

#include <iostream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "../models/physical_constants.hh"
#include "../models/chemical-species-library.hh"

#include "eq-const-from-thermo.hh"

using namespace std;

Eq_const_from_thermo::
Eq_const_from_thermo(map<int, int> &nu, Gas_model &g, int iT)
    : Equilibrium_constant( iT ), nu_(nu), g_(g)
{
    Q_ = new Gas_data(&g_);
}

Eq_const_from_thermo::
~Eq_const_from_thermo()
{
    delete Q_;
}

double
Eq_const_from_thermo::
s_eval(const Gas_data &Q)
{
    // initialise the temporary gas-data structure with the appropriate rate
    // controlling temperature
    Q_->copy_values_from(Q);
    for ( size_t imode=0; imode<Q.T.size(); ++imode ) {
    	Q_->T[imode] = Q.T[iT_];
    }
    double dG = 0.0;
    int nu_sum = 0;
    map<int, int>::const_iterator it;
    for ( it = nu_.begin(); it != nu_.end(); ++it ) {
	// We want to use G in J/mol, hence conversion based on molecular weight.
	dG += it->second * g_.Gibbs_free_energy(*Q_, it->first)*(g_.molecular_weight(it->first));
	nu_sum += it->second;
	// printf("G=%16.15e\n",g_.Gibbs_free_energy(Q_, it->first)*(g_.molecular_weight(it->first)));
    }
    double K_p = exp(-dG/(PC_R_u*Q_->T[0]));
    double K_c = K_p*pow(PC_P_atm/(PC_R_u*Q_->T[0]), nu_sum);

    // printf("K_c=%16.15e\n", K_c);
    return K_c;
}

Equilibrium_constant* create_Eq_const_from_thermo(lua_State *L,
						  std::map<int, int> &nu,
						  Gas_model &g)
{
    int iT = get_int( L, -1, "iT" );
    if ( iT<0 ) {
    	// The temperature of a specific species mode is being requested
    	string species = get_string(L,-1,"species");
    	string mode = get_string(L,-1,"mode");
    	iT = get_library_species_pointer_from_name(species)->get_mode_pointer_from_type(mode)->get_iT();
    }
    
    if ( iT<0 || iT >= g.get_number_of_modes() ) {
	ostringstream ost;
	ost << "create_Eq_const_from_thermo():\n";
	ost << "Error in specification of rate controlling temperature index.\n";
	ost << "iT = " << iT << ", nmodes = " << g.get_number_of_modes() << endl;
	input_error(ost);
    }
    
    return new Eq_const_from_thermo(nu, g, iT);
}
