// Author: Daniel F. Potter
// Date: 20-Apr-2010
// Place: Dutton Park, Brisbane, QLD
//

#include <iostream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "../models/physical_constants.hh"
#include "../models/chemical-species-library.hh"

#include "eq-const-from-CEA-curves.hh"

using namespace std;

Eq_const_from_CEA_curves::
Eq_const_from_CEA_curves(map<int, int> &nu, Gas_model &g, int iT)
    : Equilibrium_constant( iT ), nu_(nu), g_(g)
{
    // Check that the chemical species library has been initialised as required
    // for this particular equilibrium constant method
    if ( !chemical_species_library_initialised() ) {
    	ostringstream oss;
    	oss << "Eq_const_from_CEA_curves::Eq_const_from_CEA_curves()" << endl
    	    << "The chemical species library has not been initialised." << endl;
    	input_error( oss );
    }
}

Eq_const_from_CEA_curves::
~Eq_const_from_CEA_curves() {}

double
Eq_const_from_CEA_curves::
s_eval(const Gas_data &Q)
{
    double T = Q.T[iT_];
    double dG = 0.0;
    int nu_sum = 0;
    map<int, int>::const_iterator it;
    for ( it = nu_.begin(); it != nu_.end(); ++it ) {
	// We want to use G in J/mol, hence conversion based on molecular weight.
	double G = get_library_species_pointer(it->first)->eval_CEA_Gibbs_free_energy(T);
	dG += it->second * G *(g_.molecular_weight(it->first));
	nu_sum += it->second;
	// printf("G=%16.15e\n",G*(g_.molecular_weight(it->first)));
    }
    double K_p = exp(-dG/(PC_R_u*T));
    double K_c = K_p*pow(PC_P_atm/(PC_R_u*T), nu_sum);

    // printf("K_c=%16.15e\n", K_c);
    return K_c;
}

Equilibrium_constant* create_Eq_const_from_CEA_curves(lua_State *L,
						      std::map<int, int> &nu,
						      Gas_model &g)
{
    // use the lua_State to pull out the user defined iT
    int iT=get_int(L,-1,"iT");
    if ( iT<0 ) {
    	// The temperature of a specific species mode is being requested
    	string species = get_string(L,-1,"species");
    	string mode = get_string(L,-1,"mode");
    	iT = get_library_species_pointer_from_name(species)->get_mode_pointer_from_type(mode)->get_iT();
    }
    
    if ( iT<0 || iT >= g.get_number_of_modes() ) {
	ostringstream ost;
	ost << "create_Eq_const_from_CEA_curves():\n";
	ost << "Error in specification of rate controlling temperature index.\n";
	ost << "iT = " << iT << ", nmodes = " << g.get_number_of_modes() << endl;
	input_error(ost);
    }
    
    return new Eq_const_from_CEA_curves(nu, g, iT);
}
