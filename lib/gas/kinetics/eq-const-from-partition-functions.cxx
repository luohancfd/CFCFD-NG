// Author: Daniel F. Potter
// Date: 30-Oct-2012
// Place: Goettingen, Germany
//

#include <iostream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "../models/physical_constants.hh"
#include "../models/chemical-species-library.hh"

#include "eq-const-from-partition-functions.hh"

using namespace std;

Eq_const_from_partition_functions::
Eq_const_from_partition_functions(map<int, int> &nu, Gas_model &g, int iT)
    : Equilibrium_constant( iT ), nu_(nu), g_(g) {}

Eq_const_from_partition_functions::
~Eq_const_from_partition_functions()
{}

double
Eq_const_from_partition_functions::
s_eval(const Gas_data &Q)
{
    // pull out the rate controlling temperature
    double T = Q.T[iT_];

    // calculate K_c as product of partition functions
    double K_c = 1.0;
    map<int, int>::const_iterator it;
    for ( it = nu_.begin(); it != nu_.end(); ++it ) {
	Chemical_species * X = get_library_species_pointer(it->first);
	K_c *= pow( X->eval_partition_function(T), it->second);
    }

    // printf("K_c=%16.15e\n", K_c);
    return K_c;
}

Equilibrium_constant*
create_Eq_const_from_partition_functions(lua_State *L,
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
    
    return new Eq_const_from_partition_functions(nu, g, iT);
}
