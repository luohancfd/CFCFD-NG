// Author: Rowan J. Gollan
// Date: 29-Oct-2008
// Place: Hampton, Viriginia, USA
//

#include <map>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../util/source/lua_service.hh"

#include "equilibrium-constant.hh"
#include "eq-const-from-thermo.hh"
#include "eq-const-from-CEA-curves.hh"
#include "eq-const-from-partition-functions.hh"

using namespace std;

Equilibrium_constant::Equilibrium_constant( int iT )
: iT_( iT ) {}

Equilibrium_constant*
create_Equilibrium_constant(lua_State *L,
			    map<int, int> &nu,
			    Gas_model &g)
{
    // Later implement as an object factory.
    map<string, Equilibrium_constant* (*)(lua_State*, map<int, int>&, Gas_model&)> ec_models;
    ec_models.insert(pair<string, Equilibrium_constant* (*)(lua_State*, map<int,int>&, Gas_model&)>("from thermo",
												   create_Eq_const_from_thermo));
    ec_models.insert(pair<string, Equilibrium_constant* (*)(lua_State*, map<int,int>&, Gas_model&)>("from CEA curves",
												   create_Eq_const_from_CEA_curves));
    ec_models.insert(pair<string, Equilibrium_constant* (*)(lua_State*, map<int,int>&, Gas_model&)>("from partition functions",
												   create_Eq_const_from_partition_functions));
    string model = get_string(L, -1, "model");
    
    if ( ec_models.find(model) == ec_models.end() ) {
	ostringstream ost;
	ost << "create_equilibrium_constant():\n";
	ost << "Error in specification of equilibrium constant.\n";
	ost << "The selected model: " << model << " is unknown.\n";
	ost << "The available models are: " << endl;
	map<string, Equilibrium_constant* (*)(lua_State*, map<int, int>&, Gas_model&)>::const_iterator it;
	for ( it = ec_models.begin(); it != ec_models.end(); ++it ) {
	    ost << "   " << it->first << endl;
	}
    }

    Equilibrium_constant *ec = ec_models[model](L, nu, g);

    if ( ec == 0 ) {
	ostringstream ost;
	ost << "create_equilibrium_constant():\n";
	ost << "Error trying to create equilibrium constant of type: " << model << endl;
	input_error(ost);
    }

    return ec;
}
