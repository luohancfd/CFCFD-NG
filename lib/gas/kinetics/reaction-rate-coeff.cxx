// Author: Rowan J. Gollan
// Date: 29-Oct-2008
// Place: Hampton, Virginia, USA

#include <string>
#include <map>
#include <sstream>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"
#include "reaction-rate-coeff.hh"
#include "generalised-Arrhenius.hh"
#include "pressure-dependent-rate.hh"
#include "Park-nonequilibrium.hh"
#include "Macheret-dissociation.hh"
#include "MarroneTreanor-dissociation.hh"
#include "Knab-molecular-reaction.hh"

using namespace std;

Reaction_rate_coefficient* create_Reaction_rate_coefficient(lua_State *L, Gas_model &g)
{
    // Later implement as an object factory.
    map<string, Reaction_rate_coefficient* (*)(lua_State *, Gas_model &)> rr_coeff_models;
    rr_coeff_models.insert(pair<string, Reaction_rate_coefficient* (*)(lua_State *, Gas_model &)>("Arrhenius",
												  create_Generalised_Arrhenius_coefficient));
    rr_coeff_models.insert(pair<string, Reaction_rate_coefficient* (*)(lua_State *, Gas_model &)>("pressure dependent",
												  create_pressure_dependent_coefficient));
    rr_coeff_models.insert(pair<string, Reaction_rate_coefficient* (*)(lua_State *, Gas_model &)>("Park",
												  create_Park_nonequilibrium_coefficient));
    rr_coeff_models.insert(pair<string, Reaction_rate_coefficient* (*)(lua_State *, Gas_model &)>("Macheret",
												  create_Macheret_dissociation_coefficient));
    rr_coeff_models.insert(pair<string, Reaction_rate_coefficient* (*)(lua_State *, Gas_model &)>("MarroneTreanor",
												  create_MarroneTreanor_dissociation_coefficient));
    rr_coeff_models.insert(pair<string, Reaction_rate_coefficient* (*)(lua_State *, Gas_model &)>("Knab",
												  create_Knab_molecular_reaction_coefficient));
    string rmodel = get_string(L, -1, "model");
    
    if ( rr_coeff_models.find(rmodel) == rr_coeff_models.end() ) {
	ostringstream ost;
	ost << "create_Reaction_rate_coefficient():\n";
	ost << "Error in specification of reaction rate coefficient.\n";
	ost << "The selected model: " << rmodel << " is unknown.\n";
	ost << "The available models are: " << endl;
	map<string, Reaction_rate_coefficient* (*)(lua_State*, Gas_model &)>::const_iterator it;
	for ( it = rr_coeff_models.begin(); it != rr_coeff_models.end(); ++it ) {
	    ost << "   " << it->first << endl;
	}
	input_error(ost);
    }
    
    Reaction_rate_coefficient* rrc = rr_coeff_models[rmodel](L, g);

    if ( rrc == 0 ) {
	ostringstream ost;
	ost << "create_Reaction_rate_coefficient():\n";
	ost << "Error trying to create reaction rate coefficient of type: " << rmodel << endl;
	input_error(ost);
    }

    return rrc;
}

Reaction_rate_coefficient* create_Reaction_rate_coefficient(string cfile, string rate, Gas_model &g)
{
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    if( luaL_dofile(L, cfile.c_str()) != 0 ) {
	ostringstream ost;
	ost << "Reaction_rate_coefficient():\n";
	ost << "Error in reaction rate coefficient input file: " << cfile << endl;
	input_error(ost);
    }

    lua_getglobal(L, rate.c_str()); // bring table with rate to TOS

    Reaction_rate_coefficient *rrc = create_Reaction_rate_coefficient(L, g);

    lua_close(L);
    return rrc;
}

int Reaction_rate_coefficient::s_eval_from_T(const double T)
{
    cout << "Reaction_rate_coefficient::s_eval_from_T()" << endl
         << "Function not implemented for this reaction rate coefficient model." << endl
         << "Exiting program." << endl;
    exit( FAILURE );
}

