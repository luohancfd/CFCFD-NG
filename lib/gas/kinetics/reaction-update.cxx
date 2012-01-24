// Author: Rowan J. Gollan
// Date: 17-Oct-2008
// Place: Hampton, Virginia, USA

#include <cstdlib>
#include <string>
#include <map>
#include <sstream>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../util/source/lua_service.hh"
#include "reaction-update.hh"
#include "chemical-kinetic-ODE-update.hh"
#include "chemical-kinetic-ODE-MC-update.hh"

using namespace std;

Reaction_update* create_Reaction_update(string cfile, Gas_model &g)
{
    // This is a candidate piece of code to be implemented
    // as an object factory.
    map<string, Reaction_update* (*)(lua_State*, Gas_model&)> reaction_updates;
    reaction_updates.insert(pair<string, Reaction_update* (*)(lua_State*, Gas_model&)>("chemical kinetic ODE", create_Chemical_kinetic_ODE_update));
    reaction_updates.insert(pair<string, Reaction_update* (*)(lua_State*, Gas_model&)>("chemical kinetic ODE MC", create_Chemical_kinetic_ODE_MC_update));

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    // Set up species table
    lua_newtable(L);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	// This species table maps to C++ indices, because
	// it is used to setup the integer maps for
	// the reaction coefficients.
	lua_pushinteger(L, isp);
	lua_setfield(L, -2, g.species_name(isp).c_str());
    }
    // Plus add a field 'size': no of species
    lua_pushinteger(L, g.get_number_of_species());
    lua_setfield(L, -2, "size");
    lua_setglobal(L, "species");

    // Path to reaction parsing script
    string home(getenv("HOME"));
    string script_file(home);
    script_file.append("/e3bin/reaction_parser.lua");

    if ( luaL_dofile(L, script_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error in loading script file: " << script_file << endl;
	input_error(ost);
    }
    
    // Parse the input file...
    lua_getglobal(L, "main");
    lua_pushstring(L, cfile.c_str());
    if ( lua_pcall(L, 1, 0, 0) != 0 ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error trying to load reaction scheme file: " << cfile << endl;
	ost << "Lua error message: " << lua_tostring(L, -1) << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, "scheme_t");
    if ( lua_isnil(L, -1) ) {
	cout << "scheme_t not set globally.";
	exit(1);
    }
    
    lua_getfield(L, -1, "update");
    string update(luaL_checkstring(L, -1));
    
    lua_pop(L, 1);
    if ( reaction_updates.find(update) == reaction_updates.end() ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error in reaction scheme input file: " << cfile << endl;
	ost << "The reaction update type: " << update << " is unknown." << endl;
	ost << "The currently available update types are: " << endl;
	map<string, Reaction_update* (*)(lua_State*, Gas_model&)>::const_iterator it;
	for ( it = reaction_updates.begin(); it != reaction_updates.end(); ++it ) {
	    ost << "   " << it->first << endl;
	}
	input_error(ost);
    }
    lua_pop(L, 1);

    Reaction_update *ru = reaction_updates[update](L, g);

    if ( ru == 0 ) {
	ostringstream ost;
	ost << "create_reaction_update():\n";
	ost << "Error trying to reaction update of type: " << update << endl;
	input_error(ost);
    }
    
    lua_close(L);

    return ru;
}
