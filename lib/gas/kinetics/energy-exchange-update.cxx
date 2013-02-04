// Author: Daniel F. Potter
// Date: 02-Dec-2009
// Place: UQ, St Lucia, Australia

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
#include "energy-exchange-update.hh"
#include "energy-exchange-ODE-update.hh"

using namespace std;

Energy_exchange_update* create_Energy_exchange_update(string cfile, Gas_model &g)
{
    // This is a candidate piece of code to be implemented
    // as an object factory.
    map<string, Energy_exchange_update* (*)(lua_State*, Gas_model&)> energy_exchange_updates;
    energy_exchange_updates.insert(pair<string, Energy_exchange_update* (*)(lua_State*, Gas_model&)>("energy exchange ODE", create_Energy_exchange_ODE_update));


    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    // Set up a species table
    lua_newtable(L);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	// This species table maps to C++ indices, because
	// it is used to setup the integer maps for
	// the reaction coefficients.
	lua_pushinteger(L, isp);
	lua_setfield(L, -2, g.species_name(isp).c_str());
	// Also add the reverse lookup
	lua_pushinteger(L, isp);
	lua_pushstring(L, g.species_name(isp).c_str());
	lua_settable(L, -3);
    }
    // Plus add a field 'size': no of species
    lua_pushinteger(L, g.get_number_of_species());
    lua_setfield(L, -2, "size");
    lua_setglobal(L, "species");
    
    // Setup a table of thermal modes
    lua_newtable(L);
    for ( int imode = 0; imode < g.get_number_of_modes(); ++imode ) {
	lua_newtable(L);
	for ( int ic = 0; ic < g.mode_no_components(imode); ++ic ) {
	    lua_pushinteger(L, ic);
	    lua_setfield(L, -2, g.mode_component_name(imode, ic).c_str());
	}
	lua_setfield(L, -2, g.mode_name(imode).c_str());
    }
    lua_setglobal(L, "modes");

    // Setup a table to find index of a given mode name
    lua_newtable(L);
    for ( int imode = 0; imode < g.get_number_of_modes(); ++imode) {
	lua_pushinteger(L, imode);
	lua_setfield(L, -2, g.mode_name(imode).c_str());
    }
    lua_setglobal(L, "mode_idx");

    cout << "Loading parser." << endl;
    // Path to reaction parsing script
    char *e3bin = getenv("E3BIN");
    string home;
    if ( e3bin == NULL ) {
	// Assume default location of $HOME/e3bin
	home.append(getenv("HOME")); home.append("/e3bin");
    }
    else {
	home.append(e3bin);
    }
    string script_file(home);
    script_file.append("/energy_exchange_parser.lua");
    
    if ( luaL_dofile(L, script_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "create_energy_exchange_update():\n";
	ost << "Error in loading script file: " << script_file << endl;
	ost << lua_tostring(L, -1) << endl;
	input_error(ost);
    }

    cout << "Parsing the input file." << endl;
    // Parse the input file...
    lua_getglobal(L, "main");
    lua_pushstring(L, cfile.c_str());
    if ( lua_pcall(L, 1, 0, 0) != 0 ) {
	ostringstream ost;
	ost << "create_Energy_exchange_update():\n";
	ost << "Error trying to load/parse energy exchange input file: " << cfile << endl;
	ost << lua_tostring(L, -1) << endl;
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
    if ( energy_exchange_updates.find(update) == energy_exchange_updates.end() ) {
	ostringstream ost;
	ost << "create_Energy_exchange_update():\n";
	ost << "Error in energy exchange scheme input file: " << cfile << endl;
	ost << "The energy exchange update type: " << update << " is unknown." << endl;
	ost << "The currently available update types are: " << endl;
	map<string, Energy_exchange_update* (*)(lua_State*, Gas_model&)>::const_iterator it;
	for ( it = energy_exchange_updates.begin(); it != energy_exchange_updates.end(); ++it ) {
	    ost << "   " << it->first << endl;
	}
	input_error(ost);
    }
    lua_pop(L, 1);

    Energy_exchange_update *eeu = energy_exchange_updates[update](L, g);

    if ( eeu == 0 ) {
	ostringstream ost;
	ost << "create_Energy_exchange_update():\n";
	ost << "Error trying to energy exchange update of type: " << update << endl;
	input_error(ost);
    }
    
    lua_close(L);
    
    return eeu;
}
