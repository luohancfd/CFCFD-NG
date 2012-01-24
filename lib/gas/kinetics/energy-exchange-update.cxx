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

    // FIXME: may want to do some parsing with a lua script file here
    
    // Parse the input file...
    if( luaL_dofile(L, cfile.c_str()) != 0 ) {
	ostringstream ost;
	ost << "create_Energy_exchange_update():\n";
	ost << "Error in input file: " << cfile << endl;
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
