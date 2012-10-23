// Author: Rowan J. Gollan
// Date: 14-Jul-2008

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <map>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

extern "C" {
#include <zlib.h>
}

#include "useful.h"
#include "lua_service.hh"

using namespace std;

void input_error(ostringstream &ost)
{
    cout << ost.str();
    cout << "Bailing out!\n";
    exit(BAD_INPUT_ERROR);
}

int check_for_integer(lua_State *L, string fname, const char *var,
		      const char *desc)
{
    int val = 0;
    lua_getglobal(L, var);
    if( lua_isnumber(L, -1) && (val = luaL_checkinteger(L, -1)) ) {
	if( val >= 1 ) {
	    return val;
	}
	else {
	    ostringstream ost;
	    ost << "Error in user-defined gas model file: " << fname << endl;
	    ost << "The " << desc << ", " << var << ", must be greater than or equal to 1.\n";
	    ost << "The specified value is: " << val << endl;
	    input_error(ost);
	}
    }
    else {
	ostringstream ost;
	ost << "Error in user-defined gas model file: " << fname << endl;
	ost << "The " << desc << " must be declared as '" << var << " = x'\n";
	ost << "where x is a positive integer.\n";
	input_error(ost);
    }
    // This should be unreachable...
    return val;
}

bool check_for_function(lua_State *L, const char *func)
{
    lua_getglobal(L, func);
    if( lua_isfunction(L, -1) == 1 ) {
	lua_pop(L, 1);
	return true;
    }
    else {
	lua_pop(L, 1);
	return false;
    }
}

void handle_lua_error(lua_State *L, const char *fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    fprintf(stderr, "\n");
    va_end(argp);
    lua_close(L);
    exit(-1);
}

int get_int(lua_State *L, int index, const char *key)
{
    // In Lua, all numbers are floats.  The best we can
    // do is cast to an integer.  We cannot ascertain
    // if the user actually put in an integer or a float.
    lua_getfield(L, index, key);
    if ( !lua_isnumber(L, -1) ) {
	ostringstream ost;
	ost << "An integer value was expected in field '" << key << "'." << endl;
	// Return stack to how it was
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }
    int val = luaL_checkint(L, -1);
    lua_pop(L, 1);
    return val;
}

int get_positive_int(lua_State *L, int index, const char *key)
{
    int val = get_int(L, index, key);
    if ( val <= 0 ) {
	ostringstream ost;
	ost << "A positive integer value was expected in field " << key << endl;
	throw runtime_error(ost.str());
    }
    return val;
}

double get_number(lua_State *L, int index, const char *key)
{
    lua_getfield(L, index, key);
    if ( !lua_isnumber(L, -1) ) {
	ostringstream ost;
	ost << "A number value was expected in field " << key << endl;
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }
    double val = lua_tonumber(L, -1);
    lua_pop(L, 1);
    return val;
}

double get_positive_number(lua_State *L, int index, const char *key)
{
    double val = get_number(L, index, key);
    if ( val <= 0.0 ) {
	ostringstream ost;
	ost << "A positive number value was expected in field " << key << endl;
	throw runtime_error(ost.str());
    }
    return val;
}

vector<double> get_vector(lua_State *L, int index, const char *key)
{
    lua_getfield(L, index, key);
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "A table was expected for field: " << key << endl;
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }

    // there should now be a vector at the top of the stack
    double val;
    vector<double> v;

    int len = (int)lua_objlen(L, -1); // cast from size_t
    for (int i = 1; i <= len; ++i) {
	lua_rawgeti(L, -1, i);
 	val = lua_tonumber(L, -1);
	v.push_back(val);
	lua_pop(L, 1);
    }
    lua_pop(L, 1);
    return v;
}

double get_negative_number(lua_State *L, int index, const char *key)
{
    double val = get_number(L, index, key);
    if ( val >= 0.0 ) {
	ostringstream ost;
	ost << "A negative number value was expected in field " << key << endl;
	throw runtime_error(ost.str());
    }
    return val;
}

string get_string(lua_State *L, int index, const char *key)
{
    lua_getfield(L, index, key);
    if ( !lua_isstring(L, -1) ) {
	ostringstream ost;
	ost << "A string was expected in field: " << key << endl;
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }
    string val = lua_tostring(L, -1);
    lua_pop(L, 1);
    return val;
}


double get_value(lua_State *L, int index, const char *key)
{
    double val = 0.0;
    
    lua_getfield(L, index, key);
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "A table was expected for field: " << key << endl;
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }
    
    lua_getfield(L, -1, "value");
    if ( !lua_isnumber(L, -1) ) {
	ostringstream ost;
	ost << "A number was expected in the 'value' field of entry: " << key << endl;
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }

    val = lua_tonumber(L, -1);

    lua_pop(L, 1); // pop 'value' off stack
    lua_pop(L, 1); // pop 'key' off stack
    return val;
}

double get_positive_value(lua_State *L, int index, const char *key)
{
    double val = get_value(L, index, key);
    if ( val <= 0.0 ) {
	ostringstream ost;
	ost << "A positive 'value' was expected in field " << key << endl;
	throw runtime_error(ost.str());
    }
    return val;
}

double get_negative_value(lua_State *L, int index, const char *key)
{
    double val = get_value(L, index, key);
    if ( val >= 0.0 ) {
	ostringstream ost;
	ost << "A negative 'value' was expected in field " << key << endl;
	throw runtime_error(ost.str());
    }
    return val;
}

int read_table_as_map(lua_State *L, int index, const char *key, map<int, int> &m)
{
    lua_getfield(L, index, key);
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "A table was expected for field: " << key << endl;
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }
    
    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i);

	lua_rawgeti(L, -1, 1);
	int index = luaL_checkint(L, -1);
	lua_pop(L, 1);

	lua_rawgeti(L, -1, 2);
	int val = luaL_checkint(L, -1);
	lua_pop(L, 1);

	m.insert(pair<int, int>(index, val));
	lua_pop(L, 1);
    }
    
    lua_pop(L, 1);
    return SUCCESS;
}

int read_table_as_map(lua_State *L, int index, const char *key, map<int, double> &m)
{
    lua_getfield(L, index, key);
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "A table was expected for field: " << key << endl;
	lua_pop(L, 1);
	throw runtime_error(ost.str());
    }
    
    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i);

	lua_rawgeti(L, -1, 1);
	int index = luaL_checkint(L, -1);
	lua_pop(L, 1);

	lua_rawgeti(L, -1, 2);
	double val = luaL_checknumber(L, -1);
	lua_pop(L, 1);

	m.insert(pair<int, double>(index, val));
	lua_pop(L, 1);
    }
    
    lua_pop(L, 1);
    return SUCCESS;
}

int do_gzfile(lua_State *L, string fname)
{
#   define BUF_SIZE 1024
    // Process a (possibly) zipped input file.
    gzFile ifp = gzopen(fname.c_str(), "r");
    if ( ifp == 0 ) {
	cout << "Problem opening lua file: " << fname << endl;
	cout << "Bailing out!\n";
	exit(BAD_INPUT_ERROR);
    }

    char buffer[BUF_SIZE];
    string lua_chunk;

    while ( gzeof(ifp) != 1 ) {
	char *buf = gzgets(ifp, buffer, BUF_SIZE);
	if ( buf == 0 ) {
	    // There was some problem, see if we
	    // can just press on until EOF.
	    continue;
	}
	lua_chunk.append(buf);
    }

    gzclose(ifp);

    return luaL_dostring(L, lua_chunk.c_str());
#   undef BUF_SIZE
}
