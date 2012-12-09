// Author: Rowan J. Gollan
// Date: 14-Jul-2008

#ifndef LUA_SERVICE_HH
#define LUA_SERVICE_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include <sstream>
#include <string>
#include <map>
#include <vector>

void input_error(std::ostringstream &ost);
int check_for_integer(lua_State *L, std::string fname, const char *var,
		      const char *desc);
bool check_for_function(lua_State *L, const char *func);
void handle_lua_error(lua_State *L, const char *fmt, ...);

int get_int(lua_State *L, int index, const char *key);
int get_positive_int(lua_State *L, int index, const char *key);

double get_number(lua_State *L, int index, const char *key);
double get_positive_number(lua_State *L, int index, const char *key);
double get_negative_number(lua_State *L, int index, const char *key);
// Remember Lua indexes from 1 while C++/Python from 0.
std::vector<double> get_vector(lua_State *L, int index, const char *key);

std::string get_string(lua_State *L, int index, const char *key);

double get_value(lua_State *L, int index, const char *key);
double get_positive_value(lua_State *L, int index, const char *key);
double get_negative_value(lua_State *L, int index, const char *key);
bool get_boolean(lua_State *L, int index, const char *key);

int read_table_as_map(lua_State *L, int index, const char *key, std::map<int, int> &m);
int read_table_as_map(lua_State *L, int index, const char *key, std::map<int, double> &m);

int do_gzfile(lua_State *L, std::string fname);

#endif
