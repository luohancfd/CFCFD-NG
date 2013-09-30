// Author: Rowan J. Gollan
// Version: 07-July-2008 
//

#ifndef LSERVICE_GAS_DATA_HH
#define LSERVICE_GAS_DATA_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "gas-model.hh"

int push_gas_data_as_table(lua_State *L, const Gas_data &Q);
void set_gas_data_at_table(lua_State *L, int tindex, const Gas_data &Q);
void get_table_as_gas_data(lua_State *L, Gas_data &Q);
int push_vector_as_table(lua_State *L, const std::vector<double> &vec);
void get_table_as_vector(lua_State *L, std::vector<double> &vec);
int push_matrix_as_table(lua_State *L, const matrix &mat);
void get_table_as_matrix(lua_State *L, matrix &mat);
int create_empty_gas_table(lua_State *L, Gas_model &gmodel);

#endif
