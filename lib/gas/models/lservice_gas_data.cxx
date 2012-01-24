// Author: Rowan J. Gollan
// Date: 04-July-2008

#include <sstream>
#include <iostream>

#include "../../util/source/lua_service.hh"
#include "lservice_gas_data.hh"

using namespace std;

int push_gas_data_as_table(lua_State *L, const Gas_data &Q)
{
    lua_newtable(L);
    int tindex = lua_gettop(L);
    
    lua_pushnumber(L, Q.rho);
    lua_setfield(L, tindex, "rho");

    push_vector_as_table(L, Q.e);
    lua_setfield(L, tindex, "e");

    lua_pushnumber(L, Q.p);
    lua_setfield(L, tindex, "p");

    lua_pushnumber(L, Q.a);
    lua_setfield(L, tindex, "a");

    push_vector_as_table(L, Q.T);
    lua_setfield(L, tindex, "T");

    lua_pushnumber(L, Q.mu);
    lua_setfield(L, tindex, "mu");
    
    push_vector_as_table(L, Q.k);
    lua_setfield(L, tindex, "k");

    push_matrix_as_table(L, Q.D_AB);
    lua_setfield(L, tindex, "D_AB");

    push_vector_as_table(L, Q.massf);
    lua_setfield(L, tindex, "massf");

    return 1;
}

void set_gas_data_at_table(lua_State *L, int tindex, const Gas_data &Q)
{
    lua_pushnumber(L, Q.rho);
    lua_setfield(L, tindex, "rho");

    push_vector_as_table(L, Q.e);
    lua_setfield(L, tindex, "e");

    lua_pushnumber(L, Q.p);
    lua_setfield(L, tindex, "p");

    lua_pushnumber(L, Q.a);
    lua_setfield(L, tindex, "a");

    push_vector_as_table(L, Q.T);
    lua_setfield(L, tindex, "T");

    lua_pushnumber(L, Q.mu);
    lua_setfield(L, tindex, "mu");
    
    push_vector_as_table(L, Q.k);
    lua_setfield(L, tindex, "k");

    push_matrix_as_table(L, Q.D_AB);
    lua_setfield(L, tindex, "D_AB");

    push_vector_as_table(L, Q.massf);
    lua_setfield(L, tindex, "massf");
}

void get_table_as_gas_data(lua_State *L, Gas_data &Q)
{
    // Assuming table is TOS
    lua_getfield(L, -1, "rho");
    Q.rho = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, -1, "e");
    get_table_as_vector(L, Q.e);
    lua_pop(L, 1);

    lua_getfield(L, -1, "p");
    Q.p = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, -1, "a");
    Q.a = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, -1, "T");
    get_table_as_vector(L, Q.T);
    lua_pop(L, 1);

    lua_getfield(L, -1, "mu");
    Q.mu = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    lua_getfield(L, -1, "k");
    get_table_as_vector(L, Q.k);
    lua_pop(L, 1);

    lua_getfield(L, -1, "D_AB");
    get_table_as_matrix(L, Q.D_AB);
    lua_pop(L, 1);

    lua_getfield(L, -1, "massf");
    get_table_as_vector(L, Q.massf);
    lua_pop(L, 1);

}

int push_vector_as_table(lua_State *L, const vector<double> &vec)
{
    lua_newtable(L);
    int t = lua_gettop(L);
    for( size_t i = 0; i < vec.size(); ++i ) {
	lua_pushnumber(L, vec[i]);
	lua_rawseti(L, t, i+1);
    }
    return t;
}

void get_table_as_vector(lua_State *L, vector<double> &vec)
{
    for( size_t i = 0; i < vec.size(); ++i ) {
	lua_rawgeti(L, -1, i+1);
	vec[i] = luaL_checknumber(L, -1);
	lua_pop(L, 1);
    }
}

int push_matrix_as_table(lua_State *L, const matrix &mat)
{
    lua_newtable(L);
    int t = lua_gettop(L);
    for( size_t i = 0; i < mat.size(); ++i ) {
	lua_newtable(L);
	int u = lua_gettop(L);
	for( size_t j = 0; j < mat[i].size(); ++j ) {
	    lua_pushnumber(L, mat[i][j]);
	    lua_rawseti(L, u, j+1);
	}
	lua_rawseti(L, t, i+1);
    }
    return t;
}

void get_table_as_matrix(lua_State *L, matrix &mat)
{
    for( size_t i = 0; i < mat.size(); ++i ) {
	lua_rawgeti(L, -1, i+1);
	for( size_t j = 0; j < mat[i].size(); ++j ) {
	    lua_rawgeti(L, -1, j+1);
	    mat[i][j] = luaL_checknumber(L, -1);
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);
    }
}

int create_empty_gas_table(lua_State *L)
{
    int nsp = 1;
    int nmodes = 1;
    int narg = lua_gettop(L);
    if ( narg == 0 ) {
	nsp = 1;
	nmodes = 1;
    }
    else if ( narg == 1 ) {
	nsp = luaL_checkint(L, 1);
	nmodes = 1;
    }
    else if ( narg == 2 ) {
	nsp = luaL_checkint(L, 1);
	nmodes = luaL_checkint(L, 2);
    }
    else {
	ostringstream ost;
	ost << "Wrong number of arguments to function create_empty_gas_table()\n";
	ost << "0, 1 or 2 expected. " << narg << " received.\n";
	input_error(ost);
    }

    lua_newtable(L);
    int tindex = lua_gettop(L);
    
    lua_pushnumber(L, 0.0);
    lua_setfield(L, tindex, "rho");

    vector<double> mv(nmodes, 0.0);
    push_vector_as_table(L, mv);
    lua_setfield(L, tindex, "e");

    lua_pushnumber(L, 0.0);
    lua_setfield(L, tindex, "p");

    lua_pushnumber(L, 0.0);
    lua_setfield(L, tindex, "a");

    push_vector_as_table(L, mv);
    lua_setfield(L, tindex, "T");

    lua_pushnumber(L, 0.0);
    lua_setfield(L, tindex, "mu");
    
    push_vector_as_table(L, mv);
    lua_setfield(L, tindex, "k");

    matrix D_AB;
    D_AB.resize(nsp);
    for ( int i = 0; i < nsp; ++i ) D_AB[i].resize(nsp);
    push_matrix_as_table(L, D_AB);
    lua_setfield(L, tindex, "D_AB");

    vector<double> mfv(nsp, 0.0);
    push_vector_as_table(L, mfv);
    lua_setfield(L, tindex, "massf");

    return 1;
}
