// lua_service_for_e3.cxx
// Author: Rowan J. Gollan
// Date: 13-Mar-2013
// Place: The University of Queensland
//
// History: 13-Mar-2013
//          Refactored some service functions from bc_user_defined.hh
//

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "kernel.hh"
#include "block.hh"
#include "cell.hh"

int luafn_sample_flow(lua_State *L)
{
    // Get arguments from stack.
    int jb = lua_tointeger(L, 1);
    int i = lua_tointeger(L, 2);
    int j = lua_tointeger(L, 3);
    int k = lua_tointeger(L, 4);

    Gas_model *gmodel = get_gas_model_ptr();
    Block *bdp= get_block_data_ptr(jb);
    FV_Cell *cell = bdp->get_cell(i,j,k);
    FlowState &fs = *(cell->fs);

    // Return the interesting data in a table with named fields.
    lua_newtable(L); // creates a table that is now at the TOS
    lua_pushnumber(L, cell->pos.x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, cell->pos.y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, cell->pos.z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, cell->volume); lua_setfield(L, -2, "vol");
    lua_pushnumber(L, fs.gas->p); lua_setfield(L, -2, "p");
    lua_pushnumber(L, fs.gas->T[0]); lua_setfield(L, -2, "T_wall");
    lua_pushnumber(L, fs.gas->rho); lua_setfield(L, -2, "rho"); 
    lua_pushnumber(L, fs.vel.x); lua_setfield(L, -2, "u"); 
    lua_pushnumber(L, fs.vel.y); lua_setfield(L, -2, "v");
    lua_pushnumber(L, fs.vel.z); lua_setfield(L, -2, "w");
    lua_pushnumber(L, fs.gas->a); lua_setfield(L, -2, "a");
    lua_pushnumber(L, fs.mu_t); lua_setfield(L, -2, "mu_t");
    lua_pushnumber(L, fs.k_t); lua_setfield(L, -2, "k_t");
    lua_pushnumber(L, fs.tke); lua_setfield(L, -2, "tke");
    lua_pushnumber(L, fs.omega); lua_setfield(L, -2, "omega");
    lua_pushnumber(L, fs.S); lua_setfield(L, -2, "S");
    lua_newtable(L); // A table for the individual temperatures
    int nmodes = gmodel->get_number_of_modes();
    for ( int imode = 0; imode < nmodes; ++imode ) {
	lua_pushinteger(L, imode);
	lua_pushnumber(L, fs.gas->T[imode]);
	lua_settable(L, -3);
    }
    // At this point, the table of temperatures should be TOS.
    lua_setfield(L, -2, "T");
    lua_newtable(L); // Another table for the mass fractions
    int nsp = gmodel->get_number_of_species();
    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_pushinteger(L, isp);
	lua_pushnumber(L, fs.gas->massf[isp]);
	lua_settable(L, -3);
    }
    // At this point, the table of mass fractions should be TOS.
    lua_setfield(L, -2, "massf");
    // Parameters for computing density and energy residuals
    lua_pushnumber(L, cell->U_old->mass); lua_setfield(L, -2, "rho_old");
    lua_pushnumber(L, cell->U->total_energy); lua_setfield(L, -2, "rE");
    lua_pushnumber(L, cell->U_old->total_energy); lua_setfield(L, -2, "rE_old");
    
    return 1; // Just the one table is left on Lua stack
}

int luafn_locate_cell(lua_State *L)
{
    double x = lua_tonumber(L, 1);
    double y = lua_tonumber(L, 2);
    double z = 0.0;
    global_data *gdp = get_global_data_ptr();
    if ( gdp->dimensions == 3 ) {
	z = lua_tonumber(L, 3);
    }
    int jb, i, j, k, found_cell;
    found_cell = locate_cell(x, y, z, &jb, &i, &j, &k);
    if ( !found_cell ) {
	found_cell = find_nearest_cell(x, y, z, &jb, &i, &j, &k);
	if ( !found_cell ) {
	    printf("luafn_locate_cell(): no cells near pos=(%g,%g,%g)", x, y, z);
	}
    }
    lua_pushinteger(L, jb);
    lua_pushinteger(L, i);
    lua_pushinteger(L, j);
    lua_pushinteger(L, k);
    lua_pushinteger(L, found_cell);
    return 5; // We leave jb, i, j, k and found_cell on the stack.
}

int register_luafns(lua_State *L)
{
    lua_pushcfunction(L, luafn_sample_flow);
    lua_setglobal(L, "sample_flow");
    lua_pushcfunction(L, luafn_locate_cell);
    lua_setglobal(L, "locate_cell");
    return 0;
}
