/**
 * user_defined_source_terms.d
 *
 * This module handles user-defined source terms
 * that the user might specify with a Lua script.
 *
 * Authors: RG & PJ
 * Date: 2015-03-17
 */

module user_defined_source_terms;

import std.conv;
import std.stdio;
import std.string;
import util.lua;
import util.lua_service;
import lua_helper;
import gas;
import fvcell;
import globalconfig;

void addUDFSourceTermsToCell(lua_State* L, FVCell cell, size_t gtl, double t, GasModel gmodel)
{
    size_t n_species = gmodel.n_species;
    size_t n_modes = gmodel.n_modes;

    // Push user function onto TOS
    lua_getglobal(L, "soureTerms");
    // Push sim_time onto TOS
    lua_pushnumber(L, t);
    // Push cell data into an args table and onto TOS
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    pushCellToTable(L, tblIdx, cell, gtl);
    // Call sourceTerms function with (t, args)
    int number_args = 2;
    int number_results = 1;

    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	    luaL_error(L, "error running user soure terms function: %s\n",
		       lua_tostring(L, -1));
    }

    // Grab values from user-returned table at TOS
    // For any missing values, put in 0.0
    cell.Q.mass += getNumberFromTable(L, -1, "mass", false, 0.0);
    cell.Q.momentum.refx += getNumberFromTable(L, -1, "momentum_x", false, 0.0);
    cell.Q.momentum.refy += getNumberFromTable(L, -1, "momentum_y", false, 0.0);
    cell.Q.momentum.refz += getNumberFromTable(L, -1, "momentum_z", false, 0.0);
    cell.Q.total_energy += getNumberFromTable(L, -1, "total_energy",false, 0.0);
    lua_getfield(L, -1, "species");
    if ( !lua_isnil(L, -1) ) {
	foreach ( isp; 0 .. n_species ) {
	    lua_rawgeti(L, -1, to!int(isp+1));
	    cell.Q.massf[isp] += lua_tonumber(L, -1);
	    lua_pop(L, 1);
	}
    }
    lua_pop(L, 1);

    lua_getfield(L, -1, "energies");
    if ( !lua_isnil(L, -1) ) {
	foreach ( imode; 0 ..  n_modes ) {
	    lua_rawgeti(L, -1, to!int(imode+1));
	    cell.Q.energies[imode] += lua_tonumber(L, -1);
	    lua_pop(L, 1);
	}
    }
    lua_pop(L, 1);
}
