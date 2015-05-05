/**
 * user_defined_source_terms.d
 *
 * This module handles user-defined source terms
 * that the user might specify with a Lua script.
 *
 * Authors: RG & PJ
 * Date: 2015-03-17
 */

import std.stdio;
import std.string;
import luad.all;
import util.lua_service;
import lua_helper;
import gas;
import fvcell;
import globalconfig;

// Keep a hold of the LuaState so that it
// is reentrant during simulation.
static LuaState lua;

void initUDFSourceTerms(string fname)
{
    lua = initLuaState(fname);
}

void addUDFSourceTermsToCell(FVCell cell, size_t gtl, double t, GasModel gmodel)
{
    size_t n_species = gmodel.n_species;
    size_t n_modes = gmodel.n_modes;

    // Push cell data into an args table.
    auto args = lua.newTable();
    pushCellToLuaTable(cell, gtl, args);

    // Call sourceTerms function with (t, args)
    auto sourceTerms = lua.get!LuaFunction("sourceTerms");
    LuaObject[] ret = sourceTerms(t, args);
    if ( ret.length < 1 ) {
	string errMsg = "ERROR: There was a problem in the call to the user-defined source terms function.\n";
	errMsg ~= format("ERROR: A table is expected, but nothing was returned.");
	throw new Exception(errMsg);
    }
    auto tab = ret[0].to!LuaTable();
    // For any missing values, put in 0.0
    cell.Q.mass += getDouble(tab, "mass", 0.0);
    cell.Q.momentum.refx += getDouble(tab, "momentum_x", 0.0);
    cell.Q.momentum.refy += getDouble(tab, "momentum_y", 0.0);
    cell.Q.momentum.refz += getDouble(tab, "momentum_z", 0.0);
    cell.Q.total_energy += getDouble(tab, "total_energy", 0.0);
    foreach ( isp; 0..n_species ) {
	try {
	    cell.Q.massf[isp] += tab.get!double("species", isp+1);
	}
	catch (Exception) {
	    // Do nothing, add no value if we couldn't find any.
	}
    }
    foreach ( i; 0..n_modes ) {
	try {
	    cell.Q.energies[i] += tab.get!double("energies", i+1);
	}
	catch (Exception) {
	    // Do nothing, add no value if we couldn't find any.
	}
    }
}
