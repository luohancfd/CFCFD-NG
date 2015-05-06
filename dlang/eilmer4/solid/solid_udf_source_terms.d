/**
 * solid_udf_source_terms.d
 *
 * This module handles user-defined source terms
 * that the user might specify with a Lua script.
 *
 * Authors: RG & PJ
 * Date: 2015-05-06
 **/

module solid_udf_source_terms;

import std.stdio;
import std.string;
import luad.all;
import util.lua_service;

import solidfvcell;

LuaState initUDFSolidSourceTerms(string fname, int blkId)
{
    auto lua = initLuaState(fname);
    // Make blkId globally available in script
    lua["blkId"] = blkId;
    return lua;
}

void addUDFSourceTermsToSolidCell(LuaState lua, SolidFVCell cell, double t)
{
    // Push useful data into an arguments table
    auto args = lua.newTable();
    // Cell data
    args["x"] = cell.pos.x;
    args["y"] = cell.pos.y;
    args["z"] = cell.pos.z;
    args["vol"] = cell.volume;
    
    // Call solidSourceTerms function with (t, args)
    auto solidSourceTerms = lua.get!LuaFunction("solidSourceTerms");
    LuaObject[] ret = solidSourceTerms(t, args);
    if ( ret.length < 1 ) {
	string errMsg = "ERROR: There was a problem in the call to the user-defined source terms function\n";
	errMsg ~= "ERROR: for the solid domain. A single double value is expected,\n";
	errMsg ~= "ERROR: but nothing was returned.";
	throw new Exception(errMsg);
    }
    cell.Q = ret[0].to!double();
}
