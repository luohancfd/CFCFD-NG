// bc/bc.d
// Base class for boundary condition objects, for use in Eilmer4
//
// Peter J. 2014-07-20 : first cut.
// RG & PJ  2015-12-03 : Decompose boundary conditions into lists of actions
//    

module bc;

import std.conv;
import std.json;
import std.stdio;
import std.string;
import luad.all;
import util.lua_service;
import json_helper;
import geom;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import fvinterface;
import fvcell;
import block;
import sblock;
import ghost_cell_effect;
import boundary_interface_effect;
import user_defined_effects;
import lua_helper;


BoundaryCondition make_BC_from_json(JSONValue jsonData, int blk_id, int boundary,
				    size_t nicell, size_t njcell, size_t nkcell)
{
    auto newBC = new BoundaryCondition(blk_id, boundary);
    // Assemble list of preReconAction effects
    auto preReconActionList = jsonData["pre_recon_action"].array;
    foreach ( jsonObj; preReconActionList ) {
	newBC.preReconAction ~= make_GCE_from_json(jsonObj, blk_id, boundary);
	// Some extra configuration in the case of a UserDefined bc.
	// We need to connect the Lua state back to the parent BC container.
	if ( newBC.preReconAction[$-1].type == "UserDefined" ) {
	    auto gce = to!UserDefinedGhostCell(newBC.preReconAction[$-1]);
	    if ( gce.luafname !in newBC.luaStates ) {
		newBC.initUserDefinedLuaState(gce.luafname, nicell, njcell, nkcell);
	    }
	    gce.setLuaState(newBC.luaStates[gce.luafname]);
	    newBC.ghost_cell_data_available = true;
	}
    }
    auto preSpatialDerivActionList = jsonData["pre_spatial_deriv_action"].array;
    foreach ( jsonObj; preSpatialDerivActionList ) {
	newBC.preSpatialDerivAction ~= make_BIE_from_json(jsonObj, blk_id, boundary);
	// [TODO] need to think about a way to connect to the user_defined BC
	// as appropriate, or do we just have separate interpreters?
    }
    // [TODO] We need to fill out the Action lists for other hook points.
    return newBC;
} // end make_BC_from_json()


class BoundaryCondition {
    // Boundary condition is built from composable pieces.
public:
    // Location of the boundary condition.
    int blk_id; // index of the structured-grid block to which this BC is applied
    int which_boundary; // identity/index of the relevant boundary

    // Nature of the boundary condition that may be checked 
    // by other parts of the CFD code.
    // BCCode type_code = BCCode.slip_wall;
    bool is_wall = true;
    bool ghost_cell_data_available = true;
    double emissivity = 0.0;

    // Storage for various LuaStates in the case of user-defined boundary
    // conditions. We store these in the BoundaryCondition object
    // so that the composable pieces might reference the same LuaState.
    // Conversely, for each different Lua file, a different LuaState exists.
    LuaState[string] luaStates;

    this(int id, int boundary, bool isWall=true, bool ghostCellDataAvailable=true, double _emissivity=0.0)
    {
	blk_id = id;
	which_boundary = boundary;
	is_wall = isWall;
	ghost_cell_data_available = ghostCellDataAvailable;
	emissivity = _emissivity;
    }

    void initUserDefinedLuaState(string luafname, size_t nicell, size_t njcell, size_t nkcell)
    {
	luaStates[luafname] = new LuaState();
	luaStates[luafname].openLibs();
	setGlobalsInLuaState(luaStates[luafname], nicell, njcell, nkcell);
	luaStates[luafname].doFile(luafname);
    }
    // Action lists.
    // The BoundaryCondition is called at four stages in a global timestep.
    // Those stages are:
    // 1. pre reconstruction
    // 2. post convective flux evaluation
    // 3. pre spatial derivative estimate
    // 4. post diffusive flux evaluation
    // Note the object may be called more than 4 times depending
    // on the type of time-stepping used to advance the solution.
    // At each of these stages, a series of effects are applied in order
    // with the end goal to leave the boundary values in an appropriate
    // state. We will call this series of effects an action.
    GhostCellEffect[] preReconAction;
    //    BoundaryFluxEffect[] postConvFluxAction;
    BoundaryInterfaceEffect[] preSpatialDerivAction;
    //    BoundaryFluxEffect[] postDiffFluxAction;

    override string toString() const
    {
	char[] repr;
	repr ~= "BoundaryCondition(";
	if ( preReconAction.length > 0 ) {
	    repr ~= "preReconAction=[" ~ to!string(preReconAction[0]);
	    foreach (i; 1 .. preReconAction.length) {
		repr ~= ", " ~ to!string(preReconAction[i]);
	    }
	    repr ~= "]";
	}
	repr ~= ", ";
	if ( preSpatialDerivAction.length > 0 ) {
	    repr ~= "preSpatialDerivAction=[" ~ to!string(preSpatialDerivAction[0]);
	    foreach (i; 1 .. preSpatialDerivAction.length) {
		repr ~= ", " ~ to!string(preSpatialDerivAction[i]);
	    }
	    repr ~= "]";
	}
	repr ~= ")";
	return to!string(repr);
    }

    final void applyPreReconAction(double t, int gtl, int ftl)
    {
	foreach ( gce; preReconAction ) gce.apply(t, gtl, ftl);
    }
    /*
    final void applyPostConvFluxAction(double t)
    {
	foreach ( bfe; postConvFluxAction ) bfe.apply(t);
    }
    */
    final void applyPreSpatialDerivAction(double t, int gtl, int ftl)
    {
	foreach ( bie; preSpatialDerivAction ) bie.apply(t, gtl, ftl);
    }
    /*
    final void applyPostDiffFluxAction(double t)
    {
	foreach ( bfe; postSpatialDerivAction ) bfe.apply(t);
    }
    */
private:
    void setGlobalsInLuaState(LuaState lua, size_t nicell, size_t njcell, size_t nkcell)
    {
	lua["blkId"] = blk_id;
	lua["whichBoundary"] = which_boundary;
	lua["n_species"] = GlobalConfig.gmodel.n_species;
	lua["n_modes"] = GlobalConfig.gmodel.n_modes;
	lua["nicell"] = nicell;
	lua["njcell"] = njcell;
	lua["nkcell"] = nkcell;
	lua["north"] = Face.north;
	lua["east"] = Face.east;
	lua["south"] = Face.south;
	lua["west"] = Face.west;
	lua["top"] = Face.top;
	lua["bottom"] = Face.bottom;
	lua["sampleFlow"] = &luafn_sampleFlow;
    }
    

} // end class BoundaryCondition

