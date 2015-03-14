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
import user_defined_effects;

/*
enum BCCode 
{ 
    slip_wall,
    adiabatic_wall,
    fixed_t_wall,
    full_face_exchange,
    mapped_cell,
    supersonic_in,
    subsonic_in,
    static_profile_in,
    fixed_p_out,
    extrapolate_out,
    ud_ghost_cells,
    ud_convective_flux,
    ud_interface_val,
    ud_diffusive_flux
}

BCCode type_code_from_name(string name)
{
    switch ( name ) {
    case "slip_wall", "slip-wall", "SlipWall", "slipwall":
	return BCCode.slip_wall;
    case "adiabatic", "adiabatic_wall", "adiabatic-wall", "AdiabaticWall",
	"adiabaticwall":
	return BCCode.adiabatic_wall;
    case "adjacent", "Adjacent", "full_face_exchange", "full-face-exchange",
	"FullFaceExchange", "fullfaceexchange":
	return BCCode.full_face_exchange;
    case "mapped_cell", "mapped-cell", "MappedCell", "mappedcell":
	return BCCode.mapped_cell;
    case "supersonic_in", "supersonic-in", "sup_in", "sup-in", "SupersonicIn",
	"SupIn", "supin":
	return BCCode.supersonic_in;
    case "subsonic_in", "subsonic-in", "sub_in", "sub-in",
	"SubsonicIn", "subsonicin":
	return BCCode.subsonic_in;
    case "static_profile_in", "static-profile-in", "static-profile", 
	"StaticProfileIn":
	return BCCode.static_profile_in;
    case "fixed_p_out", "fixed-p-out", "FixedPOut", "fixedpout":
	return BCCode.fixed_p_out;
    case "extrapolate_out", "extrapolate-out", "ExtrapolateOut", "extrapolateout":
	return BCCode.extrapolate_out;
    default: return BCCode.slip_wall;
    } // end switch
} // end type_code_from_name()
*/

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
	    newBC.initUserDefinedLuaState(gce, nicell, njcell, nkcell);
	    newBC.ghost_cell_data_available = true;
	}
    }
    // [TODO] We need to fill out the Action list
    // for other hook points.
    return newBC;
} // end make_BC_from_json()

// Boundary condition is abstract because no one is
// allowed to instantiate this barebones abstract class.

class BoundaryCondition {
    // [TODO] we need to redesign this so that we can use unstructured-grid blocks, eventually.
    // Presently, there are a lot of assumptions built in that are specific for structured-grid blocks.
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

    this(int id, int boundary, bool isWall=true, bool ghostCellDataAvailable=true, double _emissivity=0.0)
    {
	blk_id = id;
	which_boundary = boundary;
	is_wall = isWall;
	ghost_cell_data_available = ghostCellDataAvailable;
	emissivity = _emissivity;
    }

    void initUserDefinedLuaState(UserDefinedGhostCell udgc, size_t nicell, size_t njcell, size_t nkcell)
    {
	_lua = new LuaState();
	_lua.openLibs();
	setGlobalsInLuaState(nicell, njcell, nkcell);
	_lua.doFile(udgc.luafname);
	udgc.setLuaState(_lua);
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
    //    BoundaryInterfaceEffect[] preSpatialDerivAction;
    //    BoundaryFluxEffect[] postDiffFluxAction;

    override string toString() const
    {
	char[] repr;
	repr ~= "BoundaryCondition(";
	//repr ~= "type_code=" ~ to!string(type_code);
	repr ~= ")";
	return to!string(repr);
    }

    final void applyPreReconAction(double t, int tLevel)
    {
	foreach ( gce; preReconAction ) gce.apply(t, tLevel);
    }
    /*
    final void applyPostConvFluxAction(double t)
    {
	foreach ( bfe; postConvFluxAction ) bfe.apply(t);
    }

    final void applyPreSpatialDerivAction(double t)
    {
	foreach ( bie; preSpatialDerivAction ) bie.apply(t);
    }

    final void applyPostDiffFluxAction(double t)
    {
	foreach ( bfe; postSpatialDerivAction ) bfe.apply(t);
    }
    */
private:
    // We need a place to hold a lua state in memory that might
    // possibly be used by various stages of the boundary condition
    // action.
    LuaState _lua;
    void setGlobalsInLuaState(size_t nicell, size_t njcell, size_t nkcell)
    {
	_lua["blkId"] = blk_id;
	_lua["whichBoundary"] = which_boundary;
	_lua["n_species"] = GlobalConfig.gmodel.n_species;
	_lua["n_modes"] = GlobalConfig.gmodel.n_modes;
	_lua["nicell"] = nicell;
	_lua["njcell"] = njcell;
	_lua["nkcell"] = nkcell;
	_lua["NORTH"] = Face.north;
	_lua["EAST"] = Face.east;
	_lua["SOUTH"] = Face.south;
	_lua["WEST"] = Face.west;
	_lua["TOP"] = Face.top;
	_lua["BOTTOM"] = Face.bottom;
    }
    

} // end class BoundaryCondition

