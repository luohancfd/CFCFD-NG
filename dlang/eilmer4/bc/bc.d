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

BoundaryCondition make_BC_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    auto newBC = new BoundaryCondition(blk_id, boundary);
    // Assemble list of preReconAction effects
    auto preReconActionList = jsonData["pre_recon_action"].array;
    foreach ( jsonObj; preReconActionList ) {
	newBC.preReconAction ~= make_GCE_from_json(jsonObj, blk_id, boundary);
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

    final void applyPreReconAction(double t)
    {
	foreach ( gce; preReconAction ) gce.apply(t);
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

} // end class BoundaryCondition
