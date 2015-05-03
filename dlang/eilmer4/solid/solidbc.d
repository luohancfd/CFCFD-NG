/**
 * solidbc.d
 *
 * Author: Rowan G. and Peter J.
 */

module solidbc;

import std.json;
import luad.all;
import util.lua_service;
import std.conv;

import geom;
import json_helper;
import solid_boundary_interface_effect;
import solid_boundary_flux_effect;

SolidBoundaryCondition makeSolidBCFromJson(JSONValue jsonData, int blk_id, int boundary,
					   size_t nicell, size_t njcell, size_t nkcell)
{
    auto newBC = new SolidBoundaryCondition(blk_id, boundary);
    // Assemble list of preSpatialDerivAction effects
    auto preSpatialDerivActionList = jsonData["pre_spatial_deriv_action"].array;
    foreach ( jsonObj; preSpatialDerivActionList ) {
	newBC.preSpatialDerivAction ~= makeSolidBIEfromJson(jsonObj, blk_id, boundary);
	if ( newBC.preSpatialDerivAction[$-1].type == "UserDefined" ) {
	    auto sbie = to!SolidBIE_UserDefined(newBC.preSpatialDerivAction[$-1]);
	    newBC.initUserDefinedLuaState(sbie, nicell, njcell, nkcell);
	}
    }
    return newBC;
}

class SolidBoundaryCondition {
public:
    int blkId;
    int whichBoundary;
    //    bool setsFluxDirectly;

    this(int _blkId, int boundary)
    {
	blkId = _blkId;
	whichBoundary = boundary;
    }

    final void applyPreSpatialDerivAction(double t, int tLevel)
    {
	foreach ( sie; preSpatialDerivAction ) sie.apply(t, tLevel);
    }

    final void applyPostFluxAction(double t, int tLevel)
    {
	foreach ( sfe; postFluxAction ) sfe.apply(t, tLevel);
    }

    final void initUserDefinedLuaState(SolidBIE_UserDefined sbie, size_t nicell, size_t njcell, size_t nkcell)
    {
	_lua = new LuaState();
	_lua.openLibs();
	setGlobalsInLuaState(nicell, njcell, nkcell);
	_lua.doFile(sbie.luafname);
	sbie.setLuaState(_lua);
    }

    SolidBoundaryInterfaceEffect[] preSpatialDerivAction;
    SolidBoundaryFluxEffect[] postFluxAction;

private:
    LuaState _lua;
    void setGlobalsInLuaState(size_t nicell, size_t njcell, size_t nkcell)
    {
	_lua["blkId"] = blkId;
	_lua["whichBoundary"] = whichBoundary;
	_lua["nicell"] = nicell;
	_lua["njcell"] = njcell;
	_lua["nkcell"] = nkcell;
	_lua["north"] = Face.north;
	_lua["east"] = Face.east;
	_lua["south"] = Face.south;
	_lua["west"] = Face.west;
	_lua["top"] = Face.top;
	_lua["bottom"] = Face.bottom;
    }

}
