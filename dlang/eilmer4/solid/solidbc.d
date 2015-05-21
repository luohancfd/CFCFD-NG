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
    auto setsFluxDirectly = getJSONbool(jsonData, "sets_flux_directly", false);
    auto newBC = new SolidBoundaryCondition(blk_id, boundary, setsFluxDirectly);
    // Assemble list of preSpatialDerivAction effects
    auto preSpatialDerivActionList = jsonData["pre_spatial_deriv_action"].array;
    foreach ( jsonObj; preSpatialDerivActionList ) {
	newBC.preSpatialDerivAction ~= makeSolidBIEfromJson(jsonObj, blk_id, boundary);
	if ( newBC.preSpatialDerivAction[$-1].type == "UserDefined" ) {
	    auto sbie = to!SolidBIE_UserDefined(newBC.preSpatialDerivAction[$-1]);
	    if ( sbie.luafname !in newBC.luaStates ) {
		newBC.initUserDefinedLuaState(sbie.luafname, nicell, njcell, nkcell);
	    }
	    sbie.setLuaState(newBC.luaStates[sbie.luafname]);
	}
    }
    return newBC;
}

class SolidBoundaryCondition {
public:
    int blkId;
    int whichBoundary;
    bool setsFluxDirectly;
    LuaState[string] luaStates;
    SolidBoundaryInterfaceEffect[] preSpatialDerivAction;
    //SolidBoundaryFluxEffect[] postFluxAction;

    this(int _blkId, int boundary, bool _setsFluxDirectly)
    {
	blkId = _blkId;
	whichBoundary = boundary;
	setsFluxDirectly = _setsFluxDirectly;
    }

    final void applyPreSpatialDerivAction(double t, int tLevel)
    {
	foreach ( sie; preSpatialDerivAction ) sie.apply(t, tLevel);
    }

    final void applyPostFluxAction(double t, int tLevel)
    {
	//foreach ( sfe; postFluxAction ) sfe.apply(t, tLevel);
    }

    final void initUserDefinedLuaState(string luafname, size_t nicell, size_t njcell, size_t nkcell)
    {
	luaStates[luafname] = new LuaState();
	luaStates[luafname].openLibs();
	setGlobalsInLuaState(luaStates[luafname], nicell, njcell, nkcell);
	luaStates[luafname].doFile(luafname);
    }

private:
    void setGlobalsInLuaState(LuaState lua, size_t nicell, size_t njcell, size_t nkcell)
    {
	lua["blkId"] = blkId;
	lua["whichBoundary"] = whichBoundary;
	lua["nicell"] = nicell;
	lua["njcell"] = njcell;
	lua["nkcell"] = nkcell;
	lua["north"] = Face.north;
	lua["east"] = Face.east;
	lua["south"] = Face.south;
	lua["west"] = Face.west;
	lua["top"] = Face.top;
	lua["bottom"] = Face.bottom;
    }

}
