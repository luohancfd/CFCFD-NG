/**
 * solidbc.d
 *
 * Author: Rowan G. and Peter J.
 */

module solidbc;

import std.json;

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

    SolidBoundaryInterfaceEffect[] preSpatialDerivAction;
    SolidBoundaryFluxEffect[] postFluxAction;

}
