module solid_boundary_interface_effect;

import std.stdio;
import std.json;
import std.string;
import std.format;

import json_helper;
import geom;
import globaldata;
import solidfvinterface;

SolidBoundaryInterfaceEffect makeSolidBIEfromJson(JSONValue jsonData, int blk_id, int boundary)
{
    string bieType = jsonData["type"].str;
    SolidBoundaryInterfaceEffect newBIE;
    switch (bieType) {
    case "fixed_temperature":
	double Twall = getJSONdouble(jsonData, "Twall", 300.0);
	newBIE = new SolidBIE_FixedT(blk_id, boundary, Twall);
	break;
    default:
	string errMsg = format("ERROR: The SolidBoundaryInterfaceEffect type: '%s' is unknown.", bieType);
	throw new Exception(errMsg);
    }
    return newBIE;
}

class SolidBoundaryInterfaceEffect {
public:
    int blkId;
    int whichBoundary;
    
    this(int id, int boundary) {
	blkId = id;
	whichBoundary = boundary;
    }

    void apply(double t, int tLevel) {}
}

class SolidBIE_FixedT : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary, double Twall)
    {
	super(id, boundary);
	_Twall = Twall;
    }

    override void apply(double t, int tLevel)
    {
	size_t i, j, k;
	SolidFVInterface IFace;
	auto blk = allSolidBlocks[blkId];

	final switch (whichBoundary) {
	case Face.north:
	    j = blk.jmax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.east:
	    i = blk.imax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    IFace.T = _Twall;
		}
	    }
	    break;
	case Face.top:
	    throw new Error("[TODO] FixedT bc not implemented for TOP face.");
	case Face.bottom:
	    throw new Error("[TODO] FixedT bc not implemented for BOTTOM face.");

	}

    }

private:
    double _Twall;
}
