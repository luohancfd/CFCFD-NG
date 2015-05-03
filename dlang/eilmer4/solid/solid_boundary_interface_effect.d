module solid_boundary_interface_effect;

import std.stdio;
import std.json;
import std.string;
import std.format;
import luad.all;
import luad.c.lua;


import simcore;
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
    case "user_defined":
	string fname = getJSONstring(jsonData, "filename", "none");
	newBIE = new SolidBIE_UserDefined(blk_id, boundary, fname);
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
    string type;
    
    this(int id, int boundary, string _type) {
	blkId = id;
	whichBoundary = boundary;
	type = _type;
    }

    void apply(double t, int tLevel) {}
}

class SolidBIE_FixedT : SolidBoundaryInterfaceEffect {
public:
    this(int id, int boundary, double Twall)
    {
	super(id, boundary, "FixedT");
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

class SolidBIE_UserDefined : SolidBoundaryInterfaceEffect {
public:
    string luafname;
private:
    LuaState _lua;

public:
    this(int id, int boundary, string fname)
    {
	super(id, boundary, "UserDefined");
	luafname = fname;
    }
    
    void setLuaState(LuaState lua)
    {
	_lua = lua;
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
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace);
		}
	    }
	    break;
	case Face.east:
	    i = blk.imax + 1;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace);
		}
	    }
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    IFace = blk.getIfj(i, j, k);
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace);
		}
	    }
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    IFace = blk.getIfi(i, j, k);
		    callSolidIfaceUDF(t, tLevel, i, j, k, IFace);
		}
	    }
	    break;
	case Face.top:
	    throw new Error("[TODO] FixedT bc not implemented for TOP face.");
	case Face.bottom:
	    throw new Error("[TODO] FixedT bc not implemented for BOTTOM face.");

	}
    }

    void callSolidIfaceUDF(double t, int tLevel, size_t i, size_t j, size_t k,
			   SolidFVInterface IFace)
    {
	// Set some userful values for the caller in table
	auto args = _lua.newTable(0, 19);
	args["t"] = t;
	args["dt"] = dt_global;
	args["timeStep"] = step;
	args["timeLevel"] = tLevel;
	args["x"] = IFace.pos.x;
	args["y"] = IFace.pos.y;
	args["z"] = IFace.pos.z;
	args["csX"] = IFace.n.x;
	args["csY"] = IFace.n.y;
	args["csZ"] = IFace.n.z;
	args["csX1"] = IFace.t1.x;
	args["csY1"] = IFace.t1.y;
	args["csZ1"] = IFace.t1.z;
	args["csX2"] = IFace.t2.x;
	args["csY2"] = IFace.t2.y;
	args["csZ2"] = IFace.t2.z;
	args["i"] = i;
	args["j"] = j;
	args["k"] = k;

	// Call LuaFunction and expect back a temperature value.
	auto f = _lua.get!LuaFunction("solidInterface");
	LuaObject[] ret = f(args);
	if ( ret.length < 1 ) {
	    string errMsg = "ERROR: There was a problem in the call to the user-defined solid interface boundary condition.\n";
	    errMsg ~= format("ERROR: This occurred for solidblock [%d] on the %s boundary.", blkId, face_name[whichBoundary]);
	    errMsg ~= "ERROR: A single float value for temperature was expected, but not received.\n";
	    throw new Exception(errMsg);
	}

	// Grab temperature and set interface with that value
	IFace.T = ret[0].to!double();
    }
}
