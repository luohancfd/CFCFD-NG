// user_defined_effects.d
//
// Authors: RG & PJ
// Date: 2015-03-14

import std.string;
import std.stdio;
import luad.all;
import util.lua_service;

import simcore;
import fvcore;
import fvcell;
import fvinterface;
import globalconfig;
import globaldata;
import ghost_cell_effect;

class UserDefinedGhostCell : GhostCellEffect {
public:
    string luafname;
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
	FVCell ghostCell0, ghostCell1;
	FVInterface IFace;
	auto blk = allBlocks[blk_id];

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k)  {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    // ghostCell0 is closest to domain
		    // ghostCell1 is one layer out.
		    ghostCell0 = blk.get_cell(i,j+1,k);
		    ghostCell1 = blk.get_cell(i,j+2,k);
		    IFace = ghostCell0.iface[Face.south];
		    callGhostCellUDF(t, tLevel, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end i loop
	    } // end k loop
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <- blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i+i,j,k);
		    ghostCell1 = blk.get_cell(i+2,j,k);
		    IFace = ghostCell0.iface[Face.west];
		    callGhostCellUDF(t, tLevel, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end k loop
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i=blk.imin; i <= blk.imax; ++i) {
		    ghostCell0 = blk.get_cell(i,j-1,k);
		    ghostCell1 = blk.get_cell(i,j-2,k);
		    IFace = ghostCell0.iface[Face.north];
		    callGhostCellUDF(t, tLevel, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end i loop
	    } // end j loop
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i-1,j,k);
		    ghostCell1 = blk.get_cell(i-2,j,k);
		    IFace = ghostCell0.iface[Face.east];
		    callGhostCellUDF(t, tLevel, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end k loop
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i= blk.imin; i <= blk.imax; ++i) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i,j,k+1);
		    ghostCell1 = blk.get_cell(i,j,k+2);
		    IFace = ghostCell0.iface[Face.bottom];
		    callGhostCellUDF(t, tLevel, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end i loop
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i,j,k-1);
		    ghostCell1 = blk.get_cell(i,j,k-2);
		    IFace = ghostCell0.iface[Face.top];
		    callGhostCellUDF(t, tLevel, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end i loop
	    break;
	} // end switch which boundary
    }
			
private:
    LuaState _lua;

    void putFlowStateIntoGhostCell(LuaTable t, FVCell ghostCell)
    {
	try {
	    ghostCell.fs.gas.p = t.get!double("p");
	    getArray!double(t.get!LuaTable("T"), ghostCell.fs.gas.T, "T");
	    getArray!double(t.get!LuaTable("massf"), ghostCell.fs.gas.massf, "massf");
	}
	catch (Exception e) {
	    string errMsg = "There was an error trying to read p, T or massf in user-supplied table.\n";
	    errMsg ~= "The error message from the lua state follows.\n";
	    errMsg ~= e.toString();
	    throw new Exception(errMsg);
	}
	GlobalConfig.gmodel.update_thermo_from_pT(ghostCell.fs.gas);
	ghostCell.fs.vel.refx = getDouble(t, "velx", 0.0);
	ghostCell.fs.vel.refy = getDouble(t, "vely", 0.0);
	ghostCell.fs.vel.refz = getDouble(t, "velz", 0.0);
	ghostCell.fs.tke = getDouble(t, "tke", 0.0);
	ghostCell.fs.omega = getDouble(t, "omega", 0.0);
    }

    void callGhostCellUDF(double t, int tLevel, size_t i, size_t j, size_t k,
			  in FVInterface IFace, FVCell ghostCell0, FVCell ghostCell1)
    {
	// 1. Set useful values for caller in table
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
	
	// 2. Call LuaFunction and expect two tables of ghost cell flow state
	auto f = _lua.get!LuaFunction("ghostCells");
	LuaObject[] ret = f(args);
	if ( ret.length < 2 ) {
	    string errMsg = "ERROR: There was a problem in the call to the user-defined ghost cell boundary condition.\n";
	    errMsg ~= format("ERROR: This occured for block [%d] on the %s boundary.", blk_id, face_name[which_boundary]);
	    errMsg ~= "ERROR: Two tables of flow state for the ghost cells are expected\n";
	    errMsg ~= "ERROR: but were not received.\n";
	    throw new Exception(errMsg);
	}
	// 3. Grab Flowstate data from table and populate ghost cell
	putFlowStateIntoGhostCell(ret[0].to!LuaTable(), ghostCell0);
	putFlowStateIntoGhostCell(ret[1].to!LuaTable(), ghostCell1);
    }
}
