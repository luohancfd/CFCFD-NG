// user_defined_effects.d
//
// Authors: RG & PJ
// Date: 2015-03-14

import std.string;
import std.stdio;
import luad.all;
import luad.c.lua;
import util.lua_service;

import geom;
import simcore;
import flowstate;
import fvcore;
import fvcell;
import fvinterface;
import globalconfig;
import globaldata;
import ghost_cell_effect;
import boundary_interface_effect;
import luaflowstate;
import lua_helper;

class UserDefinedGhostCell : GhostCellEffect {
public:
    string luafname;
    this(int id, int boundary, string fname)
    {
	super(id, boundary, "UserDefined");
	luafname = fname;
    }
    override string toString() const
    {
	return "UserDefinedGhostCellEffect(fname=" ~ luafname ~ ")";
    }

    void setLuaState(LuaState lua)
    {
	_lua = lua;
    }

    override void apply(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell ghostCell0, ghostCell1;
	FVInterface IFace;

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
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end i loop
	    } // end k loop
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    ghostCell0 = blk.get_cell(i+1,j,k);
		    ghostCell1 = blk.get_cell(i+2,j,k);
		    IFace = ghostCell0.iface[Face.west];
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
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
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
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
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
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
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
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
		    callGhostCellUDF(t, gtl, ftl, i, j, k, IFace, ghostCell0, ghostCell1);
		} // end j loop
	    } // end i loop
	    break;
	} // end switch which boundary
    }
			
private:
    LuaState _lua;

    void putFlowStateIntoGhostCell(LuaTable t, FVCell ghostCell)
    {
	auto gmodel = blk.myConfig.gmodel;
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
	gmodel.update_thermo_from_pT(ghostCell.fs.gas);
	gmodel.update_sound_speed(ghostCell.fs.gas);
	ghostCell.fs.vel.refx = getDouble(t, "velx", 0.0);
	ghostCell.fs.vel.refy = getDouble(t, "vely", 0.0);
	ghostCell.fs.vel.refz = getDouble(t, "velz", 0.0);
	ghostCell.fs.tke = getDouble(t, "tke", 0.0);
	ghostCell.fs.omega = getDouble(t, "omega", 0.0);
    }

    void callGhostCellUDF(double t, int gtl, int ftl, size_t i, size_t j, size_t k,
			  in FVInterface IFace, FVCell ghostCell0, FVCell ghostCell1)
    {
	// 1. Set useful values for caller in table
	auto args = _lua.newTable(0, 20);
	args["t"] = t; 
	args["dt"] = dt_global;
	args["timeStep"] = step;
	args["gridTimeLevel"] = gtl;
	args["flowTimeLevel"] = ftl;
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
	    errMsg ~= format("ERROR: This occurred for block [%d] on the %s boundary.", blk.id, face_name[which_boundary]);
	    errMsg ~= "ERROR: Two tables of flow state for the ghost cells are expected\n";
	    errMsg ~= "ERROR: but were not received.\n";
	    throw new Exception(errMsg);
	}
	// 3. Grab Flowstate data from table and populate ghost cell
	putFlowStateIntoGhostCell(ret[0].to!LuaTable(), ghostCell0);
	putFlowStateIntoGhostCell(ret[1].to!LuaTable(), ghostCell1);
    }
}

class BIE_UserDefined : BoundaryInterfaceEffect {
public:
    string luafname;
    this(int id, int boundary, string fname)
    {
	super(id, boundary, "UserDefined");
	luafname = fname;
    }

    override string toString() const
    {
	return "UserDefined(fname=" ~ luafname ~ ")";
    }

    void setLuaState(LuaState lua)
    {
	_lua = lua;
    }

    override void apply(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k)  {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end i loop
	    } // end k loop
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end k loop
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i=blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end i loop
	    } // end j loop
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end k loop
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i= blk.imin; i <= blk.imax; ++i) {
		for (j=blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end i loop
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    callInterfaceUDF(t, gtl, ftl, i, j, k, IFace);
		} // end j loop
	    } // end i loop
	    break;
	} // end switch which boundary
    }
private:
    LuaState _lua;

    void putFlowStateIntoInterface(LuaTable t, FVInterface iface)
    {
	// [TODO] It would be more elegant to iterate over
	// the key-value pairs. We'll do this when we move
	// back to the pure C LUA API.

	// Now the user might only set some of the flowstate.
	// So we need to test every possibility and only set
	// the non-nil values.
	FlowState fs = iface.fs;

	if ( !t.get!LuaObject("p").isNil ) {
	    fs.gas.p = t.get!double("p");
	}
	if ( !t.get!LuaObject("T").isNil ) {
	    // Temperature should be provided as an array.
	    getArray!double(t.get!LuaTable("T"), fs.gas.T, "T");
	}
	if ( !t.get!LuaObject("massf").isNil ) {
	    // mass fractions should be provided as an array
	    getArray!double(t.get!LuaTable("massf"), fs.gas.massf, "massf");
	}
	if ( !t.get!LuaObject("velx").isNil ) {
	    fs.vel.refx = t.get!double("velx");
	}
	if ( !t.get!LuaObject("vely").isNil ) {
	    fs.vel.refy = t.get!double("vely");
	}
	if ( !t.get!LuaObject("velz").isNil ) {
	    fs.vel.refz = t.get!double("velz");
	}
	if ( !t.get!LuaObject("tke").isNil ) {
	    fs.tke = t.get!double("tke");
	}
	if ( !t.get!LuaObject("omega").isNil ) {
	    fs.omega = t.get!double("omega");
	}
    }
	    
    void callInterfaceUDF(double t, int gtl, int ftl, size_t i, size_t j, size_t k,
			  FVInterface IFace)
    {
	// 1. Set useful values for caller in table
	auto args = _lua.newTable(0, 20);
	args["t"] = t; 
	args["dt"] = dt_global;
	args["timeStep"] = step;
	args["gridTimeLevel"] = gtl;
	args["flowTimeLevel"] = ftl;
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
	
	// 2. Call LuaFunction and expect a table of values for flow state
	auto f = _lua.get!LuaFunction("interface");
	LuaObject[] ret = f(args);
	if ( ret.length < 1 ) {
	    string errMsg = "ERROR: There was a problem in the call to the user-defined interface boundary condition.\n";
	    errMsg ~= format("ERROR: This occurred for block [%d] on the %s boundary.", blk.id, face_name[which_boundary]);
	    errMsg ~= "ERROR: One table of flow state for the interface is expected\n";
	    errMsg ~= "ERROR: but zero were not received.\n";
	    throw new Exception(errMsg);
	}
	// 3. Grab Flowstate data from table and populate interface
	putFlowStateIntoInterface(ret[0].to!LuaTable(), IFace);
    }

}
