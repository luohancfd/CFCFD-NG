// boundary_interface_effect.d
//
// Effects needed to compute viscous fluxes and the like.
//
// PJ and RG, 2015-04-28, initial code mainly from 
//    the break-up of the Fixed_T boundary condition.
//

module boundary_interface_effect;

import std.json;
import std.string;
import std.conv;
import std.stdio;

import geom;
import json_helper;
import globalconfig;
import globaldata;
import flowstate;
import fvcore;
import fvinterface;
import fvcell;
import block;
import sblock;
import gas;
import user_defined_effects;


BoundaryInterfaceEffect make_BIE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string bieType = jsonData["type"].str;
    // If we need access to a gas model in here, 
    // be sure to use GlobalConfig.gmodel_master.
    BoundaryInterfaceEffect newBIE;
    switch (bieType) {
    case "copy_cell_data":
	newBIE = new BIE_CopyCellData(blk_id, boundary);
	break;
    case "zero_velocity":
	newBIE = new BIE_ZeroVelocity(blk_id, boundary);
	break;
    case "fixed_temperature":
	double Twall = getJSONdouble(jsonData, "Twall", 300.0);
	newBIE = new BIE_FixedT(blk_id, boundary, Twall);
	break;
    case "update_thermo_trans_coeffs":
	newBIE = new BIE_UpdateThermoTransCoeffs(blk_id, boundary);
	break;
    case "wall_k_omega":
	newBIE = new BIE_WallKOmega(blk_id, boundary);
	break;
    case "user_defined":
     	string fname = getJSONstring(jsonData, "filename", "none");
	newBIE = new BIE_UserDefined(blk_id, boundary, fname);
	break;
    default:
	string errMsg = format("ERROR: The BoundaryInterfaceEffect type: '%s' is unknown.", bieType);
	throw new Exception(errMsg);
    }
    return newBIE;
}


class BoundaryInterfaceEffect {
public:
    int blk_id;
    int which_boundary;
    string type;

    this(int id, int boundary, string _type)
    {
	blk_id = id;
	which_boundary = boundary;
	type = _type;
    }
    override string toString() const
    {
	return "BoundaryInterfaceEffect()";
    }
    abstract void apply(double t, int gtl, int ftl);
} // end class BoundaryInterfaceEffect


class BIE_CopyCellData : BoundaryInterfaceEffect {
    this(int id, int boundary, double Twall=300.0)
    {
	super(id, boundary, "CopyCellData");
    }

    override string toString() const 
    {
	return "CopyCellData()";
    }

    override void apply(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	// [TODO] [FIXME] PJ
	// writeln("BoundaryInterfaceEffect.apply(): blk_id=", blk_id, " length=", allBlocks.length);
	// The following is going to be a problem, when running in parallel 
	// allBlocks is not shared, so not all threads see the filled array.
	// Some see a zero-length array.
	auto blk = allBlocks[blk_id];
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.copy_values_from(cell.fs);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_CopyCellData


class BIE_ZeroVelocity : BoundaryInterfaceEffect {
    this(int id, int boundary)
    {
	super(id, boundary, "ZeroVelocity");
    }

    override string toString() const 
    {
	return "ZeroVelocity()";
    }

    override void apply(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto blk = allBlocks[blk_id];
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.vel.refx = 0.0; fs.vel.refy = 0.0; fs.vel.refz = 0.0;
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_ZeroVelocity


class BIE_FixedT : BoundaryInterfaceEffect {
public:
    double Twall;

    this(int id, int boundary, double Twall)
    {
	super(id, boundary, "FixedT");
	this.Twall = Twall;
    }

    override string toString() const 
    {
	return "FixedT(Twall=" ~ to!string(Twall) ~ ")";
    }

    override void apply(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto blk = allBlocks[blk_id];
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    foreach(ref elem; fs.gas.T) elem = Twall;
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    foreach(ref elem; fs.gas.T) elem = Twall;
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    foreach(ref elem; fs.gas.T) elem = Twall;
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    foreach(ref elem; fs.gas.T) elem = Twall;
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    foreach(ref elem; fs.gas.T) elem = Twall;
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    foreach(ref elem; fs.gas.T) elem = Twall;
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_FixedT


class BIE_UpdateThermoTransCoeffs : BoundaryInterfaceEffect {
    this(int id, int boundary, double Twall=300.0)
    {
	super(id, boundary, "UpdateThermoTransCoeffs");
    }

    override string toString() const 
    {
	return "UpdateThermoTransCoeffs()";
    }

    override void apply(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto blk = allBlocks[blk_id];
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    gmodel.update_thermo_from_pT(fs.gas);
		    gmodel.update_trans_coeffs(fs.gas);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()
} // end class BIE_UpdateThermoTransCoeffs


class BIE_WallKOmega : BoundaryInterfaceEffect {
    // Menter's slightly-rough-surface boundary condition is described
    // in Wilcox's 2006 text, eqn 7.36.
    // For low-resolution grids, the k-omega model is reported to over-estimate
    // the magnitude of omega, well out into the boundary layer so,
    // to get reasonable values for omega close to the wall, we propagate
    // the 1/y**2 form of the omega data out a few cells from the wall.
    this(int id, int boundary)
    {
	super(id, boundary, "WallKOmega");
    }

    override string toString() const 
    {
	return "WallKOmega()";
    }

    override void apply(double t, int gtl, int ftl)
    {
	size_t i, j, k;
	FVCell cell;
	FVInterface IFace;
	auto blk = allBlocks[blk_id];
	auto gmodel = blk.myConfig.gmodel;

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.north];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.east];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.south];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end i loop
	    } // end for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.west];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.top];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    cell = blk.get_cell(i,j,k);
		    IFace = cell.iface[Face.bottom];
		    FlowState fs = IFace.fs;
		    fs.tke = 0.0;
		    fs.omega = ideal_omega_at_wall(cell);
		} // end j loop
	    } // end for i
	    break;
	} // end switch which_boundary
    } // end apply()

    @nogc
    double ideal_omega_at_wall(in FVCell cell)
    {
	auto wall_gas = cell.cell_at_nearest_wall.fs.gas;
	double d0 = cell.half_cell_width_at_wall;
	return 400.0 * wall_gas.mu / wall_gas.rho / (d0 * d0);
    }

    @nogc
    double ideal_omega(in FVCell cell)
    {
	double d0 = cell.half_cell_width_at_wall;
	double d = cell.distance_to_nearest_wall;
	return ideal_omega_at_wall(cell) * (d0 * d0) / ((d0 + d) * (d0 + d));
    }
} // end class BIE_WallKOmega
