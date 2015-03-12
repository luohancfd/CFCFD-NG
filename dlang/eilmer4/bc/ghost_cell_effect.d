// ghost_cell_effect.d
//
// RG & PJ 2015-12-03 : first hack

import std.json;
import std.string;

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

@nogc
void reflect_normal_velocity(ref FVCell cell, in FVInterface IFace)
// Reflects normal velocity with respect to the cell interface.
//
// The process is to rotate the velocity vector into the local frame of
// the interface, negate the normal (local x-component) velocity and
// rotate back to the global frame.
{
    cell.fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    cell.fs.vel.refx = -(cell.fs.vel.x);
    cell.fs.vel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

@nogc
void reflect_normal_magnetic_field(ref FVCell cell, in FVInterface IFace)
{
    cell.fs.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    cell.fs.B.refx = -(cell.fs.B.x);
    cell.fs.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
}

GhostCellEffect make_GCE_from_json(JSONValue jsonData, int blk_id, int boundary)
{
    string gceType = jsonData["type"].str;
    GhostCellEffect newGCE;
    switch (gceType) {
    case "internal_copy_then_reflect":
	newGCE = new GhostCellInternalCopyThenReflect(blk_id, boundary);
	break;
    case "flowstate_copy":
	auto flowstate = new FlowState(jsonData["flowstate"], GlobalConfig.gmodel);
	newGCE = new GhostCellFlowStateCopy(blk_id, boundary, flowstate);
	break;
    case "extrapolate_copy":
	int xOrder = getJSONint(jsonData, "x_order", 0);
	newGCE = new GhostCellExtrapolateCopy(blk_id, boundary, xOrder);
	break;
    case "full_face_exchange_copy":
	int otherBlock = getJSONint(jsonData, "other_block", -1);
	string otherFaceName = getJSONstring(jsonData, "other_face", "none");
	
	int neighbourOrientation = getJSONint(jsonData, "neighbour_orientation", 0);
	newGCE = new GhostCellFullFaceExchangeCopy(blk_id, boundary, 
						   otherBlock, face_index(otherFaceName), neighbourOrientation);
	break;
    default:
	string errMsg = format("ERROR: The GhostCellEffect type: '%s' is unknown.", gceType);
	throw new Exception(errMsg);
    }
    return newGCE;
}

class GhostCellEffect {
public:
    int blk_id;
    int which_boundary;

    this(int id, int boundary)
    {
	blk_id = id;
	which_boundary = boundary;
    }
    
    abstract void apply(double t);
}

class GhostCellInternalCopyThenReflect : GhostCellEffect {
public:
    
    this(int id, int boundary)
    {
	super(id, boundary);
    }

    override void apply(double t)
    {
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVInterface IFace;
	auto copy_opt = CopyDataOption.minimal_flow;
	auto blk = allBlocks[blk_id];

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.north];
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j-1,k);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.east];
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i-1,j,k);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.south];
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j+1,k);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.west];
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i+1,j,k);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.top];
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j,k-1);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    // ghost cell 1.
		    src_cell = blk.get_cell(i,j,k);
		    IFace = src_cell.iface[Face.bottom];
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		    // ghost cell 2.
		    src_cell = blk.get_cell(i,j,k+1);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.copy_values_from(src_cell, copy_opt);
		    reflect_normal_velocity(dest_cell, IFace);
		    if (GlobalConfig.MHD) {
			reflect_normal_magnetic_field(dest_cell, IFace);
		    }
		} // end j loop
	    } // for i
	    break;
	} // end switch which_boundary
    }
}

class GhostCellFlowStateCopy : GhostCellEffect {
public:
    FlowState fstate;

    this(int id, int boundary, in FlowState _fstate)
    {
	super(id, boundary);
	fstate = new FlowState(_fstate);
    }

    override void apply(double t)
    {
    	// Fill ghost cells with data from just inside the boundary
	// using zero-order extrapolation (i.e. just copy the data).
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVInterface dest_face;
	auto blk = allBlocks[blk_id];

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply()

}

class GhostCellExtrapolateCopy : GhostCellEffect {
public:
    int xOrder;
    
    this(int id, int boundary, int x_order)
    {
	super(id, boundary);
	xOrder = x_order;
    }
    
    override void apply(double t)
    {
	// Fill ghost cells with data from just inside the boundary
	// using zero-order extrapolation (i.e. just copy the data).
	// We assume that this boundary is an outflow boundary.
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	FVCell cell_1, cell_2;
	auto gmodel = GlobalConfig.gmodel;
	size_t nsp = gmodel.n_species;
	size_t nmodes = gmodel.n_modes;
	auto blk = allBlocks[blk_id];

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (j-1)        (j)           (j+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j-1,k);
			dest_cell = blk.get_cell(i,j+1,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (j)        (j+1)       (j+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j+2,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    }
		    else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j+1,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j+2,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    } 
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (i-1)        (i)           (i+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i-1,j,k);
			dest_cell = blk.get_cell(i+1,j,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (i)        (i+1)       (i+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i+2,j,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    }
		    else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i+1,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i+2,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (j+1)        (j)           (j-1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j+1,k);
			dest_cell = blk.get_cell(i,j-1,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (j)        (j-1)       (j-2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j-2,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j-1,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j-2,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  ---[ghost cell 2]---|--- [dest] ---|||--- [1] ---|---[2]----
			//      (i-2)                 (i-1)           (i)       (i+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i+1,j,k);
			dest_cell = blk.get_cell(i-1,j,k);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[dest]---|---[1]---|||---[2]---|------|
			//       (i-2)       (i-1)       (i)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i-2,j,k);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i-1,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i-2,j,k);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (k-1)        (k)           (k+1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j,k-1);
			dest_cell = blk.get_cell(i,j,k+1);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (k)        (k+1)       (k+2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j,k+2);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			// Zero-order extrapolation
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j,k+1);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j,k+2);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    if ( xOrder == 1 ) {
			//  |--- [2] ---|--- [1] ---|||--- [dest] ---|---[ghost cell 2]----
			//      (k+1)        (k)           (k-1)
			//  dest: ghost cell 1
			//  [1]: first interior cell
			//  [2]: second interior cell
			// This extrapolation assumes that cell-spacing between
			// cells 1 and 2 continues on in the exterior
			cell_1 = blk.get_cell(i,j,k);
			cell_2 = blk.get_cell(i,j,k+2);
			dest_cell = blk.get_cell(i,j,k-1);
			// Extrapolate on primitive variables
			// 1. First exterior point
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
			// 2. Second exterior point
			//  |---[2]---|||---[1]---|---[dest]------
			//      (k)        (k-1)       (k-2)
			cell_2 = cell_1;
			cell_1 = dest_cell;
			dest_cell = blk.get_cell(i,j,k-2);
			dest_cell.fs.gas.rho = 2.0*cell_1.fs.gas.rho - cell_2.fs.gas.rho;
			for ( size_t imode = 0; imode < nmodes; ++imode ) {
			    dest_cell.fs.gas.e[imode] = 2.0*cell_1.fs.gas.e[imode] - cell_2.fs.gas.e[imode];
			}
			if ( nsp > 1 ) {
			    for ( size_t isp = 0; isp < nsp; ++isp ) {
				dest_cell.fs.gas.massf[isp] = 2.0*cell_1.fs.gas.massf[isp] - cell_2.fs.gas.massf[isp];
			    }
			    scale_mass_fractions(dest_cell.fs.gas.massf);
			}
			else {
			    dest_cell.fs.gas.massf[0] = 1.0;
			}
			gmodel.update_thermo_from_rhoe(dest_cell.fs.gas);
			dest_cell.fs.vel.refx = 2.0*cell_1.fs.vel.x - cell_2.fs.vel.x;
			dest_cell.fs.vel.refy = 2.0*cell_1.fs.vel.y - cell_2.fs.vel.y;
			dest_cell.fs.vel.refz = 2.0*cell_1.fs.vel.z - cell_2.fs.vel.z;
			dest_cell.fs.tke = 2.0*cell_1.fs.tke - cell_2.fs.tke;
			dest_cell.fs.omega = 2.0*cell_1.fs.omega - cell_2.fs.omega;
			dest_cell.fs.mu_t = 2.0*cell_1.fs.mu_t - cell_2.fs.mu_t;
			dest_cell.fs.k_t = 2.0*cell_1.fs.k_t - cell_2.fs.k_t;
		    } else {
			src_cell = blk.get_cell(i,j,k);
			dest_cell = blk.get_cell(i,j,k-1);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
			dest_cell = blk.get_cell(i,j,k-2);
			dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    }
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply()
}

class GhostCellFullFaceExchangeCopy : GhostCellEffect {
public:
    int neighbourBlock;
    int neighbourFace;
    int neighbourOrientation;

    this(int id, int boundary,
	 int otherBlock, int otherFace, int orient)
    {
	super(id, boundary);
	neighbourBlock = otherBlock;
	neighbourFace = otherFace;
	neighbourOrientation = orient;
    }

    override void apply(double t)
    {
	// TODO Check me!  This is a work-around.
	// We should be able to directly reference the BCs block as blk.
	auto blk = allBlocks[blk_id];
	blk.copy_into_ghost_cells(which_boundary, 
				  allBlocks[neighbourBlock],
				  neighbourFace, neighbourOrientation,
				  CopyDataOption.all, true);
    }
}
