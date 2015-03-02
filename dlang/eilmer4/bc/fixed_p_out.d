// bc/fixed_p_out.d
//
// Outflow boundary condition that simply extrapolates the flow near the boundary
// into the ghost cells just outside the boundary and then sets the pressure as
// specified.
//
// [TODO] implement the high-order reconstruction, as per ExtrapolateOutBC.
// Peter J. 2014-07-26

import std.conv;

import gas;
import fvcore;
import fvinterface;
import fvcell;
import bc;
import block;
import sblock;
import globalconfig;
import globaldata;

class FixedPOutBC: BoundaryCondition {
public:
    double Pout;
    double Tout;
    bool use_Tout = false;
    int x_order = 0; // default to lowest order

    this(int id, int boundary, double Pout, double Tout,
	 bool use_Tout=false, int x_order=0) 
    {
	blk_id = id;
	which_boundary = boundary;
	type_code = BCCode.fixed_p_out;
	is_wall = false;
	ghost_cell_data_available = true;
	sets_conv_flux_directly = false;
	sets_visc_flux_directly = false;
	emissivity = 0.0;
	this.Pout = Pout;
	this.Tout = Tout;
	this.use_Tout = use_Tout;
	this.x_order = x_order;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FixedPOutBC(";
	repr ~= "Pout=" ~ to!string(Pout);
	repr ~= ", Tout=" ~ to!string(Tout);
	repr ~= ", use_Tout=" ~ to!string(use_Tout);
	repr ~= ", x_order=" ~ to!string(x_order);
	repr ~= ")";
	return to!string(repr);
    }

    override void apply_convective(double t)
    {
	// Fill ghost cells with data from just inside the boundary
	// using zero-order extrapolation (i.e. just copy the data).
	// We assume that this boundary is an outflow boundary.
	size_t i, j, k;
	FVCell src_cell, dest_cell;
	auto gmodel = GlobalConfig.gmodel;
	auto blk = allBlocks[blk_id];

	final switch (which_boundary) {
	case Face.north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end i loop
	    } // for k
	    break;
	case Face.east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for k
	    break;
	case Face.south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end i loop
	    } // for k
	    break;
	case Face.west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for k
	    break;
	case Face.top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for i
	    break;
	case Face.bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    src_cell = blk.get_cell(i,j,k);
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.copy_values_from(src_cell, CopyDataOption.minimal_flow);
		    dest_cell.fs.gas.p = Pout;
		    if ( use_Tout ) foreach(ref elem; dest_cell.fs.gas.T) elem = Tout; 
		    gmodel.update_thermo_from_pT(dest_cell.fs.gas);
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply_convective()

    // Let the base class implementations do the work.
    // apply_viscous

} // end class FixedPOutBC