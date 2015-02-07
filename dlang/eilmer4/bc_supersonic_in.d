// bc_supersonic_in.d
//
// Inflow boundary condition that simply copies the (assumed) supersonic flow 
// into the ghost cells just outside the boundary.
//
// Peter J. 2014-07-26

import std.conv;

import fvcore;
import flowstate;
import fvinterface;
import fvcell;
import bc;
import block;
import sblock;
import globaldata;

class SupersonicInBC: BoundaryCondition {
public:
    int inflow_condition_id = 0;

    this(ref SBlock blk, int which_boundary, int inflow_condition_id=0) 
    {
	type_code = BCCode.supersonic_in;
	is_wall = false;
	this.inflow_condition_id = inflow_condition_id;
	this.which_boundary = which_boundary;
	blk.bc[which_boundary] = this;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "SupersonicInBC(";
	repr ~= "inflow_condition_id=" ~ to!string(inflow_condition_id);
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
	FVInterface dest_face;
	FlowState fstate = myFlowStates[inflow_condition_id];

	final switch ( which_boundary ) {
	case north:
	    j = blk.jmax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    dest_cell = blk.get_cell(i,j+1,k);
		    dest_cell.fs.copy_values_from(fstate);
		    // Although this is principally an inviscid BC,
		    // we need the face values for derivatives.
		    dest_face = blk.get_ifj(i,j+1,k);
		    dest_face.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j+2,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end i loop
	    } // for k
	    break;
	case east:
	    i = blk.imax;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i+1,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		    // Although this is principally an inviscid BC,
		    // we need the face values for derivatives.
		    dest_face = blk.get_ifi(i+1,j,k);
		    dest_face.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i+2,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for k
	    break;
	case south:
	    j = blk.jmin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (i = blk.imin; i <= blk.imax; ++i) {
		    dest_cell = blk.get_cell(i,j-1,k);
		    dest_cell.fs.copy_values_from(fstate);
		    // Although this is principally an inviscid BC,
		    // we need the face values for derivatives.
		    dest_face = blk.get_ifj(i,j,k);
		    dest_face.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j-2,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end i loop
	    } // for k
	    break;
	case west:
	    i = blk.imin;
	    for (k = blk.kmin; k <= blk.kmax; ++k) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i-1,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		    // Although this is principally an inviscid BC,
		    // we need the face values for derivatives.
		    dest_face = blk.get_ifi(i,j,k);
		    dest_face.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i-2,j,k);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for k
	    break;
	case top:
	    k = blk.kmax;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i,j,k+1);
		    dest_cell.fs.copy_values_from(fstate);
		    // Although this is principally an inviscid BC,
		    // we need the face values for derivatives.
		    dest_face = blk.get_ifk(i,j,k+1);
		    dest_face.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j,k+2);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for i
	    break;
	case bottom:
	    k = blk.kmin;
	    for (i = blk.imin; i <= blk.imax; ++i) {
		for (j = blk.jmin; j <= blk.jmax; ++j) {
		    dest_cell = blk.get_cell(i,j,k-1);
		    dest_cell.fs.copy_values_from(fstate);
		    // Although this is principally an inviscid BC,
		    // we need the face values for derivatives.
		    dest_face = blk.get_ifk(i,j,k);
		    dest_face.fs.copy_values_from(fstate);
		    dest_cell = blk.get_cell(i,j,k-2);
		    dest_cell.fs.copy_values_from(fstate);
		} // end j loop
	    } // for i
	    break;
	} // end switch
    } // end apply_convective()

    // Let the base class implementations do the work.
    // apply_viscous

} // end class SupersonicInBC
