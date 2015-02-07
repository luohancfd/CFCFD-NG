// bc_full_face_exchange.d
//
// Cell-to-cell exchange boundary condition.
// Peter J. 2014-07-26

import std.conv;

import bc;
import sblock;

class FullFaceExchangeBC: BoundaryCondition {
public:
    int neighbour_block;
    int neighbour_face;
    int neighbour_orientation;

    this(ref SBlock blk_, int which_boundary_, 
	 int other_block, int other_face, int orient)
    {
	blk = blk_;
	which_boundary = which_boundary_;
	type_code = BCCode.full_face_exchange;
	is_wall = false;
	neighbour_block = other_block;
	neighbour_face = other_face;
	neighbour_orientation = orient;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FullFaceExchangeBC(";
	repr ~= "neighbour_block=" ~ to!string(neighbour_block);
	repr ~= ", neighbour_face=" ~ to!string(neighbour_face);
	repr ~= ", neighbour_orientation=" ~ to!string(neighbour_orientation);
	repr ~= ")";
	return to!string(repr);
    }

    override void apply_convective(double t)
    {
	throw new Error("TODO Not implemented yet.");
    } // end apply_convective

    override void apply_viscous(double t)
    {
	throw new Error("TODO Not implemented yet.");
    }  // end apply_viscous

} // end class FullFaceExchangeBC
