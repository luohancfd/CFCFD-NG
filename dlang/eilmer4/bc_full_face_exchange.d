// bc_full_face_exchange.d
//
// Cell-to-cell exchange boundary condition.
// Peter J. 2014-07-26

import std.conv;

import bc;

class FullFaceExchangeBC: BoundaryCondition {
public:
    int neighbour_block;
    int neighbour_face;
    int neighbour_orientation;

    this()
    {
	type_code = BCCode.full_face_exchange;
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
	throw new Error("Not implemented yet.");
    } // end apply_convective

    override void apply_viscous(double t)
    {
	throw new Error("Not implemented yet.");
    }  // end apply_viscous

} // end class FullFaceExchangeBC
