// bc_full_face_exchange.d
//
// Cell-to-cell exchange boundary condition.
// Peter J. 2014-07-26

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

    override void apply_convective(double t)
    {
	throw new Error("Not implemented yet.");
    } // end apply_convective

    override void apply_viscous(double t)
    {
	throw new Error("Not implemented yet.");
    }  // end apply_viscous

} // end class FullFaceExchangeBC
