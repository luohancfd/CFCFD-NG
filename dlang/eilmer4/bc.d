// bc.d
// Base class for boundary condition objects, for use in Eilmer4
//
// Peter J. 2014-07-20 first cut.

enum
    SLIP_WALL = 0;

class BoundaryCondition {
    int type_code;
    bool is_wall() 
    {
	return true; // [TODO] port C++ version
    }

} // end class BoundaryCondition
