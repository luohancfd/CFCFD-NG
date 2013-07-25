// bc_adjacent.hh

#include "bc.hh"

class AdjacentBC : public BoundaryCondition {
public:
    AdjacentBC( Block *bdp, int which_boundary, int other_block, 
		int other_face, int neighbour_orientation=0 );
    AdjacentBC( const AdjacentBC &bc );
    AdjacentBC();
    AdjacentBC & operator=(const AdjacentBC &bc);
    virtual ~AdjacentBC();
    virtual int apply_convective( double t ); // does nothing
    // default apply_viscous() (does nothing)
};
