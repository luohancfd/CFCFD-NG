// bc_slip_wall.hh

#include "bc.hh"

class SlipWallBC : public BoundaryCondition {
public:
    SlipWallBC( Block &bdp, int which_boundary );
    SlipWallBC( const SlipWallBC &bc );
    virtual ~SlipWallBC();
    // default apply_inviscid() is just to reflect normal velocity
    // default apply_viscous() (does nothing)
};
