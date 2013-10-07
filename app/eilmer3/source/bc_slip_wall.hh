// bc_slip_wall.hh

#include "bc.hh"

class SlipWallBC : public BoundaryCondition {
public:
    SlipWallBC(Block *bdp, int which_boundary, double emissivity);
    SlipWallBC(const SlipWallBC &bc);
    SlipWallBC();
    SlipWallBC & operator=(const SlipWallBC &bc);
    virtual ~SlipWallBC();
    // default apply_convective() is just to reflect normal velocity
    // default apply_viscous() (does nothing)
};
