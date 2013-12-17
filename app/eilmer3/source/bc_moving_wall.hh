// bc_moving_wall.hh

#include "../../../lib/geometry2/source/geom.hh"
#include "bc.hh"

class MovingWallBC : public BoundaryCondition {
public:
    Vector3 r_omega;  // angular velocity of rotating surface
    Vector3 centre;   // a point on the axis of rotation
    Vector3 v_trans;  // translational velocity to superimpose
public:
    MovingWallBC(Block *bdp, int which_boundary, Vector3 omega,
		 Vector3 centre, Vector3 v_trans, double emissivity);
    MovingWallBC(const MovingWallBC &bc);
    MovingWallBC();
    MovingWallBC & operator=(const MovingWallBC &bc);
    virtual ~MovingWallBC();
    virtual void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    virtual int apply_viscous(double t); 
};
