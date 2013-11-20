// bc_moving_wall.hh

#include "bc.hh"

class MovingWallBC : public BoundaryCondition {
public:
    double r_omega;
public:
    MovingWallBC(Block *bdp, int which_boundary, double r_omega, double emissivity);
    MovingWallBC(const MovingWallBC &bc);
    MovingWallBC();
    MovingWallBC & operator=(const MovingWallBC &bc);
    virtual ~MovingWallBC();
    virtual void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    virtual int apply_viscous(double t); 
};
