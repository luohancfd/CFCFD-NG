// bc_jump_wall.hh

#include "bc.hh"

class JumpWallBC : public BoundaryCondition {
public:
    double Twall; // specified wall temperature
    double sigma_jump; // accommodation coefficient
public:
    JumpWallBC(Block *bdp, int which_boundary, double Twall, double sigma_jump);
    JumpWallBC(const JumpWallBC &bc);
    JumpWallBC();
    JumpWallBC & operator=(const JumpWallBC &bc);
    virtual ~JumpWallBC();
    virtual void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    virtual int apply_viscous(double t); // sets wall temperature and jump velocity
};
