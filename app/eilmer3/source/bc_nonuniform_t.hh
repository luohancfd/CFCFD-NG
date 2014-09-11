// bc_nonuniform_t.hh

#include "bc.hh"

class NonuniformTBC : public BoundaryCondition {
public:
    std::vector<double> T_non;     // specificed nonuniform temperature
    int starting_blk;              // starting block id for this boundary
    std::vector<double> no_blk;       // number of blocks for this boundary
    Vector3 r_omega;  // angular velocity of rotating surface
    Vector3 centre;   // a point on the axis of rotation
    Vector3 v_trans;  // translational velocity to superimpose
public:
    NonuniformTBC(Block *bdp, int which_boundary, vector<double>& T_non_, 
                              int starting_blk, vector<double>& no_blk_, 
                              Vector3 omega, Vector3 centre, Vector3 v_trans,
                              double emissivity);
    NonuniformTBC(const NonuniformTBC &bc);
    NonuniformTBC();
    NonuniformTBC & operator=(const NonuniformTBC &bc);
    virtual ~NonuniformTBC();
    virtual void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    virtual int apply_viscous(double t); // sets wall temperature
};
