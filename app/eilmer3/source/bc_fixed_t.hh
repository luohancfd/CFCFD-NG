// bc_fixed_t.hh

#include "bc.hh"

class FixedTBC : public BoundaryCondition {
public:
    double Twall;
public:
    FixedTBC(Block *bdp, int which_boundary, double Twall);
    FixedTBC(const FixedTBC &bc);
    FixedTBC();
    FixedTBC & operator=(const FixedTBC &bc);
    virtual ~FixedTBC();
    virtual void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    virtual int apply_viscous(double t); // sets wall temperature
};
