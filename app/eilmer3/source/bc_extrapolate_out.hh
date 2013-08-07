// bc_extrapolate_out.hh

#include "bc.hh"

class ExtrapolateOutBC : public BoundaryCondition {
public:
    ExtrapolateOutBC(Block *bdp, int which_boundary, int x_order, int sponge_flag=0);
    ExtrapolateOutBC(const ExtrapolateOutBC &bc);
    ExtrapolateOutBC();
    ExtrapolateOutBC & operator=(const ExtrapolateOutBC &bc);
    virtual ~ExtrapolateOutBC();
    virtual int apply_convective(double t); 
    // default apply_viscous() (does nothing)
};
