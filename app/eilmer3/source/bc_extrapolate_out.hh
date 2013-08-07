// bc_extrapolate_out.hh

#include "bc.hh"

class ExtrapolateOutBC : public BoundaryCondition {
public:
    int x_order;
    int sponge_flag;

public:
    ExtrapolateOutBC(Block *bdp, int which_boundary, int x_order, int sponge_flag=0);
    ExtrapolateOutBC(const ExtrapolateOutBC &bc);
    ExtrapolateOutBC();
    ExtrapolateOutBC & operator=(const ExtrapolateOutBC &bc);
    virtual ~ExtrapolateOutBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t); 
    // default apply_viscous() (does nothing)
};
