// bc_fixed_p_out.hh

#include "bc.hh"

class FixedPOutBC : public BoundaryCondition {
public:
    double Pout;
    double Tout;
    bool use_Tout;
    int x_order;

public:
    FixedPOutBC(Block *bdp, int which_boundary, double Pout,
		double Tout, bool use_Tout, int x_order);
    FixedPOutBC(const FixedPOutBC &bc);
    FixedPOutBC();
    FixedPOutBC & operator=(const FixedPOutBC &bc);
    virtual ~FixedPOutBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    // default apply_viscous() (does nothing)
};
