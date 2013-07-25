// bc_fixed_p_out.hh

#include "bc.hh"

class FixedPOutBC : public BoundaryCondition {
public:
    double Pout;
public:
    FixedPOutBC( Block *bdp, int which_boundary, double Pout, int x_order );
    FixedPOutBC( const FixedPOutBC &bc );
    FixedPOutBC();
    FixedPOutBC & operator=(const FixedPOutBC &bc);
    virtual ~FixedPOutBC();
    virtual int apply_convective( double t );
    // default apply_viscous() (does nothing)
};
