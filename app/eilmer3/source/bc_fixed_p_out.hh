// bc_fixed_p_out.hh

#include "bc.hh"

class FixedPOutBC : public BoundaryCondition {
public:
    double Pout;
public:
    FixedPOutBC( Block &bdp, int which_boundary, double Pout );
    FixedPOutBC( const FixedPOutBC &bc );
    virtual ~FixedPOutBC();
    virtual int apply_inviscid( double t );
    // default apply_viscous() (does nothing)
};
