// bc_fixed_t.hh

#include "bc.hh"

class FixedTBC : public BoundaryCondition {
public:
    double Twall;
public:
    FixedTBC( Block *bdp, int which_boundary, double Twall );
    FixedTBC( const FixedTBC &bc );
    virtual ~FixedTBC();
    // default apply_inviscid() is just to reflect normal velocity
    virtual int apply_viscous( double t ); // sets wall temperature
};
