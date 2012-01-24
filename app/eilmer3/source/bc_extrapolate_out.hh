// bc_extrapolate_out.hh

#include "bc.hh"

class ExtrapolateOutBC : public BoundaryCondition {
public:
    ExtrapolateOutBC( Block &bdp, int which_boundary, int sponge_flag=0 );
    ExtrapolateOutBC( const ExtrapolateOutBC &bc );
    virtual ~ExtrapolateOutBC();
    virtual int apply_inviscid( double t ); // copies interior flow to ghost cells
    // default apply_viscous() (does nothing)
};
