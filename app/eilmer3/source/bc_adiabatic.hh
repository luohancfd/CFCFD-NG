// bc_adiabatic.hh

#include "bc.hh"

class AdiabaticBC : public BoundaryCondition {
public:
    AdiabaticBC( Block *bdp, int which_boundary );
    AdiabaticBC( const AdiabaticBC &bc );
    AdiabaticBC();
    AdiabaticBC & operator=(const AdiabaticBC &bc);
    virtual ~AdiabaticBC();
    // default apply_inviscid() is just to reflect normal velocity
    virtual int apply_viscous( double t ); // sets wall T to interior T
};
