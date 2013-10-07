// bc_adiabatic.hh

#include "bc.hh"

class AdiabaticBC : public BoundaryCondition {
public:
    AdiabaticBC(Block *bdp, int which_boundary, double emissivity);
    AdiabaticBC(const AdiabaticBC &bc);
    AdiabaticBC();
    AdiabaticBC & operator=(const AdiabaticBC &bc);
    virtual ~AdiabaticBC();
    // default apply_convective() is just to reflect normal velocity
    virtual int apply_viscous(double t); // sets wall T to interior T
};
