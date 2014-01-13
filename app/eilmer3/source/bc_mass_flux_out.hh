// bc_mass_flux_out.hh

#include "bc.hh"

class MassFluxOutBC : public BoundaryCondition {
public:
    double mass_flux;
    double external_pressure;
    double relax_factor;

public:
    MassFluxOutBC(Block *bdp, int which_boundary,
		  double mass_flux, double p_init, double relax_factor);
    MassFluxOutBC(const MassFluxOutBC &bc);
    MassFluxOutBC();
    MassFluxOutBC & operator=(const MassFluxOutBC &bc);
    virtual ~MassFluxOutBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    // default apply_viscous() (does nothing)
};
