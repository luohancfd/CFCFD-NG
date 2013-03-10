// bc_surface_energy_balance.hh

#include "bc.hh"

class SurfaceEnergyBalanceBC : public BoundaryCondition {
public:
    SurfaceEnergyBalanceBC( Block *bdp, int which_boundary, double epsilon );
    SurfaceEnergyBalanceBC( const SurfaceEnergyBalanceBC &bc );
    SurfaceEnergyBalanceBC();
    SurfaceEnergyBalanceBC & operator=(const SurfaceEnergyBalanceBC &bc);
    virtual ~SurfaceEnergyBalanceBC();
    // default apply_inviscid() is just to reflect normal velocity
    virtual int apply_viscous( double t ); // sets wall temperature
private:
    int solve_for_wall_temperature( FV_Interface * IFace,
    				    FV_Cell * cell_one,
    				    int index );
    void update_interface_properties( FV_Interface * IFace );
private:
    Gas_data *Q;
    double epsilon;
    double tol;
    int max_iterations;
    double f_relax;
};
