// bc_shock_in.hh

#include "bc.hh"

class ShockInBC : public BoundaryCondition {
public:
    int inflow_condition_id; // index into the collection of inflow_conditions
public:
    ShockInBC( Block &bdp, int which_boundary, int inflow_condition_id );
    ShockInBC( const ShockInBC &bc );
    virtual ~ShockInBC();
    virtual int apply_inviscid( double t ); // copies from FlowCondition to ghost cells
    // default apply_viscous() (does nothing)
private:
    int shock_inflow_fluxes(FV_Interface *IFace, double omegaz);
};
