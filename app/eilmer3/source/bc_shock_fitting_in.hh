// bc_shock_fitting_in.hh

#include "bc.hh"

class ShockFittingInBC : public BoundaryCondition {
public:
    int inflow_condition_id; // index into the collection of inflow_conditions
public:
    ShockFittingInBC( Block &bdp, int which_boundary, int inflow_condition_id );
    ShockFittingInBC( const ShockFittingInBC &bc );
    virtual ~ShockFittingInBC();
    virtual int apply_inviscid( double t ); // copies from FlowCondition to ghost cells
    // default apply_viscous() (does nothing)
private:
    int calculate_shock_speed(FV_Cell *cL1, FV_Cell *cL0, FV_Cell *cR0, FV_Cell *cR1, FV_Cell *cR2, 
				           double lenL1, double lenL0, double lenR0, double lenR1, double lenR2, 
				           const FV_Interface *qL, FV_Interface *qR);
    int set_inflow_fluxes(FV_Interface *IFaceL, FV_Interface *IFaceR, double omegaz);
    inline double velocity_weighting_factor(double M);
};
