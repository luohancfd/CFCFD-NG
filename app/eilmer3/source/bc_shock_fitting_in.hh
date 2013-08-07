// bc_shock_fitting_in.hh

#include "bc.hh"

class ShockFittingInBC : public BoundaryCondition {
public:
    int inflow_condition_id; // index into the collection of inflow_conditions
public:
    ShockFittingInBC(Block *bdp, int which_boundary, int inflow_condition_id);
    ShockFittingInBC(const ShockFittingInBC &bc);
    ShockFittingInBC();
    ShockFittingInBC & operator=(const ShockFittingInBC &bc);
    virtual ~ShockFittingInBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    virtual int apply_viscous(double t);
private:
    int calculate_shock_speed(const FV_Cell &cL0, const FV_Cell &cR0,
			      const FV_Cell &cR1, const FV_Cell &cR2, 
			      double lenL0, double lenR0, double lenR1, double lenR2, 
			      FV_Interface &IFaceR);
 };
