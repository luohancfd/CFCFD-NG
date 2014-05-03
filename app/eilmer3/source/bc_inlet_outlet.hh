// bc_inletoutlet.hh
// This boundary condition is derived from Fixed Pressure Outlet
// boundary condition, The tke and omega value will be switched 
// between FixedValue amd ZeroGradient depending from the direction
// of mass flux of boundary faces
#include "bc.hh"

class InletOutletBC : public BoundaryCondition {
public:
    double Pout;
    double I_turb;     // turbulence intensity
    double u_turb_lam; // turbulence to laminar viscosity ratio
    double Tout;
    bool use_Tout;
    int x_order;
public:
    InletOutletBC(Block *bdp, int which_boundary, double Pout,
                double I_turb, double u_turb_lam,
		double Tout, bool use_Tout, int x_order);
    InletOutletBC(const InletOutletBC &bc);
    InletOutletBC();
    InletOutletBC & operator=(const InletOutletBC &bc);
    virtual ~InletOutletBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    virtual int apply_viscous(double t); 
};
