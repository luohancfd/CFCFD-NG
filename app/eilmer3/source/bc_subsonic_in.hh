// bc_subsonic_in.hh

#include "bc.hh"

class SubsonicInBC : public BoundaryCondition {
public:
    int inflow_condition_id; // index into the collection of inflow_conditions
    double mass_flux;
    double relax_factor;
    bool use_ideal_gas_relations;
public:
    SubsonicInBC(Block *bdp, int which_boundary, int inflow_condition_id, 
		 double mass_flux=0.0, double relax_factor=0.05, bool assume_ideal=false);
    SubsonicInBC(const SubsonicInBC &bc);
    SubsonicInBC();
    SubsonicInBC & operator=(const SubsonicInBC &bc);
    virtual ~SubsonicInBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    // default apply_viscous() (does nothing)
private:
    int subsonic_inflow_properties( 
        const CFlowCondition &stagnation,
	double dir_x, double dir_y, double dir_z,
	CFlowCondition &inflow_state,
	double inflow_pressure );

    // Our own copy of the stagnation condition that may get its pressure updated
    // if the mass-flux is in control.
    CFlowCondition gstagp; 
    double s0; // stagnation entropy
    double h0; // stagnation enthalpy
    void setup_stagnation_condition();
};
