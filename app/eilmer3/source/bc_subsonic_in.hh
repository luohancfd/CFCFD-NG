// bc_subsonic_in.hh

#include "bc.hh"

class SubsonicInBC : public BoundaryCondition {
public:
    int inflow_condition_id; // index into the collection of inflow_conditions
    double mass_flux;
    double relax_factor;
    double p0_min;
    double p0_max;
    subsonic_in_direction_t direction_type; 
    std::vector<double> direction_vector;
    double direction_alpha;
    double direction_beta;
    bool use_ideal_gas_relations;
public:
    SubsonicInBC(Block *bdp, int which_boundary, int inflow_condition_id, 
		 double mass_flux=0.0, double relax_factor=0.05,
		 double p0_min=0.25e5, double p0_max=2.0e5,
		 subsonic_in_direction_t direction_type_=SUBSONIC_IN_NORMAL,
		 std::vector<double> direction_vector_=std::vector<double>(3, 0.0),
		 double direction_alpha_=0.0, double direction_beta_=0.0,
		 bool assume_ideal=false);
    SubsonicInBC(const SubsonicInBC &bc);
    SubsonicInBC();
    SubsonicInBC & operator=(const SubsonicInBC &bc);
    virtual ~SubsonicInBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    // default apply_viscous() (does nothing)
private:
    double subsonic_inflow_properties(const CFlowCondition &stagnation,
				      CFlowCondition &inflow_state,
				      double inflow_pressure );

    // Our own copy of the stagnation condition that may get its pressure updated
    // if the mass-flux is in control.
    CFlowCondition gstagp; 
    double s0; // stagnation entropy
    double h0; // stagnation enthalpy
    void setup_stagnation_condition();
};
