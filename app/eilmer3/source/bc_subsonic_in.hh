// bc_subsonic_in.hh

#include "bc.hh"

class SubsonicInBC : public BoundaryCondition {
public:
    int inflow_condition_id; // index into the collection of inflow_conditions
    int use_ideal_gas_relations;
public:
    SubsonicInBC( Block *bdp, int which_boundary, int inflow_condition_id,
		  int assume_ideal=0 );
    SubsonicInBC( const SubsonicInBC &bc );
    virtual ~SubsonicInBC();
    virtual int apply_inviscid( double t );
    // default apply_viscous() (does nothing)
private:
    int subsonic_inflow_properties( 
        const CFlowCondition *stagnation,
	CFlowCondition *inflow_state,
	double inflow_velocity );
};
