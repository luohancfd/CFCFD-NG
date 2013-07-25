// bc_supersonic_in.hh

#include "bc.hh"

class SupersonicInBC : public BoundaryCondition {
public:
    int inflow_condition_id; // index into the collection of inflow_conditions
public:
    SupersonicInBC( Block *bdp, int which_boundary, int inflow_condition_id );
    SupersonicInBC( const SupersonicInBC &bc );
    SupersonicInBC();
    SupersonicInBC & operator=(const SupersonicInBC &bc);
    virtual ~SupersonicInBC();
    virtual int apply_convective( double t ); // copies from FlowCondition to ghost cells
    // default apply_viscous() (does nothing)
};
