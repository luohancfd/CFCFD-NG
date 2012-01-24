// bc_static_profile.hh

#include "bc.hh"

class StaticProfileBC : public BoundaryCondition {
private:
    std::string filename;
    int nsp, nmodes;
    int ncell_for_profile;
    std::vector<CFlowCondition*> flow_profile;
public:
    StaticProfileBC( Block &bdp, int which_boundary, 
		     const std::string filename="profile.dat" );
    StaticProfileBC( const StaticProfileBC &bc );
    virtual ~StaticProfileBC();
    virtual int apply_inviscid( double t );
    // default apply_viscous() (does nothing)
};
