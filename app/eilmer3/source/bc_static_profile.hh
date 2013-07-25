// bc_static_profile.hh

#include "bc.hh"

class StaticProfileBC : public BoundaryCondition {
private:
    std::string filename;
    size_t n_profile;
    size_t nsp, nmodes;
    size_t ncell_for_profile;
    std::vector<CFlowCondition*> flow_profile;
public:
    StaticProfileBC( Block *bdp, int which_boundary, 
		     const std::string filename="profile.dat", size_t n_profile=1 );
    StaticProfileBC( const StaticProfileBC &bc );
    StaticProfileBC();
    StaticProfileBC & operator=(const StaticProfileBC &bc);
    virtual ~StaticProfileBC();
    virtual int apply_convective( double t );
    // default apply_viscous() (does nothing)
};
