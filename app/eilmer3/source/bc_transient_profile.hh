// bc_transient_profile.hh

#include "bc.hh"

class TransientProfileBC : public BoundaryCondition {
private:
    std::string filename;
    size_t n_profile;
    size_t nsp, nmodes;
    size_t ncell_for_profile;
    std::vector<CFlowCondition*> flow_profile;
public:
    TransientProfileBC(Block *bdp, int which_boundary, 
		    const std::string filename="profile.dat", size_t n_profile=1);
    TransientProfileBC(const TransientProfileBC &bc);
    TransientProfileBC();
    TransientProfileBC & operator=(const TransientProfileBC &bc);
    virtual ~TransientProfileBC();
    virtual void print_info(std::string lead_in);
    virtual int apply_convective(double t);
    // default apply_viscous() (does nothing)
};
