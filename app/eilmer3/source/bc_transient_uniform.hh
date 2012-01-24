// bc_transient_uniform.hh

#include "bc.hh"

class TransientUniformBC : public BoundaryCondition {
private:
    std::string filename;
    std::vector<double> tta, pa, ua;
    std::vector<double> va, wa, Ta;
    std::vector<double> tkea, omegaa;
    std::vector<std::vector<double> > massfa;
    Gas_model *gmodel;
    int nsample, nsp, nmodes;
public:
    TransientUniformBC( Block &bdp, int which_boundary, 
			std::string filename="transient_uniform.dat" );
    TransientUniformBC( const TransientUniformBC &bc );
    virtual ~TransientUniformBC();
    virtual int apply_inviscid( double t );
    // default apply_viscous() (does nothing)
};
