// bc_transient_profile.hh

#include "bc.hh"

class TransientProfileBC : public BoundaryCondition {
private:
    std::string filename;
    size_t n_profile;
    size_t nsp, nmodes;
    size_t ncell_for_profile;
    std::ifstream fstrm;
    std::string text;
    double t0, t1;
    std::vector<CFlowCondition*> flow_profile_t0;
    std::vector<CFlowCondition*> flow_profile_t1;
    CFlowCondition *cf_interp;
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
private:
    void read_profile_into_t0();
    void read_profile_into_t1();
    void update_profile_t0_with_t1();
    void interpolate_condition(CFlowCondition &cf0, CFlowCondition &cf1, double w, CFlowCondition &cf_intep);
};
