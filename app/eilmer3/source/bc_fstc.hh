// bc_fstc.hh

#include "bc.hh"

class fstcBC : public BoundaryCondition {
private:
    std::string filename;
    std::vector<double> fstc_TProfile;
    unsigned int ncell_for_profile;
public:
    fstcBC(Block *bdp, int which_boundary, const std::string filename="fstc_temp.txt");
    fstcBC(const fstcBC &bc);
    fstcBC();
    fstcBC & operator=(const fstcBC &bc);
    virtual ~fstcBC();
    virtual void print_info(std::string lead_in);
    // default apply_convective() is just to reflect normal velocity
    virtual int apply_viscous(double t); // sets wall temperature
};
