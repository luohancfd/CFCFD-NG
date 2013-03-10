// bc_fstc.hh

#include "bc.hh"

class fstcBC : public BoundaryCondition {
private:
    std::string filename;
    std::vector<double> fstc_TProfile;
    int ncell_for_profile;
public:
    fstcBC( Block *bdp, int which_boundary, const std::string filename="fstc_temp.txt" );
    fstcBC( const fstcBC &bc );
    virtual ~fstcBC();
    // default apply_inviscid() is just to reflect normal velocity
    virtual int apply_viscous( double t ); // sets wall temperature
};
