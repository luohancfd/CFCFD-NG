#ifndef BC_USERDEF_T_HH
#define BC_USERDEF_T_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_USERDEF_T: public SolidBoundaryCondition
{
public:
    
    std::vector<double> Twall;

    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
    
    BC_USERDEF_T();
    BC_USERDEF_T(std::vector<double> Twall);
};

#endif // BC_USERDEF_T_HH
