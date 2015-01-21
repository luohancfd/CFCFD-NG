#ifndef BC_FIXED_T_HH
#define BC_FIXED_T_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_FIXED_T : public SolidBoundaryCondition
{

public:

    double Twall;

    BC_FIXED_T(double Twall);

    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
};

#endif // BC_FIXED_T_HH
