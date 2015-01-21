#ifndef BC_FIXED_Q_HH
#define BC_FIXED_Q_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_FIXED_Q: public SolidBoundaryCondition
{
public:

    double qwall;

    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();

    BC_FIXED_Q(double qwall);
};

#endif // BC_FIXED_Q_HH
