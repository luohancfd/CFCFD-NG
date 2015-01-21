#ifndef BC_ADIABATIC_HH
#define BC_ADIABATIC_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_ADIABATIC : public SolidBoundaryCondition
{
public:

    BC_ADIABATIC();
    
    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
};

#endif // BC_ADIABATIC_HH
