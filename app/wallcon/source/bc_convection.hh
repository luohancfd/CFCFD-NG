#ifndef BC_CONVECTION_HH
#define BC_CONVECTION_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_CONVECTION: public SolidBoundaryCondition
{
public:
    
    double h;
    double T_inf;
    
    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
    
    BC_CONVECTION(double h, double T_inf);
    
};

#endif // BC_CONVECTION_HH
