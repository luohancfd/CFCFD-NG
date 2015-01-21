#ifndef BC_TIMEVARYING_HH
#define BC_TIMEVARYING_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_TIMEVARYING: public SolidBoundaryCondition
{
public:
    
    std::vector<double> qwall;
    
    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
    
    BC_TIMEVARYING();
    BC_TIMEVARYING(std::vector<double> qwall);
};

#endif // BC_TIMEVARYING_HH
