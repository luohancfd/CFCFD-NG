#ifndef BC_TIMEVARYING_Q_HH
#define BC_TIMEVARYING_Q_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_TIMEVARYING_Q: public SolidBoundaryCondition
{
public:
    
    std::vector<double> qwall;
    
    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
    
    BC_TIMEVARYING_Q();
    BC_TIMEVARYING_Q(std::vector<double> qwall);
};

#endif // BC_TIMEVARYING_Q_HH
