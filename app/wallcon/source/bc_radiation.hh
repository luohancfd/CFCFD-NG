#ifndef BC_RADIATION_HH
#define BC_RADIATION_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_RADIATION: public SolidBoundaryCondition
{
public:

    double eps;
    double T_inf;

    double Tcell, distance, k;
    double sigma;
    double tol;

    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();

    BC_RADIATION(double eps, double T_inf);

};

#endif // BC_RADIATION_HH
