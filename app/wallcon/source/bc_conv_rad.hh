#ifndef BC_CONV_RAD_HH
#define BC_CONV_RAD_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_CONV_RAD: public SolidBoundaryCondition
{
public:
    
    double eps; //Emissivity 
    double h; //Convective coefficient W.m^-2
    double T_inf; //Free stream temperature K
    
    double Tcell, distance, k;
    double sigma;
    double tol;
    
    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
    
    BC_CONV_RAD(double eps, double h,  double T_inf);

};

#endif // BC_CONV_RAD_HH
