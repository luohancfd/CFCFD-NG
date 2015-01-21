#ifndef BC_E3CONNECTION_HH
#define BC_E3CONNECTION_HH

#include "solid_bc.hh"
#include "solid_block.hh"

class BC_E3CONNECTION : public SolidBoundaryCondition
{
public:

    //Determine whether wall recievs flux(default 0) or temperature(1) for bc
    std::vector<double> connection_vector;
    
    BC_E3CONNECTION(); 
    BC_E3CONNECTION(std::vector<double> connection_vector);
    
    void apply_bc(SolidBlock *blk, int which_boundary);
    void print_type();
};

class BC_E3_FLUX : public BC_E3CONNECTION
{
public:    
    std::vector<double> flux_vector;

};
#endif // BC_E3CONNECTION_HH
