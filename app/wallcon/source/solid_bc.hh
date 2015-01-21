 #ifndef SOLID_BC_HH
#define SOLID_BC_HH

#include "solid_block.hh"
#include "solid_cell.hh"

//wanted boundary condtions
// FIXED_T
// ADIABATIC
// CONVECTIVE
// RADIATIVE
// CONV_RAD

class SolidBoundaryCondition {

public:
    virtual void print_type() {}
    virtual void apply_bc() {}
    virtual void apply_bc(SolidBlock *blk, int i, int j, int which_boundary) {}
    virtual void apply_bc(SolidBlock *blk, int which_boundary) {}

    //Constructors
    SolidBoundaryCondition() {}
};




#include "bc_adiabatic.hh"
#include "bc_fixed_t.hh"
#include "bc_fixed_q.hh"
#include "bc_convection.hh"
#include "bc_e3connection.hh"
#include "bc_userdef_q.hh"
#include "bc_userdef_t.hh"
#include "bc_radiation.hh"
#include "bc_conv_rad.hh"
#include "bc_timevarying_q.hh"
#include "bc_timevarying.hh"

#endif // SOLID_BC_HH
