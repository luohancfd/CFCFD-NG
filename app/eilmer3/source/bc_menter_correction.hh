// bc_menter_correction.hh

#include "cell.hh"
#include "block.hh"

double ideal_omega_at_wall(FV_Cell *cell);
double ideal_omega(FV_Cell *cell);
int apply_menter_boundary_correction(Block &bdp);

# define WILSON_OMEGA_FILTER 0
# if WILSON_OMEGA_FILTER == 1
int apply_wilson_omega_correction(Block &bdp);
# endif
