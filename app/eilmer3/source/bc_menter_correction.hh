// bc_menter_correction.hh

#include "cell.hh"
#include "block.hh"

double ideal_omega_at_wall(FV_Cell *cell);
double ideal_omega(FV_Cell *cell);
int apply_menter_boundary_correction(Block &bd, size_t ftl);
int apply_wilson_omega_correction(Block &bd, size_t ftl);

