// bc_standdard_wall_function.hh

#include "cell.hh"
#include "block.hh"

int wall_function_correction(Block &bd, size_t ftl);
int apply_turbulent_model_for_wall_function(Block &bd);

// Helper functions
void correction_adiabatic_wall(FV_Cell *cell, FV_Interface *IFace, size_t bc_type);
void correction_fixedt_wall(FV_Cell *cell, FV_Interface *IFace, size_t bc_type);
