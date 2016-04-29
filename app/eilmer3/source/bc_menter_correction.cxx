// bc_menter_correction.cxx
/// \brief Apply Menter boundary correction to the cells near solid walls.
///
/// Apply Menter's correction for the omega values at the wall, as described
/// in Menter's 1994 AIAA Journal paper, v.32, n.8, pp.1598-1605. 
///
/// PJ, October 2007
/// Wilson C, March 2015

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"


double ideal_omega_at_wall(FV_Cell *cell)
{
    Gas_data *wall_gas = cell->cell_at_nearest_wall->fs->gas;
    double d0 = cell->half_cell_width_at_wall;
    return 60.0 * (wall_gas->mu / wall_gas->rho) / (0.075 * d0 * d0);
}
