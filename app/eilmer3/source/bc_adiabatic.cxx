// bc_adiabatic.cxx

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_adiabatic.hh"
#include "bc_menter_correction.hh"

//------------------------------------------------------------------------

AdiabaticBC::AdiabaticBC( Block &bdp, int which_boundary )
    : BoundaryCondition(bdp, which_boundary, ADIABATIC, "AdiabaticBC",
			0, true, false, -1, -1, 0) 
{}

AdiabaticBC::AdiabaticBC( const AdiabaticBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.x_order, bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation) 
{}

AdiabaticBC::~AdiabaticBC() {}

int AdiabaticBC::apply_viscous( double t )
// Notes:
// We make the wall non-catalytic to ensure no heat transfer.
// This is not strictly correct to set the species here,
// rather qx and qy should be set to 0, however,
// it gives the identical effect on flow variables.
//
// Menter's slightly-rough-surface boundary condition as described
// in Wilcox 2006 text, eqn 7.36.
// We assume that the y2 in eqn 7.16 is the same as
// the height of our finite-volume cell.
{
    int i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;

    switch ( which_boundary ) {
    case NORTH:
	j = bdp.jmax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		cell = bdp.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
	    } // end i loop
	} // end for k
	break;
    case EAST:
	i = bdp.imax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		cell = bdp.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
	    } // end j loop
	} // end for k
	break;
    case SOUTH:
	j = bdp.jmin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		cell = bdp.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
	    } // end i loop
	} // end for k
	break;
    case WEST:
	i = bdp.imin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		cell = bdp.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
	    } // end j loop
	} // end for k
 	break;
    case TOP:
	k = bdp.kmax;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		cell = bdp.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
	    } // end j loop
	} // end for i
	break;
    case BOTTOM:
	k = bdp.kmin;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		cell = bdp.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
	    } // end j loop
	} // end for i
 	break;
    default:
	printf( "Error: apply_viscous not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
} // end AdiabaticBC::apply_viscous()

