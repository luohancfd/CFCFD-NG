// bc_moving_wall.cxx

// Top and Bottom surfaces: rotate with z-axis k 
// West and East surfaces: rotate with x-axis i
// North and South Surfaces: rotate with y-axis j

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "math.h"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_moving_wall.hh"
#include "bc_catalytic.hh"
#include "bc_menter_correction.hh"

//------------------------------------------------------------------------

MovingWallBC::MovingWallBC(Block *bdp, int which_boundary, double _r_omega, double _emissivity)
    : BoundaryCondition(bdp, which_boundary, MOVING_WALL),
      r_omega(_r_omega)
{
    is_wall_flag = true;
    emissivity = _emissivity;
}

MovingWallBC::MovingWallBC(const MovingWallBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      r_omega(bc.r_omega)
{
    is_wall_flag = bc.is_wall_flag;
    emissivity = bc.emissivity;
}

MovingWallBC::MovingWallBC()
    : BoundaryCondition(bdp, which_boundary, MOVING_WALL),
      r_omega(0.0)
{
    is_wall_flag = true;
    emissivity = 1.0;
}

MovingWallBC & MovingWallBC::operator=(const MovingWallBC &bc)
{
    BoundaryCondition::operator=(bc); 
    r_omega = bc.r_omega;
    return *this;
}

MovingWallBC::~MovingWallBC() {}

void MovingWallBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "r_omega= " << r_omega << endl;
    return;
}

int MovingWallBC::apply_viscous(double t)
{
    double r_angle ; 
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    Block & bd = *bdp;
   // size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
    
    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                cell->copy_values_from(*cell, COPY_ALL_CELL_DATA, 0);
                r_angle = atan(cell->pos[0].x/cell->pos[0].z);
                fs.vel.z = -r_omega*sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].z*cell->pos[0].z)*sin(r_angle);
                fs.vel.x = r_omega*sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].z*cell->pos[0].z)*cos(r_angle);
                fs.vel.y = 0.0; 
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[NORTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                cell->copy_values_from(*cell, COPY_ALL_CELL_DATA, 0);
                r_angle = atan(cell->pos[0].z/cell->pos[0].y);
                fs.vel.y = -r_omega*sqrt(cell->pos[0].y*cell->pos[0].y+cell->pos[0].z*cell->pos[0].z)*sin(r_angle);
                fs.vel.z = r_omega*sqrt(cell->pos[0].y*cell->pos[0].y+cell->pos[0].z*cell->pos[0].z)*cos(r_angle);
                fs.vel.x = 0.0; 
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[EAST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                cell->copy_values_from(*cell, COPY_ALL_CELL_DATA, 0);
                r_angle = atan(cell->pos[0].x/cell->pos[0].z);
                fs.vel.z = -r_omega*sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].z*cell->pos[0].z)*sin(r_angle);
                fs.vel.x = r_omega*sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].z*cell->pos[0].z)*cos(r_angle);
                fs.vel.y = 0.0; 
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[SOUTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                cell->copy_values_from(*cell, COPY_ALL_CELL_DATA, 0);
                r_angle = atan(cell->pos[0].z/cell->pos[0].y);
                fs.vel.y = -r_omega*sqrt(cell->pos[0].y*cell->pos[0].y+cell->pos[0].z*cell->pos[0].z)*sin(r_angle);
                fs.vel.z = r_omega*sqrt(cell->pos[0].y*cell->pos[0].y+cell->pos[0].z*cell->pos[0].z)*cos(r_angle);
                fs.vel.x = 0.0; 
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[WEST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                cell->copy_values_from(*cell, COPY_ALL_CELL_DATA, 0);
                r_angle = atan(cell->pos[0].y/cell->pos[0].x);
                fs.vel.x = -r_omega*sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].y*cell->pos[0].y)*sin(r_angle);
                fs.vel.y = r_omega*sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].y*cell->pos[0].y)*cos(r_angle);
                fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[TOP]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
	for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
                cell->copy_values_from(*cell, COPY_ALL_CELL_DATA, 0);
                r_angle = atan(cell->pos[0].y/cell->pos[0].x);
                fs.vel.x = -r_omega * sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].y*cell->pos[0].y)*sin(r_angle);
                fs.vel.y = r_omega * sqrt(cell->pos[0].x*cell->pos[0].x+cell->pos[0].y*cell->pos[0].y)*cos(r_angle);
                fs.vel.z = 0.0;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bd.bcp[BOTTOM]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_viscous not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
} // end MovingWallTBC::apply_viscous()
