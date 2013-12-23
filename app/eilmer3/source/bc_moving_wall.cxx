// bc_moving_wall.cxx
// Rotating or translating surface for Kan Qin's gas bearing study.
// KQ and PJ, November 2013.

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "../../../lib/geometry2/source/geom.hh"
#include "math.h"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_moving_wall.hh"
#include "bc_catalytic.hh"
#include "bc_menter_correction.hh"

//------------------------------------------------------------------------

MovingWallBC::MovingWallBC(Block *bdp, int which_boundary, Vector3 _r_omega,
			   Vector3 _centre, Vector3 _v_trans, double Twall, double _emissivity)
    : BoundaryCondition(bdp, which_boundary, MOVING_WALL),
      r_omega(_r_omega), centre(_centre), v_trans(_v_trans), Twall(Twall)
{
    is_wall_flag = true;
    emissivity = _emissivity;
}

MovingWallBC::MovingWallBC(const MovingWallBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      r_omega(bc.r_omega), centre(bc.centre), v_trans(bc.v_trans), Twall(bc.Twall)
{
    is_wall_flag = bc.is_wall_flag;
    emissivity = bc.emissivity;
}

MovingWallBC::MovingWallBC()
    : BoundaryCondition(bdp, which_boundary, MOVING_WALL),
      r_omega(Vector3(0.0,0.0,0.0)), centre(Vector3(0.0,0.0,0.0)), v_trans(Vector3(0.0,0.0,0.0)), Twall(300.0)
{
    is_wall_flag = true;
    emissivity = 1.0;
}

MovingWallBC & MovingWallBC::operator=(const MovingWallBC &bc)
{
    BoundaryCondition::operator=(bc); 
    r_omega.x = bc.r_omega.x; r_omega.y = bc.r_omega.y; r_omega.z = bc.r_omega.z;
    centre.x = bc.centre.x; centre.y = bc.centre.y; centre.z = bc.centre.z;
    v_trans.x = bc.v_trans.x; v_trans.y = bc.v_trans.y; v_trans.z = bc.v_trans.z;
    Twall = bc.Twall;
    return *this;
}

MovingWallBC::~MovingWallBC() {}

void MovingWallBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "r_omega= " << r_omega << endl;
    cout << lead_in << "centre= " << centre << endl;
    cout << lead_in << "v_trans= " << v_trans << endl;
    cout << lead_in << "Twall= " << Twall << endl;
    return;
}

int MovingWallBC::apply_viscous(double t)
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    Block & bd = *bdp;
    // Implement the boundary condition as a spinning disc with superimposed
    // translation velocity.  It's a bit odd but we intend to use either the 
    // rotational velocity OR the translational velocity (not both) in any
    // particular application.
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		// The following vector expression will involve the creation of temporaries
		// the might consume significant CPU time.  If that becomes a problem, we'll
		// expand the computation into component form.
		fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
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
		fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
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
		fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
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
		fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
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
		fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
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
		fs.vel = cross(r_omega, IFace->pos - centre) + v_trans;
                for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
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
} // end MovingWallBC::apply_viscous()
