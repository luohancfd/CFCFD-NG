// bc_jump_wall.cxx

#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_jump_wall.hh"
#include "bc_catalytic.hh"
#include "bc_menter_correction.hh"

//------------------------------------------------------------------------

JumpWallBC::JumpWallBC(Block *bdp, int which_boundary, double Twall, double sigma)
    : BoundaryCondition(bdp, which_boundary, FIXED_T),
      Twall(Twall), sigma(sigma)
{
    is_wall_flag = true;
}

JumpWallBC::JumpWallBC(const JumpWallBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      Twall(bc.Twall), sigma(bc.sigma)
{
    is_wall_flag = bc.is_wall_flag;
}

JumpWallBC::JumpWallBC()
    : BoundaryCondition(0, 0, FIXED_T),
      Twall(300.0), sigma(0.0)
{
    is_wall_flag = true;
}

JumpWallBC & JumpWallBC::operator=(const JumpWallBC &bc)
{
    BoundaryCondition::operator=(bc);
    Twall = bc.Twall; // Ok for self-assignment.
    sigma = bc.sigma;
    return *this;
}

JumpWallBC::~JumpWallBC() {}

void JumpWallBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "Twall= " << Twall << endl;
    cout << lead_in << "sigma= " << sigma << endl;
    return;
}

int JumpWallBC::apply_viscous(double t)
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    Gas_model & gm = *get_gas_model_ptr();
    int gas_status = 0;
    size_t nmodes = gm.get_number_of_modes();
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		// [TODO] here's where we need to evaluate the velocity and temperature
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
		// [TODO] maybe should we re-evaluate the thermo and transport coeffs?
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
		// [TODO] here's where we need to evaluate the velocity and temperature.
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
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
		// Evaluate the slip-velocity.
		double cell_half_width = cell->jLength / 2.0;
		fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		double v_cell_tangent = hypot(fs.vel.y, fs.vel.z);
		double dvdn = v_cell_tangent / cell_half_width;
		double c_bar = sqrt(8.0 * gm.R(*(fs.gas), gas_status) * fs.gas->T[0] / M_PI);
		double lambda_v = 2.0 * fs.gas->mu / (fs.gas->rho * c_bar);
		double v_slip = 2.0 * lambda_v * dvdn / sigma;
		// Now, recover the slip-velocity components.
		fs.vel.y *= v_slip/v_cell_tangent;
		fs.vel.z *= v_slip/v_cell_tangent;
		fs.vel.transform_to_local(IFace->n, IFace->t1, IFace->t2);
		// Evaluate effective (jump) temperature that the gas feels.
		double dTdn = (cell->fs->gas->T[0] - Twall) / cell_half_width;
		double lambda_T = 4.0 / (gm.gamma(*(fs.gas), gas_status) + 1.0) * 
		    fs.gas->k[0] / (fs.gas->rho * c_bar * gm.Cv(*(fs.gas), gas_status));
		double T_effective = Twall + 2.0 * lambda_T * dTdn / sigma;
		for ( size_t imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = T_effective;
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
		// [TODO] here's where we need to evaluate the velocity and temperature.
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
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
		// [TODO] here's where we need to evaluate the velocity and temperature.
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
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
		// [TODO] here's where we need to evaluate the velocity and temperature.
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
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
} // end JumpWallBC::apply_viscous()
