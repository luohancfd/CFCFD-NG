// bc_surface_energy_balance.cxx

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <valarray>

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas_data.hh"
#include "../../../lib/gas/models/gas-model.hh"
#include "../../../lib/gas/models/physical_constants.hh"
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "bc_catalytic.hh"
#include "bc_surface_energy_balance.hh"
#include "bc_menter_correction.hh"

//------------------------------------------------------------------------

SurfaceEnergyBalanceBC::SurfaceEnergyBalanceBC(Block *bdp, int which_boundary,
                                                           double emissivity)
    : BoundaryCondition(bdp, which_boundary, SEB, emissivity),
      tol(1.0e-4), max_iterations(100), f_relax(0.05)
{
    is_wall_flag = true;
    // Ensure that emissivity is between 0 and 1 as required
    if ( emissivity < 0.0 || emissivity > 1.0 ) {
    	cerr << "SurfaceEnergyBalanceBC::SurfaceEnergyBalanceBC()" << endl
    	     << "Given value of emissivity (" << emissivity << ") "
    	     << "is not between 0 and 1 as required." << endl
    	     << "Exiting program." << endl;
    	exit( BAD_INPUT_ERROR );
    }
    
    // Initialise the local gas-data structure (used for EOS calls)
    Gas_model * gmodel = get_gas_model_ptr();
    Q = new Gas_data(gmodel);
}

SurfaceEnergyBalanceBC::SurfaceEnergyBalanceBC(const SurfaceEnergyBalanceBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      tol(bc.tol), max_iterations(bc.max_iterations),
      f_relax(bc.f_relax)
{
    is_wall_flag = bc.is_wall_flag;
    emissivity = bc.emissivity;
    // Initialise the local gas-data structure (used for EOS calls)
    Gas_model * gmodel = get_gas_model_ptr();
    Q = new Gas_data(gmodel);
}

SurfaceEnergyBalanceBC & 
SurfaceEnergyBalanceBC::operator=(const SurfaceEnergyBalanceBC &bc)
{
    if ( this != &bc ) {
	is_wall_flag = bc.is_wall_flag;
	emissivity = bc.emissivity;
	tol = bc.tol;
	max_iterations = bc.max_iterations;
	f_relax = bc.f_relax;
	Q = new Gas_data(*(bc.Q));
    }
    return *this;
}

SurfaceEnergyBalanceBC::SurfaceEnergyBalanceBC()
    : BoundaryCondition(0, 0, SEB),
      tol(1.0e-4), max_iterations(100), f_relax(0.05)
{
    is_wall_flag = true;
    emissivity = 1.0;
    // Initialise the local gas-data structure (used for EOS calls)
    Gas_model * gmodel = get_gas_model_ptr();
    Q = new Gas_data(gmodel);
}

SurfaceEnergyBalanceBC::~SurfaceEnergyBalanceBC() 
{
    delete Q;
}

void SurfaceEnergyBalanceBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "emissivity= " << emissivity << endl;
    cout << lead_in << "tol= " << tol << endl;
    cout << lead_in << "max_iterations= " << max_iterations << endl;
    cout << lead_in << "f_relax= " << f_relax << endl;
    return;
}

int SurfaceEnergyBalanceBC::apply_viscous(double t)
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    size_t index;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
	for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		index = (bd.jmax-jmin+1)*(imax-imin+1)*(k-kmin) + \
			(imax-imin+1)*(j-jmin) + (i-imin);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[NORTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		if ( solve_for_wall_temperature( IFace, cell, index ) ) {
		    cerr << "SurfaceEnergyBalanceBC::apply_viscous()" << endl
		    	 << "solve_for_wall_temperature() failed for index:"
		    	 << index << " of block: " << bd.id << ", boundary: " 
		    	 << which_boundary << endl;
		    return FAILURE;
		}
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
	    	index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + 
			(imax-imin+1)*(j-jmin) + (i-imin);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		if ( solve_for_wall_temperature(IFace, cell, index) ) {
		    cerr << "SurfaceEnergyBalanceBC::apply_viscous()" << endl
		    	 << "solve_for_wall_temperature() failed for index:"
		    	 << index << " of block: " << bd.id << ", boundary: " 
		    	 << which_boundary << endl;
		    return FAILURE;
		}
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
	    	index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + 
			(imax-imin+1)*(j-jmin) + (i-imin);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		if ( solve_for_wall_temperature(IFace, cell, index) ) {
		    cerr << "SurfaceEnergyBalanceBC::apply_viscous()" << endl
		    	 << "solve_for_wall_temperature() failed for index:"
		    	 << index << " of block: " << bd.id << ", boundary: " 
		    	 << which_boundary << endl;
		    return FAILURE;
		}
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
	    	index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + 
			(imax-imin+1)*(j-jmin) + (i-imin);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		if ( solve_for_wall_temperature(IFace, cell, index) ) {
		    cerr << "SurfaceEnergyBalanceBC::apply_viscous()" << endl
		    	 << "solve_for_wall_temperature() failed for index:"
		    	 << index << " of block: " << bd.id << ", boundary: " 
		    	 << which_boundary << endl;
		    return FAILURE;
		}
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
	    	index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + 
			(imax-imin+1)*(j-jmin) + (i-imin);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[TOP];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		if ( solve_for_wall_temperature(IFace, cell, index) ) {
		    cerr << "SurfaceEnergyBalanceBC::apply_viscous()" << endl
		    	 << "solve_for_wall_temperature() failed for index:"
		    	 << index << " of block: " << bd.id << ", boundary: " 
		    	 << which_boundary << endl;
		    return FAILURE;
		}
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
	    	index = (jmax-jmin+1)*(imax-imin+1)*(k-kmin) + 
			(imax-imin+1)*(j-jmin) + (i-imin);
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[BOTTOM];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		if ( solve_for_wall_temperature(IFace, cell, index) ) {
		    cerr << "SurfaceEnergyBalanceBC::apply_viscous()" << endl
		    	 << "solve_for_wall_temperature() failed for index:"
		    	 << index << " of block: " << bd.id << ", boundary: " 
		    	 << which_boundary << endl;
		    return FAILURE;
		}
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
} // end SurfaceEnergyBalanceBC::apply_viscous()

int
SurfaceEnergyBalanceBC::
solve_for_wall_temperature( FV_Interface * IFace, FV_Cell * cell_one, size_t index )
{
    size_t iteration;
    double q_total, q_total_prev = 0.0;
    double Twall = IFace->fs->gas->T[0];
    double Twall_prev = IFace->fs->gas->T[0];
    
    // 0. Check that given interface temperature is a reasonable initial guess
    if ( IFace->fs->gas->T[0] >= cell_one->fs->gas->T[0] ) {
//    	cout << "SurfaceEnergyBalanceBC::solve_for_wall_temperature()" << endl
//    	     << "Given initial wall temperature gradient is zero or negative." << endl
//    	     << "Trying Twall_prev = Tcell / 2" << endl;
    	Twall_prev = cell_one->fs->gas->T[0] / 2.0;
    	for ( size_t itm=0; itm<IFace->fs->gas->T.size(); ++itm )
    	    IFace->fs->gas->T[itm] = Twall_prev;
    }
    
    // NOTE: convergence test is on heat-flux not temperature as q is very 
    //       sensitive to small Twall perturbations
    for ( iteration=0; iteration<max_iterations; ++iteration ) {
    	// FIXME: need to re-evalate the thermo, diffusion and transport properties on the IFace first
    	// 0. Update IFace
    	this->update_interface_properties( IFace );
    	// 1. Fill in q vector for this interface
    	if ( compute_cell_interface_surface_heat_flux( IFace, cell_one, index ) ) {
    	    cerr << "SurfaceEnergyBalanceBC::solve_for_wall_temperature()" << endl
    	         << "compute_cell_interface_surface_heat_flux() failed for index: "
    	         << index << endl;
    	    return FAILURE;
    	}
    	// 2. calculate total heat flux incident on this interface
    	q_total = q_rad[index] + q_cond[index] + q_diff[index];
    	// 3. calculate Twall by assuming radiative equilibrium at the wall
    	Twall = pow( q_total / ( emissivity * PC_sigma_SI ), 0.25 );
        // 4. calculate new guess as a weighted average so we don't take too big a step
        Twall = f_relax * Twall + ( 1.0 - f_relax ) * Twall_prev;
// cout << "iteration = " << iteration << ", q_total = " << q_total << ", Twall = " << Twall << ", error = " << fabs(Twall-Twall_prev)/Twall_prev << endl;
    	// 5. set new Twall on the interface
    	for ( size_t itm=0; itm<IFace->fs->gas->T.size(); ++itm )
    	    IFace->fs->gas->T[itm] = Twall;
    	// 6. check for convergence
    	if ( fabs(q_total-q_total_prev)/q_total < tol ) break;
    	// 7. set Twall_prev and q_total_prev
    	Twall_prev = Twall;
    	q_total_prev = q_total;
    }
    
    if ( iteration==max_iterations ) {
    	cerr << "SurfaceEnergyBalanceBC::solve_for_wall_temperature()" << endl
    	     << "Maximum iteration limit reached before convergence." << endl
    	     << "Twall = " << Twall << ", max_iterations = " << max_iterations
    	     << ", tol = " << tol << endl;
    	return FAILURE;
    }
    
    return SUCCESS;
}

void
SurfaceEnergyBalanceBC::
update_interface_properties(FV_Interface * IFace)
{
    // Just the temperatures and mass-fractions (if catalytic) should have changed
    Gas_model * gm = get_gas_model_ptr();
    
    // 1. Thermo properties
    gm->eval_thermo_state_rhoT(*(IFace->fs->gas));
    
    // 2. Transport coefficients
    gm->eval_transport_coefficients(*(IFace->fs->gas));
    
    // 3. Diffusion coefficients
    gm->eval_diffusion_coefficients(*(IFace->fs->gas));
    
    return;
}
