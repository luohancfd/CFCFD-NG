/** bc_ablating.cxx
 *
 * Created by Daniel Potter
 * Edited by Elise Fahy, beginning June 2013.
 * objective: species injection from the wall, by calculating
 * mass flow rates of species participating in surface reactions.
 *
 * For a flowfield with nsp species, there will be nsp+2 unknowns
 * to solve for: nsp species mass fractions, and density and
 * velocity at the wall. Hence, nsp+2 equations are required.
 *
 * For more details, see Park's 2001 paper.
 */

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>

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
#include "bc_ablating.hh"
#include "bc_menter_correction.hh"

//------------------------------------------------------------------------

// default constructor

AblatingBC::AblatingBC()
    : BoundaryCondition(0, 0, ABLATING)
      // Default-initialise everything else since we really can't use
      // this default-initialised BC.
{}

// normal constructor

AblatingBC::AblatingBC(Block *bdp, int which_boundary, double Twall, double _emissivity)
    : BoundaryCondition(bdp, which_boundary, ABLATING),
      Twall(Twall), mdot(mdot), max_iterations(1000000), tol(1.0e-6)
{
    is_wall_flag = true;
    sets_conv_flux_flag = true;
    ghost_cell_data_available = false;
    sets_visc_flux_flag = true;

    double T = Twall;          // how to get this from input script?
    TProfile.push_back(T);
    
    emissivity = _emissivity;

    // 0. Get gas-model pointer
    gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();

    // 2. Initialise the local gas-data structure (used for EOS calls)
    Q = new Gas_data(gmodel);

    // 3. Size the CFD cell mass-fraction vector
    cell_massf.resize(nsp);		// not sure if this is the right size...

    // 4. initialise the zero system components
    // total length of y-vectors = (nsp+1) because:
    // species from index 0 to (nsp-1), u_wall at index nsp, rho_wall at index (nsp+1)
    u0_index = nsp;
    rho_index = (nsp+1);
    y_guess.resize(nsp+1);
    y_out.resize(nsp+1);
    zero_solver = new NewtonRaphsonZF(nsp+1, tol, max_iterations, true);
    f_jac = 1.0;    // Jacobian scale factor
}

//copy constructor for AblatingBC class

AblatingBC::AblatingBC(const AblatingBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      Twall(bc.Twall), gmodel(bc.gmodel), u0_index(bc.u0_index), 
      max_iterations(bc.max_iterations), tol(bc.tol), f_jac(bc.f_jac)
{
    is_wall_flag = bc.is_wall_flag;
    emissivity = bc.emissivity;
    // 0. Get gas-model pointer
    // -> copied explicitly above
    int nsp = gmodel->get_number_of_species();
    // 1. Calculate the total mass flux from the given species-specific components
    // -> copied explicitly above
    // 2. Initialise the local gas-data structure (used for EOS calls)
    Q = new Gas_data(gmodel);
    // 3. Size the CFD cell mass-fraction vector
    cell_massf.resize( nsp );
    // 4. initialise the zero system components
    // -> u0_index already copied
    y_guess.resize(nsp+1);
    y_out.resize(nsp+1);
    zero_solver = new NewtonRaphsonZF(nsp+1, tol, max_iterations, true);
}


// Assignment operator

AblatingBC & AblatingBC::operator=(const AblatingBC &bc)
{
    BoundaryCondition::operator=(bc);
    if ( this != &bc ) {
	Twall = bc.Twall;
	emissivity = bc.emissivity;
	TProfile = bc.TProfile;
	ncell_for_profile = bc.ncell_for_profile;
	mdot_total = bc.mdot_total;
	// 0. Get gas-model pointer
	gmodel = bc.gmodel;
	size_t nsp = gmodel->get_number_of_species();
	// 1. Calculate the total mass flux from the given species-specific components
	// -> copied explicitly above
	// 2. Initialise the local gas-data structure (used for EOS calls)
	Q = new Gas_data(gmodel);
	// 3. Size the CFD cell mass-fraction vector
	cell_massf = std::vector<double>(bc.cell_massf);
	// 4. initialise the zero system components
	// -> u0_index already copied
	y_guess = std::vector<double>(bc.y_guess);
	y_out = std::vector<double>(bc.y_guess);
	zero_solver = new NewtonRaphsonZF(nsp+1, tol, max_iterations, true);
    }
    return *this;
}

// Default destructor

AblatingBC::~AblatingBC()
{}

void AblatingBC::print_info(std::string lead_in)
{
    BoundaryCondition::print_info(lead_in);
    cout << lead_in << "Twall= " << Twall << endl;
    cout << lead_in << "*** FIX-ME *** Elise. More configuration data should be written." << endl;
    return;
}

// apply_convective function definition

int AblatingBC::apply_convective(double t)
{
    // The default convective boundary condition is to reflect
    // the normal component of the velocity at the ghost-cell
    // centres -- slip-wall. 
    global_data &G = *get_global_data_ptr();
    size_t i, j, k;
    FV_Cell *src_cell, *dest_cell;
    FV_Interface *IFace;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[NORTH];
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i,j-1,k);
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[EAST];
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i-1,j,k);
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[SOUTH];
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i,j+1,k);
		dest_cell = bd.get_cell(i,j-2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i+1,j,k);
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[TOP];
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
		// ghost cell 2.
		src_cell = bd.get_cell(i,j,k-1);
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		// ghost cell 1.
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[BOTTOM];
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		// ghost cell 2.
		src_cell = bd.get_cell(i,j,k+1);
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		reflect_normal_velocity(dest_cell, IFace);
		if ( G.MHD ) {
		    reflect_normal_magnetic_field(dest_cell, IFace);
		}
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
} // end AblatingBC::apply_convective()

// calculate viscous fluxes

int AblatingBC::apply_viscous(double t)
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    size_t nmodes = get_gas_model_ptr()->get_number_of_modes();
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
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
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
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
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
} // end AblatingBC::apply_viscous()

// the functions required by the ZeroSystem are f and Jac.

int AblatingBC::f(const vector<double> &y_guess, vector<double> &G)
{
    /* Create the equation system for the ZeroSystem for a given y vector.
       - Species mass conservation: unknowns are massfs, rho_w, v_w
       - Total mass conservation: unknowns are rho_w, v_w
       - Momentum conservation: unknowns are rho_w, v_w                
       Need diffusion coefficients and species molar masses. */
    
    size_t iG=0;
    double mdot_c = 1.0;
    double R = gmodel->R(*Q);
 
    //for ( iG=0; iG<cell_massf.size(); ++iG ) {
    for ( iG=0; iG<2; ++iG ) {
	G[iG]=0.0;
    }

    // line 1: total mass conservation equation
    G[0] = y_guess[0]*y_guess[1] - mdot_c;

    // line 2: total momentum conservation equation
    G[1] = y_guess[0]*R*Twall + y_guess[0]*y_guess[1]*y_guess[1] - cell_momentum_flux;

    // line 3 -> (nsp+2): species mass conservation equations
    G[2] = 2.0*y_guess[0];

    return SUCCESS;
}

int AblatingBC::Jac(const vector<double> &y_guess, Valmatrix &dGdy)
{
    /* Create the Jacobian matrix for the ZeroSystem for a given y vector */
    double R = gmodel->R(*Q);

    //  Clear the Jacobian matrix
    for ( size_t i=0; i<nsp; ++i ) {
	for ( size_t j=0; j<nsp; ++j ) {
	    dGdy.set(i,j,0.0);
	}
    }
    
    // line 1: total mass conservation equation derivatives
    dGdy.set(0,0,y_guess[1]);
    dGdy.set(0,1,y_guess[0]);

    // line 2: total momentum conservation equation derivatives
    dGdy.set(1,0,(R*Twall + y_guess[1]*y_guess[1]));
    dGdy.set(1,1,(2*y_guess[0]*y_guess[1]));

    // line 3 -> (nsp+2): species mass conservation equations derivatives

    return SUCCESS;
}

// the other functions needed to solve the system:
// create y vector for G[y]=0

int AblatingBC::char_mass_flow()
{
    /* Calculate the char mass flow rate term for each species.
       This will require the reaction rates (alpha) - from a
       reaction file? - and the species mass, from the species
       files. Then these can feed straight into the species
       mass conservation equations. */

    return SUCCESS;
}

int AblatingBC::create_y_guess(FV_Cell *cell1, FV_Interface *wall, FV_Cell *cell0)
{
    /* line 1: rho_w
       line 2: v_w
       lines 3->(nsp+2): massf[isp] */

    size_t iy;

    // normal velocity 
    for ( iy=0; iy<2; ++iy ){
	    y_guess[iy] = 0;
	}
	      y_guess[0] = Q->rho;
	  y_guess[1] = cell0->fs->vel.x;
	  //species mass densities
	  for ( iy=0; iy<Q->massf.size(); ++iy ){
	      y_guess[iy+2] = Q->massf[iy];
	  }

	      cell_un = cell1->fs->vel.x;
	  cell_local_vel = cell1->fs->vel;  // FIX-ME copy components instead ??
	  cell1->fs->vel.transform_to_global(wall->n, wall->t1, wall->t2);
	  cell_mass_flux = cell1->fs->gas->rho * cell_un;
	  cell_momentum_flux = cell1->fs->gas->p + cell1->fs->gas->rho * cell_un * cell_un;

	  return SUCCESS;
	  }

    // solve the system with a zero-solving method
    int AblatingBC::solve_system(FV_Cell *cell1, FV_Interface *wall, FV_Cell *cell)
    {
	// how do I pull in the y_guess from the other function?

	// solve system
	if ( zero_solver->solve( *this, y_guess, y_out ) ) {
	    cout << "AblatingBC::solve_system()" << endl
		 << "Zero solver has failed, bailing out!" << endl;
	    exit( FAILURE );
	}

	// map results
	// mass fractions, density, velocity, evaluate new thermo state

	return SUCCESS;
    }
