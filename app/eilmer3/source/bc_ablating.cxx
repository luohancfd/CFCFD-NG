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

AblatingBC::AblatingBC(Block *bdp, int which_boundary, double Twall, 
		       vector<double> &mdot)
    : BoundaryCondition(bdp, which_boundary, ABLATING),
      Twall(Twall), mdot(mdot), max_iterations(1000000), tol(1.0e-6)
{
    is_wall_flag = true;

    double T = Twall;
    TProfile.push_back(T);
    
    // 0. Get gas-model pointer
    gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    // 1. Calculate the total mass flux from the given species-specific components
    mdot_total = 0.0;
    
    for ( size_t isp=0; isp<mdot.size(); ++isp )
    	mdot_total += mdot[isp];
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
    // f_jac is the Jacobian scale factor
    f_jac = 1.0;
}

//copy constructor for AblatingBC class

AblatingBC::AblatingBC(const AblatingBC &bc)
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code),
      Twall(bc.Twall), mdot(bc.mdot), mdot_total(bc.mdot_total),
      gmodel(bc.gmodel), u0_index(bc.u0_index), 
      max_iterations(bc.max_iterations), tol(bc.tol), f_jac(bc.f_jac)
{
    is_wall_flag = bc.is_wall_flag;
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
	mdot = bc.mdot;
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
	y_guess = std::valarray<double>(bc.y_guess);
	y_out = std::valarray<double>(bc.y_guess);
	zero_solver = new NewtonRaphsonZF(nsp+1, tol, max_iterations, true);
    }
    return *this;
}

// Default destructor

AblatingBC::~AblatingBC()
{
    delete zero_solver;
    delete Q;
}

// apply_convective function definition

int AblatingBC::apply_convective(double t)
{
    // Calculate the ghost cell flow-states from the given wall temperature
    // and species-specific mass flux.
    size_t i, j, k;
    FV_Cell *src_cell, *dest_cell;
    FV_Interface *IFace;
    Block & bd = *bdp;

    switch ( which_boundary ) {
    case NORTH:
	j = bd.jmax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[NORTH];
		// ghost cell 1
		dest_cell = bd.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bd.get_cell(i,j+1,k);
		dest_cell = bd.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bd.imax;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[EAST];
		// ghost cell 1
		dest_cell = bd.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bd.get_cell(i+1,j,k);
		dest_cell = bd.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bd.jmin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (i = bd.imin; i <= bd.imax; ++i) {
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[SOUTH];
		// ghost cell 1
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bd.get_cell(i,j-1,k);
		dest_cell = bd.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bd.imin;
        for (k = bd.kmin; k <= bd.kmax; ++k) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		// ghost cell 1
		dest_cell = bd.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bd.get_cell(i-1,j,k);
		dest_cell = bd.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bd.kmax;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		// ghost cell 1
		dest_cell = bd.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bd.get_cell(i,j,k+1);
		dest_cell = bd.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bd.kmin;
        for (i = bd.imin; i <= bd.imax; ++i) {
	    for (j = bd.jmin; j <= bd.jmax; ++j) {
		src_cell = bd.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		// ghost cell 1
		dest_cell = bd.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bd.get_cell(i,j,k-1);
		dest_cell = bd.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE, 0);
	    } // end j loop
	} // for i
 	break;
    default:
	printf( "Error: apply_inviscid not implemented for boundary %d\n", 
		which_boundary );
	return NOT_IMPLEMENTED_ERROR;
    }
    return SUCCESS;
}

// apply_viscous function definition

int AblatingBC::apply_viscous(double t)
{
    size_t i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    size_t nmodes = gmodel->get_number_of_modes();
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
		for ( size_t imode=0; imode < nmodes; ++imode ) { 
		    fs.gas->T[imode] = TProfile[i-bd.imin];
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
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[EAST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
		    fs.gas->T[imode] = TProfile[j-bd.jmin];
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
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[SOUTH];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
		    fs.gas->T[imode] = TProfile[i-bd.imin];
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
		cell = bd.get_cell(i,j,k);
		IFace = cell->iface[WEST];
		FlowState &fs = *(IFace->fs);
		fs.copy_values_from(*(cell->fs));
		fs.vel.x = 0.0; fs.vel.y = 0.0; fs.vel.z = 0.0;
		for ( size_t imode=0; imode < nmodes; ++imode ) {
		    fs.gas->T[imode] = TProfile[j-bd.jmin];
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

// perfect gas equation for calculating the wall (ghost cell) density

#define GHOST_DENSITY_QUADRATIC_SOLUTION \
 ( cell_momentum_flux + sqrt( cell_momentum_flux*cell_momentum_flux - \
     4.0 * R * T * pow( 2.0 * mdot_total - cell_mass_flux, 2 ) ) ) / ( 2.0 * R * T )


// calculate_ghost_cell_flow_state function definition

int AblatingBC::
calculate_ghost_cell_flow_state(FV_Cell *cell1, FV_Interface *wall, FV_Cell *cell0)
{
    // Perform Newton-Raphson iterations to solve for ghost cell flow state
    // based on mass and momentum balance
    // cell1 -> CFD cell bounding wall
    //  wall -> ablating wall interface element
    // cell0 -> first ghost cell
    
    // 0. Set variables that are constant for the calculation
    // NOTE: - using CFD cell density as initial guess
    //       - setting ghost cell velocity vector to transformed cfd cell values
    //         then overwriting x component from calculation, then transforming back
    cell_rho = cell1->fs->gas->rho;
    cell1->fs->vel.transform_to_local(wall->n, wall->t1, wall->t2);
    cell_un = cell1->fs->vel.x;
    cell_local_vel = cell1->fs->vel;  // FIX-ME copy components instead ??
    cell1->fs->vel.transform_to_global(wall->n, wall->t1, wall->t2);
    cell_mass_flux = cell1->fs->gas->rho * cell_un;
    cell_momentum_flux = cell1->fs->gas->p + cell1->fs->gas->rho * cell_un * cell_un;
    for ( size_t isp=0; isp<cell1->fs->gas->massf.size(); ++isp )
    	cell_massf[isp] = cell1->fs->gas->massf[isp];
    // Ghost-cell temperatures for working gas-data structure
    for ( size_t itm=0; itm<Q->T.size(); ++itm ) {
        Q->T[itm] = cell1->fs->gas->T[itm];
        //Q->T[itm] = 2.0 * Twall - cell1->fs->gas->T[itm];
    }
    // y_guess
    if ( cell0->fs->gas->rho > 0.0 ) {
    	// just use the current ghost cell state as the guess
    	size_t iy;
    	// species mass densities
    	for ( iy=0; iy<Q->massf.size(); ++iy )
    	    y_guess[iy] = cell0->fs->gas->massf[iy] * cell0->fs->gas->rho;
    	// normal velocity 
    	cell0->fs->vel.transform_to_local(wall->n, wall->t1, wall->t2);
    	y_guess[iy] = cell0->fs->vel.x;
    	// cout << "using existing solution for initial guess" << endl;
    }
    else {
    	// This is probably the first CFD step - we need an educated guess
    	// NOTE: it is worth putting even a lot of effort here as this is only 
    	//       used at simulation start and the non-linear system is very sensitive
        for ( size_t isp=0; isp<cell1->fs->gas->massf.size(); ++isp )
    	    Q->massf[isp] = cell1->fs->gas->massf[isp];	// use CFD cell mass-fracs for first guess
    	double R = gmodel->R(*Q);
    	double T = Q->T[0];
    	double R_new, rho_0, u_0;
    	for ( size_t i=0; i<10; ++i ) {
            //previous
            //rho_0 = GHOST_DENSITY_QUADRATIC_SOLUTION;
            rho_0 = fabs(mdot_total)/cell_rho;

           // from eq. 17, solved for v0 - not sure about this...
            u_0 = ( 2.0 * mdot_total - cell_mass_flux ) / rho_0;
           // cout << "rho_0 = " << rho_0 << ", u_0 = " << u_0 << endl;

    	    // Reculate mass-fractions and R
    	    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	    	Q->massf[isp] = ( 2.0 * mdot[isp] - cell1->fs->gas->massf[isp] * cell_un ) / ( rho_0 * u_0 );
    	    R_new = gmodel->R(*Q);
    	    if ( ( R_new - R ) / R_new < 1.0e-6 ) break;
    	    R = R_new;
    
    	}
        if ( rho_0 < 0.0 ) {
            cout << "rho in ghost cell is negative!" << endl
            << "cell position: " << endl
            << "(" <<  cell0->pos[0].x << "," << cell0->pos[0].y << "," <<  cell0->pos[0].z << ")" << endl
	     << "Bailing out!" << endl;
	    exit( FAILURE );

            //previous
	    // cout << "AblatingBC::calculate_ghost_cell_flow_state()" << endl
	    //	 << "No physical solution exists" << endl;
            // use approx solution with CFD cell mass-fracs
            for ( size_t isp=0; isp<cell1->fs->gas->massf.size(); ++isp )
            	Q->massf[isp] = cell1->fs->gas->massf[isp];	// use CFD cell mass-fracs for first guess
            R = gmodel->R(*Q);
            T = Q->T[0];
            rho_0 = GHOST_DENSITY_QUADRATIC_SOLUTION;
            u_0 = ( 2.0 * mdot_total - cell_mass_flux ) / rho_0;
            // cout << "rho_0 = " << rho_0 << ", u_0 = " << u_0 << endl;
	}
    	size_t iy;
    	// species mass densities
    	for ( iy=0; iy<Q->massf.size(); ++iy )
    	    y_guess[iy] = Q->massf[iy] * rho_0;
    	// normal velocity 
    	y_guess[iy] = u_0;
    	// cout << "calculating an initial guess" << endl;
    	// cout << "mdot_total = " << mdot_total << ", rho_0 = " << rho_0 << endl;
    	// cout << "cell_un = " << cell_un << ", cell_mass_flux = " << cell_mass_flux << ", cell_momentum_flux = " << cell_momentum_flux << endl;
    }
    
    //cout << "CFD cell:" << " rho = " << cell_rho << ", un = " << cell_un << ", p = " << cell1->fs->gas->p << ", T = " << cell1->fs->gas.T[0] << endl;
    //cout << "cell_mass_flux = " << cell_mass_flux << endl;
    //cout << "cell_momentum_flux = " << cell_momentum_flux << endl;
    //cout << "mdot[0] = " << mdot[0] << ", mdot[1] = " << mdot[1] << endl;
    
    // 1. Solve the system
    if ( zero_solver->solve( *this, y_guess, y_out ) ) {
	cout << "AblatingBC::calculate_ghost_cell_flow_state()" << endl
	     << "Zero solver has failed, bailing out!" << endl;
	exit( FAILURE );
    }
    // 2.  map results
    //  a. temperatures, density and mass-fractions
    for ( size_t itm=0; itm<Q->T.size(); ++itm )
    	cell0->fs->gas->T[itm] = Q->T[itm];
    cell0->fs->gas->rho = 0.0;
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	cell0->fs->gas->rho += y_out[isp];
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	cell0->fs->gas->massf[isp] = y_out[isp] / cell0->fs->gas->rho;
    //  b. eval thermo state
    gmodel->eval_thermo_state_rhoT(*(cell0->fs->gas));
    //  d. non-normal and normal velocity component
    cell0->fs->vel = cell_local_vel; // FIX-ME, maybe copy components
    cell0->fs->vel.x = y_out[u0_index];
    
    valarray<double> G; G.resize( y_out.size() );
    f(y_out,G);
    //just output for checking
    //    for ( size_t i=0; i<G.size(); ++i )
    //    	cout << "value of the function after inserting the root = G[" << i << "] = " << G[i] << endl;
    Valmatrix dGdy;
    dGdy.resize( y_out.size(), y_out.size() );
    Jac(y_out,dGdy);
    //just output for checking
    //cout << "Jacobian" << dGdy.str() << endl;
    
    //  e. transform back to global coordinate system
    // cout << "A. cell0->fs->vel.x = " << cell0->fs->vel.x << endl;
    cell0->fs->vel.transform_to_global(wall->n, wall->t1, wall->t2);
    // cout << "B. cell0->fs->vel.x = " << cell0->fs->vel.x << endl;
    
    // 3. Check that values are permitted
    if ( !cell0->check_flow_data() ) {
    	cout << "AblatingBC::calculate_ghost_cell_flow_state()" << endl
    	     << "Computed ghost cell flow state contains bad data." << endl
    	     << "Exiting program." << endl;
    	exit( NUMERICAL_ERROR );
    }
    
    // cell0->print(1);
    // exit(1);
    // cout << "Success!" << endl;

    return SUCCESS;
}

// some constant definitions

const size_t WITH_TOTAL_MASS_CONSERVATION = 1;
const int NORMALISE_G = 1;


// creation of source terms N, following the formulation in Beerman et al
// and Chen & Milos.
// Structure of Jacobian based on that in chemical-equilibrium-system.cxx.

int AblatingBC::compute_source_terms(vector<double> &massf)
{
    // Initialise N vector with zeros
    for ( size_t isp=0; isp<nsp; ++isp ) {
	 N[isp] = 0.0;
    }

    //int iQ=0;	// current element index for the Q valarray


    return SUCCESS;
}

// 'f' function definition - creating the non-linear system
// of equations to be solved using the Newton-Raphson method.
// Structure of f based on that in chemical-equilibrium-system.cxx.

int AblatingBC::f(const valarray<double> &y, valarray<double> &G)
{

    /* Create the equation system for the ZeroSystem for a given y vector */

    // 0.  Apply source terms (as negatives as we are creating a zero system)
    for ( size_t isp=0; isp<nsp; ++isp ) {
    }

    // 0. unpack the y valarray
    double u0 = y[u0_index];
    //double rho_wall = y[rho_index];
    Q->rho = 0.0;
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->rho += y[isp];
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->massf[isp] = y[isp] / Q->rho;
    gmodel->eval_thermo_state_rhoT(*Q); 
    
    size_t iG=0;
    
    // 1. species mass conservation
    for ( iG=0; iG<cell_massf.size(); ++iG ) {
    	G[iG] = 0.5 * ( cell_mass_flux * cell_massf[iG] + y[iG] * u0 ) - mdot[iG];
    }
    
    // 2. total momentum conservation
    G[iG] = cell_momentum_flux - Q->p - Q->rho * u0 * u0;


// TESTING ONLY
// for ( size_t iG=0; iG<G.size(); ++iG )
//    cout << "G[" << iG << "] = " << G[iG] << endl;
// cout << "u0 = " << u0 << endl;
// Q.print_values();
    	
    return SUCCESS;
}


// 'Jac' function definition - create the Jacobian for use in
// N-R iterations. Manual definition required as specified by
// zero_finders files.
// Structure of Jacobian based on that in chemical-equilibrium-system.cxx.

int AblatingBC::Jac(const valarray<double> &y, Valmatrix &dGdy)
{
	/* Create the Jacobian matrix for the ZeroSystem for a given y vector */

	// 0.  Clear the jacobian matrix
	for ( size_t i=0; i<nsp; ++i ) {
		for ( size_t j=0; j<nsp; ++j ) {
		    dGdy.set(i,j,0.0);
	}
	}

	int iG = 0;		// current matrix line

    // 0. unpack the y valarray
    double u0 = y[u0_index];
    Q->rho = 0.0;
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->rho += y[isp];
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->massf[isp] = y[isp] / Q->rho;
    gmodel->eval_thermo_state_rhoT(*Q);


    // 1. species mass conservation derivaties

    
    // 2. total momentum conservation derivatives

// cout << Gdy.str() << endl;
// cout << "Gdy.det() = " << eval_matrix_determinant(Gdy) << endl;
    
    return SUCCESS;
}
