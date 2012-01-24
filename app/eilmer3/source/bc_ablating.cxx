// bc_ablating.cxx

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

AblatingBC::AblatingBC( Block &bdp, int which_boundary, double Twall, 
			vector<double> &mdot, const std::string filename )
    : BoundaryCondition(bdp, which_boundary, ABLATING, "AblatingBC",
			true, false, -1, -1, 0),
      Twall(Twall), mdot(mdot), filename(filename), max_iterations(1000000), tol(1.0e-6)
{
    // Reads the temperature profile from the solid solver output.

    char line[512], token[512];
    double T;
    FILE *fp;
    int ncell, nread;
    bool with_input_T_file = 1;

    // Get number of cells for this boundary
    if ( which_boundary == NORTH || which_boundary == SOUTH ) {
	    ncell = bdp.nni;
    } else {
	    ncell = bdp.nnj;
    }
    // Use fixed boundary temperature if there is no input file
    fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "AblatingBC() constructor:"
	    << " cannot open file " << filename << endl;
        cerr << " Using fixed wall temperature = " << Twall << endl; 
        with_input_T_file = 0;        
    }

    if ( with_input_T_file == 1 ) {
        if ( fgets(line, sizeof(line), fp) == NULL ) {
	    cerr << "AblatingBC(): failure of fgets()" << endl;
	    cerr << "Quitting program." << endl;
	    exit(FILE_ERROR);
        }

        nread = sscanf(line, "%d", &ncell_for_profile);
        if ( nread != 1 ) {
            cerr << "AblatingBC() constructor:"
	        << "Could not read ncell_for_profile from line:" << endl
	        << line;
            exit(BAD_INPUT_ERROR);
        }
    
        // Check for that the number of cells is appropriate for this boundary
    
        if ( ncell != ncell_for_profile ) {
            cerr << "AblatingBC() constructor:" << endl
	        << "    Inconsistent numbers of cells: ncell=" << ncell
	        << ", ncell_for_profile=" << ncell_for_profile << endl;
            exit(BAD_INPUT_ERROR);
        }
    }

    /* For each line in the file, store the temperature data. */
    for ( int ii = 0; ii < ncell; ++ii ) {
        if ( with_input_T_file == 1 ) {
	        if ( fgets(line, sizeof(line), fp) == NULL ) {
	            cerr << "AblatingBC(): failure of fgets()" << endl;
	            cerr << "Quitting program." << endl;
	            exit(FILE_ERROR);
	        }
            /* Pull the line apart with the string tokenizer. */
            strcpy( token, strtok(line, " ") );
            sscanf( token, "%lf", &T );
        } else {
            T = Twall;
        }
        TProfile.push_back(T);
    } // end for
    
    // 0. Get gas-model pointer
    gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    // 1. Calculate the total mass flux from the given species-specific components
    mdot_total = 0.0;
    
    for ( size_t isp=0; isp<mdot.size(); ++isp )
    	mdot_total += mdot[isp];
    // 2. Initialise the local gas-data structure (used for EOS calls)
    Q = new Gas_data(gmodel);
    // 3. Size the CFD cell mass-fraction vector
    cell_massf.resize(nsp);
    // 4. initialise the zero system components
    u0_index = nsp;
    y_guess.resize(nsp+1);
    y_out.resize(nsp+1);
    zero_solver = new NewtonRaphsonZF(nsp+1, tol, max_iterations, true);
    // f_jac is the Jacobian scale factor
    f_jac = 1.0;
}

AblatingBC::AblatingBC( const AblatingBC &bc )
    : BoundaryCondition(bc.bdp, bc.which_boundary, bc.type_code, bc.name_of_BC,
			bc.is_wall_flag, bc.use_udf_flux_flag,
			bc.neighbour_block, bc.neighbour_face,
			bc.neighbour_orientation),
      Twall(bc.Twall), mdot(bc.mdot),
      filename(bc.filename), mdot_total(bc.mdot_total),
      gmodel(bc.gmodel), u0_index(bc.u0_index), 
      max_iterations(bc.max_iterations), tol(bc.tol), f_jac(bc.f_jac)
{
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

AblatingBC::~AblatingBC()
{
    delete zero_solver;
    delete Q;
}

int AblatingBC::apply_inviscid( double t )
{
    // Calculate the ghost cell flow-states from the given wall temperature
    // and species-specific mass flux.
    int i, j, k;
    FV_Cell *src_cell, *dest_cell;
    FV_Interface *IFace;

    switch ( which_boundary ) {
    case NORTH:
	j = bdp.jmax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		src_cell = bdp.get_cell(i,j,k);
		IFace = src_cell->iface[NORTH];
		// ghost cell 1
		dest_cell = bdp.get_cell(i,j+1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bdp.get_cell(i,j+1,k);
		dest_cell = bdp.get_cell(i,j+2,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
	    } // end i loop
	} // for k
	break;
    case EAST:
	i = bdp.imax;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		IFace = src_cell->iface[EAST];
		// ghost cell 1
		dest_cell = bdp.get_cell(i+1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bdp.get_cell(i+1,j,k);
		dest_cell = bdp.get_cell(i+2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
	    } // end j loop
	} // for k
	break;
    case SOUTH:
	j = bdp.jmin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (i = bdp.imin; i <= bdp.imax; ++i) {
		src_cell = bdp.get_cell(i,j,k);
		IFace = src_cell->iface[SOUTH];
		// ghost cell 1
		dest_cell = bdp.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bdp.get_cell(i,j-1,k);
		dest_cell = bdp.get_cell(i,j-1,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
	    } // end i loop
	} // for k
	break;
    case WEST:
	i = bdp.imin;
        for (k = bdp.kmin; k <= bdp.kmax; ++k) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		// ghost cell 1
		dest_cell = bdp.get_cell(i-1,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bdp.get_cell(i-1,j,k);
		dest_cell = bdp.get_cell(i-2,j,k);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
	    } // end j loop
	} // for k
 	break;
    case TOP:
	k = bdp.kmax;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		// ghost cell 1
		dest_cell = bdp.get_cell(i,j,k+1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bdp.get_cell(i,j,k+1);
		dest_cell = bdp.get_cell(i,j,k+2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
	    } // end j loop
	} // for i
	break;
    case BOTTOM:
	k = bdp.kmin;
        for (i = bdp.imin; i <= bdp.imax; ++i) {
	    for (j = bdp.jmin; j <= bdp.jmax; ++j) {
		src_cell = bdp.get_cell(i,j,k);
		IFace = src_cell->iface[WEST];
		// ghost cell 1
		dest_cell = bdp.get_cell(i,j,k-1);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
		calculate_ghost_cell_flow_state( src_cell, IFace, dest_cell );
		// ghost cell 2 (copy of 1)
		src_cell = bdp.get_cell(i,j,k-1);
		dest_cell = bdp.get_cell(i,j,k-2);
		dest_cell->copy_values_from(*src_cell, COPY_FLOW_STATE);
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

int AblatingBC::apply_viscous( double t )
{
    int i, j, k;
    FV_Cell *cell;
    FV_Interface *IFace;
    int nmodes = gmodel->get_number_of_modes();

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
		for ( int imode=0; imode < nmodes; ++imode ) { 
		    fs.gas->T[imode] = TProfile[i-bdp.imin];
	    }
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bdp.bcp[NORTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
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
		for ( int imode=0; imode < nmodes; ++imode ) {
		    fs.gas->T[imode] = TProfile[j-bdp.jmin];
		}
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bdp.bcp[EAST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
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
		for ( int imode=0; imode < nmodes; ++imode ) {
		    fs.gas->T[imode] = TProfile[i-bdp.imin];
	    }
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bdp.bcp[SOUTH]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end i loop
	} // for k
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
		for ( int imode=0; imode < nmodes; ++imode ) {
		    fs.gas->T[imode] = TProfile[j-bdp.jmin];
	    }
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bdp.bcp[WEST]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for k
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
		for ( int imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bdp.bcp[TOP]->wc_bc != NON_CATALYTIC) {
		    cw->apply(*(cell->fs->gas), fs.gas->massf);
		}
	    } // end j loop
	} // for i
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
		for ( int imode=0; imode < nmodes; ++imode ) fs.gas->T[imode] = Twall;
		fs.tke = 0.0;
		fs.omega = ideal_omega_at_wall(cell);
		if (bdp.bcp[BOTTOM]->wc_bc != NON_CATALYTIC) {
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

#define GHOST_DENSITY_QUADRATIC_SOLUTION \
 ( cell_momentum_flux + sqrt( cell_momentum_flux*cell_momentum_flux - \
     4.0 * R * T * pow( 2.0 * mdot_total - cell_mass_flux, 2 ) ) ) / ( 2.0 * R * T )



int AblatingBC::
calculate_ghost_cell_flow_state( FV_Cell *cell1, FV_Interface *wall, FV_Cell *cell0 )
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
    	int iy;
    	// species mass densities
    	for ( iy=0; iy<(int)Q->massf.size(); ++iy )
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
    	for ( int i=0; i<10; ++i ) {
            //previous
            //rho_0 = GHOST_DENSITY_QUADRATIC_SOLUTION;
            rho_0 = fabs(mdot_total)/cell_rho;

           // from eq. 17, solved for v0
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
            << "(" <<  cell0->pos.x << "," << cell0->pos.y << "," <<  cell0->pos.z << ")" << endl
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
    	int iy;
    	// species mass densities
    	for ( iy=0; iy<(int)Q->massf.size(); ++iy )
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
    
    // double ghost_cell_momentum_flux = cell0->fs->gas->p + cell0->fs->gas->rho * cell0->fs->vel.x * cell0->fs->vel.x;
    // cout << "cell_momentum_flux = " << cell_momentum_flux << ", ghost_cell_momentum_flux = " << ghost_cell_momentum_flux << endl;
    // double ghost_cell_mass_flux = cell0->fs->gas->rho * cell0->fs->vel.x;
    //just output for checking
    //cout << "massf[0]: 05*(ghost_cell_mass_flux + cell_mass_flux) = " << 0.5*(cell0->fs->gas->massf[0]*ghost_cell_mass_flux + cell1->fs->gas->massf[0]*cell_mass_flux) << ", mdot = " << mdot[0] << endl;
    //cout << "massf[1]: 05*(ghost_cell_mass_flux + cell_mass_flux) = " << 0.5*(cell0->fs->gas->massf[1]*ghost_cell_mass_flux + cell1->fs->gas->massf[1]*cell_mass_flux) << ", mdot = " << mdot[1] << endl;
    //cout << "massf[2]: 05*(ghost_cell_mass_flux + cell_mass_flux) = " << 0.5*(cell0->fs->gas->massf[2]*ghost_cell_mass_flux + cell1->fs->gas->massf[2]*cell_mass_flux) << ", mdot = " << mdot[2] << endl;
    
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

#define WITH_TOTAL_MASS_CONSERVATION 1
#define NORMALISE_G 1

int AblatingBC::f( const valarray<double> &y, valarray<double> &G )
{
    // 0. unpack the y valarray
    double u0 = y[u0_index];
    Q->rho = 0.0;
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->rho += y[isp];
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->massf[isp] = y[isp] / Q->rho;
    gmodel->eval_thermo_state_rhoT(*Q); 
    
#   if NORMALISE_G
    double mass_flux_norm = fabs(mdot_total) + fabs(cell_mass_flux) + fabs(u0*Q->rho);
    double momentum_flux_norm = cell_momentum_flux + Q->p + Q->rho*u0*u0;
#   else
    double mass_flux_norm = 1.0;
    double momentum_flux_norm = 1.0;
#   endif
    
    int iG;
    
#   if WITH_TOTAL_MASS_CONSERVATION
    // 1a. species mass conservation
    for ( iG=0; iG<(int)cell_massf.size()-1; ++iG ) {
    	G[iG] = 0.5 * ( cell_mass_flux * cell_massf[iG] + y[iG] * u0 ) - mdot[iG];
    	// Normalize
    	G[iG] /= mass_flux_norm;
    }
    // 1b. Total mass conservation
    G[iG] = 0.5 * ( cell_mass_flux + Q->rho * u0 ) - mdot_total;
    // Normalize
    G[iG] /= mass_flux_norm;
    
    ++iG;
#   else
    // 1. species mass conservation
    for ( iG=0; iG<(int)cell_massf.size(); ++iG ) {
    	G[iG] = 0.5 * ( cell_mass_flux * cell_massf[iG] + y[iG] * u0 ) - mdot[iG];
    	// Normalize
    	G[iG] /= mass_flux_norm;
    }
#   endif
    
    // 2. total momentum conservation
    G[iG] = cell_momentum_flux - Q->p - Q->rho * u0 * u0;
    // Normalize
    G[iG] /= momentum_flux_norm;

// TESTING ONLY
// for ( size_t iG=0; iG<G.size(); ++iG )
//    cout << "G[" << iG << "] = " << G[iG] << endl;
// cout << "u0 = " << u0 << endl;
// Q.print_values();
    	
    return SUCCESS;
}

int AblatingBC::Jac( const valarray<double> &y, Valmatrix &Gdy )
{
    // 0. unpack the y valarray
    double u0 = y[u0_index];
    Q->rho = 0.0;
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->rho += y[isp];
    for ( size_t isp=0; isp<Q->massf.size(); ++isp )
    	Q->massf[isp] = y[isp] / Q->rho;
    gmodel->eval_thermo_state_rhoT(*Q);
    
#   if NORMALISE_G
    double mass_flux_norm = fabs(mdot_total) + fabs(cell_mass_flux) + fabs(u0*Q->rho);
    double momentum_flux_norm = cell_momentum_flux + Q->p + Q->rho*u0*u0;
#   else
    double mass_flux_norm = 1.0;
    double momentum_flux_norm = 1.0;
#   endif

    int iG;

#if WITH_TOTAL_MASS_CONSERVATION
    // 1a. species mass conservation derivaties
    for ( iG=0; iG<(int)cell_massf.size()-1; ++iG ) {
    	//  a. wrt species mass densities (off-diagonals are zero)
    	Gdy.set(iG,iG,0.5*u0/mass_flux_norm*f_jac);
    	//  b. wrt velocity
    	Gdy.set(iG,u0_index,0.5*y[iG]/mass_flux_norm*f_jac);
    }
    
    // 1a. total mass conservation derivaties
    for ( int jsp=0; jsp<(int)cell_massf.size(); ++jsp ) {
    	//  a. wrt species mass densities (off-diagonals are zero)
    	Gdy.set(iG,jsp,0.5*u0/mass_flux_norm*f_jac);
    }
    //  b. wrt velocity
    Gdy.set(iG,u0_index,0.5*Q->rho/mass_flux_norm*f_jac);
    ++iG;
#else
    // 1. species mass conservation derivaties
    for ( iG=0; iG<(int)cell_massf.size(); ++iG ) {
    	//  a. wrt species mass densities (off-diagonals are zero)
    	Gdy.set(iG,iG,0.5*u0/mass_flux_norm*f_jac);
    	//  b. wrt velocity
    	Gdy.set(iG,u0_index,0.5*y[iG]/mass_flux_norm*f_jac);
    }
#endif
    
    // 2. total momentum conservation derivatives
    //  a. wrt species mass densities
    int status;
    for ( size_t jsp=0; jsp<cell_massf.size(); ++jsp ) {
    	Gdy.set(iG,jsp,(-gmodel->dpdrho_i_const_T(*Q,jsp,status)-u0*u0)/momentum_flux_norm*f_jac);
    }
    //  b. wrt velocity
    Gdy.set(iG,u0_index,-2.0*Q->rho*u0/momentum_flux_norm*f_jac);
    
// cout << Gdy.str() << endl;
// cout << "Gdy.det() = " << eval_matrix_determinant(Gdy) << endl;
    
    return SUCCESS;
}
