/** \file implicit.cxx
 * \ingroup eilmer3
 * \brief Functions to compute point-implicit viscous/inviscid updates for Eilmer3.
 *
 * \author OJ / ESIL 
 * \version October 2009 - May 2010
 * \version June 2010 updated by PJ to use the new data structures.
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "block.hh"
#include "implicit.hh"
#include "kernel.hh"
#include "exch2d.hh"
#include "visc.hh"
#include "visc3D.hh"
#include "init.hh"
#include "bc.hh"
#include "cell.hh"
#include "diffusion.hh"
#include "main.hh"

const double VERY_SMALL = 1.0e-10;


int gasdynamic_point_implicit_inviscid_increment(void)
{
#if WITH_IMPLICIT == 1
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    int most_bad_cells;
    int attempt_number, step_failed;
    using std::swap;
	
    cout << "=== Debut gasdynamic_point_implicit_inviscid_increment ===" << endl;
    // Record the current values of the conserved variables
    // in preparation for applying the predictor and corrector
    // stages of the time step.
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
        if ( bdp->active != 1 ) continue;
	for ( FV_Cell *cp: bdp->active_cells ) cp->record_conserved();
    }

    attempt_number = 0;
    do {
#       ifdef _MPI
	MPI_Barrier( MPI_COMM_WORLD );
#       endif
	++attempt_number;
	step_failed = 0;

	//  Predictor Stage for gas-dynamics
#       ifdef _MPI
	mpi_exchange_boundary_data(G.my_mpi_rank, COPY_FLOW_STATE);
#       else
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( bdp->active != 1 ) continue;
	    exchange_shared_boundary_data( jb, COPY_FLOW_STATE );
	}
#       endif
	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( bdp->active != 1 ) continue;
	    apply_inviscid_bc( *bdp, G.sim_time, G.dimensions );
	}

	if ( get_flux_calculator() == FLUX_ADAPTIVE ) {
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		bdp->detect_shock_points( G.dimensions );
	    }
	}
	
	// Non-local radiation transport needs to be performed a-priori for parallelization.
	// Note that Q_rad is not re-evaluated for corrector step.
	if ( get_radiation_flag() ) {
	    RadiationTransportModel * rtm = get_radiation_transport_model_ptr();
	    global_data &G = *get_global_data_ptr();
	    Block * bdp;
	    
	    // Determine if a scaled or complete radiation call is required
	    if ( ( (G.step / get_radiation_update_frequency()) * 
		   get_radiation_update_frequency() == G.step) ) {
		// recompute
		rtm->compute_Q_rad_for_flowfield();
		// store the radiation scaling parameters for each cell
#		ifdef _OPENMP
#		pragma omp parallel for private(jb) schedule(runtime)
#		endif
		for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		    bdp = G.my_blocks[jb];
		    if ( bdp->active != 1 ) continue;
		    for ( FV_Cell *cp: bdp->active_cells ) cp->store_rad_scaling_params();
		}
	    }
	    else {
		// rescale
#		ifdef _OPENMP
#		pragma omp parallel for private(jb) schedule(runtime)
#		endif
		for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		    bdp = G.my_blocks[jb];
		    if ( bdp->active != 1 ) continue;
		    for ( FV_Cell *cp: bdp->active_cells ) cp->rescale_Q_rE_rad();
		}
	    }
	}

	for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	    bdp = G.my_blocks[jb];
	    if ( bdp->active != 1 ) continue;
	    bdp->inviscid_flux( G.dimensions );
	    for ( FV_Cell *cp: bdp->active_cells ) {
		cp->inviscid_source_vector(bdp->omegaz);
		if ( G.udf_source_vector_flag == 1 ) 
		    cp->udf_source_vector_for_cell(G.dt_global);
		cp->inviscid_point_implicit_update_for_cell();
		swap(cp->U[0], cp->U[1]); 
		cp->decode_conserved(0, 1, bdp->omegaz);
	    }
	} // end of for jb...

	// 2d. Check the record of bad cells and if any cells are bad, 
	//     fail this attempt at taking a step,
	//     set everything back to the initial state and
	//     reduce the time step for the next attempt
	most_bad_cells = do_bad_cell_count();
	if ( ADJUST_INVALID_CELL_DATA == 0 && most_bad_cells > 0 ) {
	    step_failed = 1;
	    G.dt_global = G.dt_reduction_factor * G.dt_global;
	    printf( "Attempt %d failed: reducing dt to %e.\n",
		    attempt_number, G.dt_global);
	    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
		bdp = G.my_blocks[jb];
		if ( bdp->active != 1 ) continue;
		for ( FV_Cell *cp: bdp->active_cells ) cp->decode_conserved(0, 0, bdp->omegaz);
	    }
	}

    } while (attempt_number < 3 && step_failed == 1);
	
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	if ( bdp->active != 1 ) continue;
	for ( FV_Cell *cp: bdp->active_cells ) {
	    swap(cp->U[0], cp->U[1]); 
	}
    }
    cout << "=== Fin gasdynamic_point_implicit_inviscid_increment ===" << endl;
    return step_failed;
#else
    return 0;
#endif
} //int gasdynamic_point_implicit_inviscid_increment

int gasdynamic_fully_implicit_inviscid_increment(void)
{
    cout << "gasdynamic_fully_implicit_inviscid_increment()" << endl
	 << "Not implemented yet!" << endl;
	
    return FAILURE;
} //int gasdynamic_fully_implicit_inviscid_increment

int inviscid_point_implicit_update_for_cell(FV_Cell *cell)
{
#if WITH_IMPLICIT == 1
    int aa, bb;
    global_data &G = *get_global_data_ptr();  // set up a reference
    double dt_global;
    int dimensions;
    dimensions = G.dimensions;
    dt_global = G.dt_global;
    //cell->encode_conserved();
    //cell->decode_conserved();
    calculate_M_inviscid(cell, dimensions);	    
			
    calculate_h_inviscid(cell, dimensions); 
	
    for ( aa = 1; aa <= 5; ++aa ) {
	for ( bb = 1; bb <= 5; ++bb ) {
	    cell->piM[aa][bb] = dt_global / cell->volume * cell->piM[aa][bb];
	}
	cell->pir[aa][1] = -dt_global / cell->volume * cell->pir[aa][1];
    }	
		
    /* Now solve the equation set eq 5.1 */
    /* 1. Increment M by the identity matrix */
    cell->piM[1][1] += 1.0;
    cell->piM[2][2] += 1.0;
    cell->piM[3][3] += 1.0;
    cell->piM[4][4] += 1.0;
    cell->piM[5][5] += 1.0;
	
    /* 2. Solve for increment in flow properties */
    gaussj(cell, 5, 5);
	
    /* 3. Update cell centre convserved variables for next iteration */
    cell->U[1]->mass = cell->U[0]->mass + cell->pir[1][1];
    cell->U[1]->momentum.x = cell->U[0]->momentum.x + cell->pir[2][1];
    cell->U[1]->momentum.y = cell->U[0]->momentum.y + cell->pir[3][1];
    cell->U[1]->momentum.z = cell->U[0]->momentum.z + cell->pir[4][1];
    cell->U[1]->total_energy = cell->U[0]->total_energy + cell->pir[5][1];
	
    // Species densities: mass of species isp per unit volume.
    for ( size_t isp = 0; isp < cell->U[1]->massf.size(); ++isp ) {
	cell->U[1]->massf[isp] = cell->U[1]->mass * cell->fs->gas->massf[isp];
    }
#endif		
    return SUCCESS;
} //int inviscid_point_implicit_update_for_cell

int calculate_M_inviscid(FV_Cell *cell, int dimensions) 
{
#if WITH_IMPLICIT == 1
    int aa, bb;
    FV_Interface *IFn = cell->iface[NORTH];
    FV_Interface *IFe = cell->iface[EAST];
    FV_Interface *IFs = cell->iface[SOUTH];
    FV_Interface *IFw = cell->iface[WEST];
    FV_Interface *IFt = cell->iface[TOP];
    FV_Interface *IFb = cell->iface[BOTTOM];
	
    calculate_inviscid_jacobian(cell, IFn);
    calculate_inviscid_jacobian(cell, IFe);
    calculate_inviscid_jacobian(cell, IFs);
    calculate_inviscid_jacobian(cell, IFw);
    if (dimensions == 3) {
	calculate_inviscid_jacobian(cell, IFt);
	calculate_inviscid_jacobian(cell, IFb);
    }

    /* 	Initialize matrix piM to zero */
    for (aa = 1; aa <=5; ++aa ) {
	for (bb = 1; bb <=5; ++bb ) {
	    cell->piM[aa][bb] = 0.0;
	}
    }
	
    // Calculate M_inviscid as per equation .....4.18
    for ( aa = 1; aa <= 5; ++aa ) {
	for ( bb = 1; bb <= 5; ++bb ) {
	    cell->piM[aa][bb] = 0.5 * (abs(IFn->J[aa][bb])*IFn->area + abs(IFs->J[aa][bb])*IFs->area +
				       abs(IFe->J[aa][bb])*IFe->area + abs(IFw->J[aa][bb])*IFw->area);
	    if (dimensions == 3) 
		cell->piM[aa][bb] += 0.5 * (abs(IFt->J[aa][bb])*IFt->area + abs(IFb->J[aa][bb])*IFb->area);
	}	
    }
#endif
    return SUCCESS;
} //int calculate_M_inviscid

int calculate_inviscid_jacobian(FV_Cell *cell, FV_Interface *iface)
{
#if WITH_IMPLICIT == 1
    Gas_model *gmodel = get_gas_model_ptr();
    int aa, bb, statusf;
    double rho, gamma;
    //double mu_eff, mu_lam, mu_t, Pr, e_int;
    Vector3 n;
    double alpha, beta, H;
    double Iu,Iv,Iw;
    //double Cu,Cv,Cw;	
    double Ul;
    //double viscous_factor; // so that we can scale down the viscous effects
    //double factor;
	
    Iu = iface->fs->vel.x;
    Iv = iface->fs->vel.y;
    Iw = iface->fs->vel.z;
     
    //Cu = cell->fs->vel.x;
    //Cv = cell->fs->vel.y;
    //Cw = cell->fs->vel.z;
    n = iface->n;
    //UL = Cu*n.x + Cv*n.y + Cw*n.z;
    gamma = gmodel->gamma(*(iface->fs->gas), statusf);
    rho = cell->fs->gas->rho;
    //Pr = iface->Prandtl;

    alpha = 0.5 * (Iu*Iu + Iv*Iv + Iw*Iw);
    beta = gamma -1;
    Ul = Iu * n.x + Iv * n.y + Iw * n.z;
    //H = cell->fs->gas->a * cell->fs->gas->a / beta + alpha;
    H = gmodel->Cp(*(iface->fs->gas), statusf) * iface->fs->gas->T[0] + alpha;


    //viscous_factor = get_viscous_factor();

    //mu_lam = viscous_factor * iface->mu;
    //mu_t = viscous_factor * iface->mu_t;
    //mu_eff = mu_lam + mu_t;
    //e_int = cell->rE - dot(cell->vel, cell->vel);
    	
    /* 	Initalise jacobian matrix J */
    
    for ( aa = 1; aa <=5; ++aa ) {
	for ( bb = 1; bb <=5; ++bb ) {
	    iface->J[aa][bb] = 0.0;
	}
    }
    
    /* 	Fill jacobian matrix as of Appendix B */
    
    iface->J[1][1] = 0.0;
    iface->J[2][1] = alpha * beta * n.x - Ul * Iu;
    iface->J[3][1] = alpha * beta * n.y - Ul * Iv;
    iface->J[4][1] = alpha * beta * n.z - Ul * Iw;
    iface->J[5][1] = alpha * beta * Ul - Ul * H;
    
    iface->J[1][2] = n.x;
    iface->J[2][2] = -beta * Iu * n.x + Iu * n.x + Ul;
    iface->J[3][2] = -beta * Iu * n.y + Iv * n.x;
    iface->J[4][2] = -beta * Iu * n.z + Iw * n.x;
    iface->J[5][2] = -beta * Iu * Ul + H * n.x;
    
    iface->J[1][3] = n.y;
    iface->J[2][3] = -beta * Iv * n.x + Iu * n.y;
    iface->J[3][3] = -beta * Iv * n.y + Iv * n.y + Ul; 
    iface->J[4][3] = -beta * Iv * n.z + Iw * n.y;
    iface->J[5][3] = -beta * Iv * Ul + H * n.y;
    
    iface->J[1][4] = n.z;
    iface->J[2][4] = -beta * Iw * n.x + Iu * n.z; 
    iface->J[3][4] = -beta * Iw * n.y + Iv * n.z;
    iface->J[4][4] = -beta * Iw * n.z + Iw * n.z + Ul;
    iface->J[5][4] = -beta * Iw * Ul + H * n.z;
    
    iface->J[1][5] = 0.0;
    iface->J[2][5] = beta * n.x;
    iface->J[3][5] = beta * n.y;
    iface->J[4][5] = beta * n.z;
    iface->J[5][5] = beta * Ul + Ul;
	
#if 0
    // Note: Presently cut out.
    factor = mu_eff * ( iface->area / cell->volume ) / rho;

    for ( aa = 1; aa <= 5; ++aa ) {
	for ( bb = 1; bb <= 5; ++bb ) {
	    iface->J[aa][bb] = factor * iface->J[aa][bb];
	}
    }
#endif
#endif	
    return SUCCESS;
} //int calculate_inviscid_jacobian

int calculate_h_inviscid(FV_Cell *cell, int dimensions)
{
#if WITH_IMPLICIT == 1
    FV_Interface *IFn = cell->iface[NORTH];
    FV_Interface *IFe = cell->iface[EAST];
    FV_Interface *IFs = cell->iface[SOUTH];
    FV_Interface *IFw = cell->iface[WEST];
    FV_Interface *IFt = cell->iface[TOP];
    FV_Interface *IFb = cell->iface[BOTTOM];
	
    /* initialize r to zero */
    for (int aa = 1; aa <= 5; ++aa ) {
	cell->pir[aa][1] = 0.0;
    }
    
    /* Calculates the components of h as per equation 4.11 */
    cell->pir[1][1] = IFe->F->mass * IFe->area - IFw->F->mass * IFw->area +
	IFn->F->mass * IFn->area - IFs->F->mass * IFs->area;; 
    cell->pir[2][1] = IFe->F->momentum.x * IFe->area - IFw->F->momentum.x * IFw->area +
	IFn->F->momentum.x * IFn->area - IFs->F->momentum.x * IFs->area;
    cell->pir[3][1] = IFe->F->momentum.y * IFe->area - IFw->F->momentum.y * IFw->area +
	IFn->F->momentum.y * IFn->area - IFs->F->momentum.y * IFs->area;
    cell->pir[4][1] = 0.0;
    cell->pir[5][1] = IFe->F->total_energy * IFe->area - IFw->F->total_energy * IFw->area +
	IFn->F->total_energy * IFn->area - IFs->F->total_energy * IFs->area;
	
    if (dimensions == 3){
	cell->pir[1][1] += (IFt->F->mass * IFt->area - IFb->F->mass * IFb->area);
	cell->pir[2][1] += IFt->F->momentum.x * IFt->area - IFb->F->momentum.x * IFb->area;
	cell->pir[3][1] += IFt->F->momentum.y * IFt->area - IFb->F->momentum.y * IFb->area;
	cell->pir[4][1] += IFe->F->momentum.z * IFe->area - IFw->F->momentum.z * IFw->area +
	    IFn->F->momentum.z * IFn->area - IFs->F->momentum.z * IFs->area +
	    IFt->F->momentum.z * IFt->area - IFb->F->momentum.z * IFb->area;
	cell->pir[5][1] += IFt->F->total_energy * IFt->area - IFb->F->total_energy * IFb->area;
    }
#endif
    return SUCCESS;
} //int calculate_h_inviscid

int gasdynamic_point_implicit_viscous_increment(void)
{
#if WITH_IMPLICIT == 1
    global_data &G = *get_global_data_ptr();  // set up a reference
    Block *bdp;
    // Record the current values of the conserved variables
    // in preparation for applying the predictor and corrector
    // stages of the time step.
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	if ( bdp->active != 1 ) continue;
	for ( FV_Cell *cp: bdp->active_cells ) cp->record_conserved();
    }
    for ( int jb = 0; jb < G.my_blocks.size(); ++jb ) {
	bdp = G.my_blocks[jb];
	if ( bdp->active != 1 ) continue;
	bdp->clear_fluxes_of_conserved_quantities( G.dimensions );
	apply_viscous_bc( *bdp, G.sim_time, G.dimensions );
	if ( get_k_omega_flag() ) apply_menter_boundary_correction( *bdp );
	if ( G.dimensions == 2 ) {
	    viscous_derivatives_2D( bdp );
	} else {
	    viscous_derivatives_3D( bdp );
	}
	estimate_turbulence_viscosity( &G, bdp );
	if ( G.dimensions == 2 ) {
	    viscous_flux_2D( bdp );
	} else {
	    viscous_flux_3D( bdp );
	}
	for ( FV_Cell *cp: bdp->active_cells ) {
	    cp->viscous_source_vector();
	    point_implicit_update_for_cell(cp);
	    cp->decode_conserved(0, 1, bdp->omegaz);
	}
    } // end of for jb...
    cout << "=== Fin gasdynamic_point_implicit_viscous_increment ===" << endl;
#endif
    return SUCCESS;
} //int gasdynamic_point_implicit_viscous_increment

int gasdynamic_fully_implicit_viscous_increment(void)
{
    cout << "gasdynamic_fully_implicit_viscous_increment()" << endl
	 << "Not implemented yet!" << endl;
	
    return FAILURE;
} //int gasdynamic_fully_implicit_viscous_increment

int point_implicit_update_for_cell(FV_Cell *cell) 
{
#if WITH_IMPLICIT == 1
    int aa, bb;
    global_data &G = *get_global_data_ptr();  // set up a reference
    double dt_global;
    int dimensions;
    dimensions = G.dimensions;
    dt_global = G.dt_global;
    //cell->encode_conserved();
    //cell->decode_conserved();
    calculate_M(cell, dimensions);	    
			
    calculate_h(cell, dimensions); 
	
    for ( aa = 1; aa <= 5; ++aa ) {
	for ( bb = 1; bb <= 5; ++bb ) {
	    cell->piM[aa][bb] = dt_global / cell->volume * cell->piM[aa][bb];
	}
	cell->pir[aa][1] = -dt_global / cell->volume * cell->pir[aa][1];
    }	
		
    /* Now solve the equation set eq 5.1 */
    /* 1. Increment M by the identity matrix */
    cell->piM[1][1] += 1.0;
    cell->piM[2][2] += 1.0;
    cell->piM[3][3] += 1.0;
    cell->piM[4][4] += 1.0;
    cell->piM[5][5] += 1.0;
	
    /* 2. Solve for increment in flow properties */
    gaussj(cell, 5, 5);
	
    /* 3. Update cell centre convserved variables for next iteration */
    cell->U[1]->mass = cell->U[0]->mass + cell->pir[1][1];
    cell->U[1]->momentum.x = cell->U[0]->momentum.x + cell->pir[2][1];
    cell->U[1]->momentum.y = cell->U[0]->momentum.y + cell->pir[3][1];
    cell->U[1]->momentum.z = cell->U[0]->momentum.z + cell->pir[4][1];
    cell->U[1]->total_energy = cell->U[0]->total_energy + cell->pir[5][1];
	
    /* Update single species mass fraction */
    cell->U[1]->massf[0] = cell->U[0]->massf[0] + cell->pir[1][1];
#endif
    return SUCCESS;
} //int point_implicit_update_for_cell

int calculate_M(FV_Cell *cell, int dimensions) 
{
#if WITH_IMPLICIT == 1
    int aa, bb;
    FV_Interface *IFn = cell->iface[NORTH];
    FV_Interface *IFe = cell->iface[EAST];
    FV_Interface *IFs = cell->iface[SOUTH];
    FV_Interface *IFw = cell->iface[WEST];
    FV_Interface *IFt = cell->iface[TOP];
    FV_Interface *IFb = cell->iface[BOTTOM];
	
    calculate_viscous_jacobian(cell, IFn);
    calculate_viscous_jacobian(cell, IFe);
    calculate_viscous_jacobian(cell, IFs);
    calculate_viscous_jacobian(cell, IFw);
    if (dimensions == 3) {
	calculate_viscous_jacobian(cell, IFt);
	calculate_viscous_jacobian(cell, IFb);
    }

    /* 	Initialize matrix piM to zero */
    for (aa = 1; aa <=5; ++aa ) {
	for (bb = 1; bb <=5; ++bb ) {
	    cell->piM[aa][bb] = 0.0;
	}
    }
	
    // Calculate M_viscous as per equation 4.18
    for ( aa = 1; aa <= 5; ++aa ) {
	for ( bb = 1; bb <= 5; ++bb ) {
	    cell->piM[aa][bb] = IFn->J[aa][bb]*IFn->area + IFs->J[aa][bb]*IFs->area +
		IFe->J[aa][bb]*IFe->area + IFw->J[aa][bb]*IFw->area;
	    if (dimensions == 3) 
		cell->piM[aa][bb] += IFt->J[aa][bb]*IFt->area + IFb->J[aa][bb]*IFb->area;
	}	
    }
#endif
    return SUCCESS;
} //int calculate_M

int calculate_viscous_jacobian(FV_Cell *cell, FV_Interface *iface)
{
#if WITH_IMPLICIT == 1
    Gas_model *gmodel = get_gas_model_ptr();
    int aa, bb, statusf;
    double rho, gamma, mu_eff, mu_lam, mu_t, Pr, e_int;
    Vector3 n;
    double Iu,Iv,Iw;
    double Cu,Cv,Cw;	
    double UL;
    double viscous_factor; // so that we can scale down the viscous effects
    double factor;
	
    Iu = iface->fs->vel.x;
    Iv = iface->fs->vel.y;
    Iw = iface->fs->vel.z;
    Cu = cell->fs->vel.x;
    Cv = cell->fs->vel.y;
    Cw = cell->fs->vel.z;
    n = iface->n;
    UL = Cu*n.x + Cv*n.y + Cw*n.z;
    gamma = gmodel->gamma(*(iface->fs->gas), statusf);
    rho = cell->fs->gas->rho;
    Pr =  iface->fs->gas->mu * gmodel->Cp(*(iface->fs->gas), statusf) / iface->fs->gas->k[0];

    viscous_factor = get_viscous_factor();

    mu_lam = viscous_factor * iface->fs->gas->mu;
    mu_t = viscous_factor * iface->fs->mu_t;
    mu_eff = mu_lam + mu_t;
    e_int = cell->U[1]->total_energy - dot(cell->fs->vel, cell->fs->vel);
    	
    /* 	Initalise jacobian matrix J */
    
    for ( aa = 1; aa <=5; ++aa ) {
	for ( bb = 1; bb <=5; ++bb ) {
	    iface->J[aa][bb] = 0.0;
	}
    }
    
    /* 	Fill jacobian matrix as of Appendix C */
    
    iface->J[1][1] = 0.0;
    iface->J[2][1] = -Cu - ( UL * n.x / 3 );
    iface->J[3][1] = -Cv - ( UL * n.y / 3 );
    iface->J[4][1] = -Cw - ( UL * n.z / 3 );
    iface->J[5][1] = Iu * iface->J[2][1] + Iv * iface->J[3][1] + Iw * iface->J[4][1] - gamma / Pr * e_int;
    
    iface->J[1][2] = 0.0;
    iface->J[2][2] = 1 + (n.x*n.x / 3);
    iface->J[3][2] = n.x * n.y / 3;
    iface->J[4][2] = n.x * n.z / 3;
    iface->J[5][2] = Iu * iface->J[2][2] + Iv * iface->J[3][2] + Iw * iface->J[4][2] - gamma * Cu / Pr;
    
    iface->J[1][3] = 0.0;
    iface->J[2][3] = iface->J[3][2];
    iface->J[3][3] =  1 + (n.y*n.y / 3 );
    iface->J[4][3] = n.y * n.z / 3;
    iface->J[5][3] = Iu * iface->J[2][3] + Iv * iface->J[3][3] + Iw * iface->J[4][3] - gamma * Cv / Pr;
    
    iface->J[1][4] = 0.0;
    iface->J[2][4] = iface->J[4][2]; 
    iface->J[3][4] = iface->J[4][3];
    iface->J[4][4] = 1 + (n.z*n.z / 3 );
    iface->J[5][4] = Iu * iface->J[2][4] + Iv * iface->J[3][4] + Iw * iface->J[4][4] - gamma * Cw / Pr;
    
    iface->J[1][5] = 0.0;
    iface->J[2][5] = 0.0;
    iface->J[3][5] = 0.0;
    iface->J[4][5] = 0.0;
    iface->J[5][5] = gamma / Pr;
	
    factor = mu_eff * ( iface->area / cell->volume ) / rho;

    for ( aa = 1; aa <= 5; ++aa ) {
	for ( bb = 1; bb <= 5; ++bb ) {
	    iface->J[aa][bb] = factor * iface->J[aa][bb];
	}
    }	
#endif
    return SUCCESS;
} //int viscous_jacobian

int calculate_h(FV_Cell *cell, int dimensions) 
{
#if WITH_IMPLICIT == 1
    FV_Interface *IFn = cell->iface[NORTH];
    FV_Interface *IFe = cell->iface[EAST];
    FV_Interface *IFs = cell->iface[SOUTH];
    FV_Interface *IFw = cell->iface[WEST];
    FV_Interface *IFt = cell->iface[TOP];
    FV_Interface *IFb = cell->iface[BOTTOM];
	
    /* initialize r to zero */
    for (int aa = 1; aa <= 5; ++aa ) {
	cell->pir[aa][1] = 0.0;
    }
    
    /* Calculates the components of h as per equation 4.11 */
	
    cell->pir[1][1] = 0.0; 
    cell->pir[2][1] = IFe->F->momentum.x * IFe->area - IFw->F->momentum.x * IFw->area +
	IFn->F->momentum.x * IFn->area - IFs->F->momentum.x * IFs->area;
    cell->pir[3][1] = IFe->F->momentum.y * IFe->area - IFw->F->momentum.y * IFw->area +
	IFn->F->momentum.y * IFn->area - IFs->F->momentum.y * IFs->area;
    cell->pir[4][1] = 0.0;
    cell->pir[5][1] = IFe->F->total_energy * IFe->area - IFw->F->total_energy * IFw->area +
	IFn->F->total_energy * IFn->area - IFs->F->total_energy * IFs->area;
	
    if (dimensions == 3){
	cell->pir[2][1] += IFt->F->momentum.x * IFt->area - IFb->F->momentum.x * IFb->area;
	cell->pir[3][1] += IFt->F->momentum.y * IFt->area - IFb->F->momentum.y * IFb->area;
	cell->pir[4][1] += IFe->F->momentum.z * IFe->area - IFw->F->momentum.z * IFw->area +
	    IFn->F->momentum.z * IFn->area - IFs->F->momentum.z * IFs->area +
	    IFt->F->momentum.z * IFt->area - IFb->F->momentum.z * IFb->area;
	cell->pir[5][1] += IFt->F->total_energy * IFt->area - IFb->F->total_energy * IFb->area;
    }
#endif
    return SUCCESS;
} //int calculate_h

void gaussj(FV_Cell *cell, int n, int m)
/* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
 * is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
 * output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
 * vectors. 
 */
{    
    /*                    i  
     *                  ---->
     *         [                     ]
     *         [                     ]  
     *      |  [                     ]
     *   j  |  [                     ] 
     *      v  [                     ]
     *         [                     ]
     *         [                     ]
     *
     */
#if WITH_IMPLICIT == 1
    int pivot_index_i, pivot_index_j;
    int i,j, row, col;
    float dum,pivinv;
    for (j = 1; j <= 5; j++) {
	for (i = 1; i <= 5; i++ ) {
	    cell->gjmtx[i][j] = cell->piM[i][j];
	}
    }
    for (j = 1; j <=5; j++) {
	cell->gjvec[j][1] = cell->pir[j][1];
    }
    for (j=1;j<=n;j++) { /* This is the main loop over the rows.*/
	pivot_index_i = j;
	pivot_index_j = j;
	/* We are now ready to divide the pivot row by the */
	/* pivot element, located at irow and icol.*/
	if (cell->gjmtx[pivot_index_i][pivot_index_j] == 0.0) {
	    fprintf(stderr,"Numerical run-time error...\n");
	    fprintf(stderr,"gaussj: Singular Matrix\n");
	    fprintf(stderr,"...now exiting to system...\n");
	    exit(1);
	}
	pivinv=1.0/cell->gjmtx[pivot_index_i][pivot_index_j];
	for (i=1;i<=n;i++) {
	    cell->gjmtx[pivot_index_j][i] *= pivinv;
	}
	for (i=1;i<=5;i++) {
	    cell->gjvec[pivot_index_j][i] *= pivinv;
	}
	for (row=1;row<=n;row++) {/*  Next, we reduce the rows...*/
	    if (row != pivot_index_j) { /* ...except for the pivot one, of course.*/
		dum=cell->gjmtx[row][pivot_index_i];
		for (col=1;col<=n;col++) {
		    cell->gjmtx[row][col] -= cell->gjmtx[pivot_index_j][col]*dum;
		}
		for (col=1;col<=1;col++) {
		    cell->gjvec[row][col] -= cell->gjvec[pivot_index_j][col]*dum;
		}
	    }
	}
    }
    for (j = 1; j <=5; j++) {
	cell->pir[j][1] = cell->gjvec[j][1];
    }
#endif
} //int gaussj
