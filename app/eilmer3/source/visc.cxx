/** \file visc.cxx
 * \ingroup eilmer3
 *
 * \brief 2D VISCOUS FLUX ROUTINES for the flow solver cns4u.c.
 *
 * 18-nov-92 : the viscous_derivatives function has been split into
 *             7 parts to get around the limitations of the TopSpeed
 *             C compiler.
 * 11-Feb-96 : Moved viscous_coefficients() from cns_misc.c to here.
 * 17-Nov-96 : deleted viscous_coefficients() as the work of computing
 *             the fluid properties is now done (via a call to EOS())
 *             while computing the inviscid fluxes.
 * 17-Jun-97 : Move the fudge for estimating the flow derivatives at
 *             the corner vertices to a separate function that is called
 *             after each of the edges is handled.
 * 16-Oct-97 : Clean up vectorisation for Cray C90
 * 23-Jul-2006 : last remains of vectorization removed.
 * 02-Mar-2008: Eilmer3 port
 * 16-Nov-2012: Added function to check velocity direction and get
 *              derivative from upwind direction.
 */

/*-----------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <stdexcept>
#include "block.hh"
#include "kernel.hh"
#include "bc.hh"
#include "visc.hh"
#include "diffusion.hh"
#include "baldwin_lomax.hh"

constexpr double VERY_SMALL = 1.0e-10;

// Working arrays for species derivatives
static std::vector<double> dfdx, dfdy, dfdz, jx, jy, jz;
// Working arrays for thermal derivatives
static std::vector<double> dTdx, dTdy, qx, qy, k_eff, TA, TB, TC, TD;
// Although we don't actually use dfdz and jz, they are needed
// as place holders because the calculation of diffusion fluxes
// treats a general 3D problem.
static std::vector<double> fA, fB, fC, fD;


/*=================================================================*/

/// \brief Decide which turbulence model to apply.
/// 
/// The turbulence interacts with the average-flow field
/// primarily through the diffusive terms.
int estimate_turbulence_viscosity(global_data *gdp, Block *bdp)
{
    switch ( gdp->turbulence_model ) {
    case TM_NONE:
	for ( FV_Cell *cp: bdp->active_cells ) cp->turbulence_viscosity_zero();
	return SUCCESS;
    case TM_BALDWIN_LOMAX:
	baldwin_lomax_turbulence_model(*gdp, *bdp, 0);
	break;
    case TM_K_OMEGA:
	for ( FV_Cell *cp: bdp->active_cells ) cp->turbulence_viscosity_k_omega();
	break;
    default:
	throw std::runtime_error("Turbulence model requested but not available.");
    }
    for ( FV_Cell *cp: bdp->active_cells ) {
	cp->turbulence_viscosity_factor(gdp->transient_mu_t_factor);
	cp->turbulence_viscosity_limit(gdp->max_mu_t_factor);
	cp->turbulence_viscosity_zero_if_not_in_zone();
    }
    return SUCCESS;
}

/// \brief Compute the viscous contribution to the cell interface fluxes.
///
/// This contribution is added to the flux variables so make sure that
/// the inviscid (Euler) contributions have already been computed and stored.
///
int viscous_flux_2D(Block *A)
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    FV_Vertex *vtx1, *vtx2;
    FV_Interface *IFace;
    double nx, ny;
    double dudx, dudy, dvdx, dvdy;
    double mu_lam; // molecular-scale transport properties
    double lmbda; // second-coefficient of viscosity
    double mu_t;  // turbulence viscosity
    double k_t;   // turbulence thermal conductivity
    double viscous_factor; // so that we can scale down the viscous effects
    double diffusion_factor; // so that we can scale down the diffusion effects
    double mu_eff; // combined laminar and turbulent viscosity after scaling
    double dtkedx, dtkedy, domegadx, domegady;
    double sigma = 0.5;
    double sigma_star = 0.6;
    double mu_effective;
    double tau_kx, tau_ky, tau_wx, tau_wy;
    double tau_xx, tau_yy, tau_xy;
    double ybar;
    double tau_wall_x, tau_wall_y;

    size_t nsp = gmodel->get_number_of_species();
    if( dfdx.size() == 0 ) {
	dfdx.resize(nsp); 
	dfdy.resize(nsp);
	dfdz.resize(nsp);
	jx.resize(nsp);
	jy.resize(nsp);
	jz.resize(nsp);
    }
    
    size_t ntm = gmodel->get_number_of_modes();
    if( dTdx.size() == 0 ) {
	dTdx.resize(ntm);
	dTdy.resize(ntm);
	qx.resize(ntm);
	qy.resize(ntm);
	k_eff.resize(ntm);
    }
    
    viscous_factor = G.viscous_factor;
    diffusion_factor = G.diffusion_factor;

    // East-facing interfaces.
    for ( size_t i = A->imin; i <= A->imax+1; ++i ) {
        for ( size_t j = A->jmin; j <= A->jmax; ++j ) {
	    if ( (i == A->imin && A->bcp[WEST]->sets_visc_flux()) ||
		 (i == A->imax+1 && A->bcp[EAST]->sets_visc_flux()) ) {
		// Retain the b.c. set flux by doing nothing, just continue
		continue;
	    }
            IFace = A->get_ifi(i,j);
	    FlowState &fs = *(IFace->fs);
            vtx1 = A->get_vtx(i,j+1);
            vtx2 = A->get_vtx(i,j);
	    // Determine some of the interface properties.
            ybar = IFace->Ybar;
            dudx = 0.5*(vtx1->dudx + vtx2->dudx);
            dudy = 0.5*(vtx1->dudy + vtx2->dudy);
            dvdx = 0.5*(vtx1->dvdx + vtx2->dvdx);
            dvdy = 0.5*(vtx1->dvdy + vtx2->dvdy);
            for ( size_t itm=0; itm<ntm; ++itm ) {
                dTdx[itm] = 0.5*(vtx1->dTdx[itm] + vtx2->dTdx[itm]);
                dTdy[itm] = 0.5*(vtx1->dTdy[itm] + vtx2->dTdy[itm]);
                k_eff[itm] = viscous_factor * fs.gas->k[itm];
            }
            mu_lam = viscous_factor * fs.gas->mu;
            mu_t = viscous_factor * fs.mu_t;
            k_t = viscous_factor * fs.k_t;
            // CHECK-ME: Assume turbulent conductivity only applies to primary mode?
            k_eff[0] += k_t;
            mu_eff = mu_lam + mu_t;
            lmbda = -2.0/3.0 * mu_eff;
            if ( G.diffusion ) {
                for ( size_t isp = 0; isp < nsp; ++isp ) {
                    dfdx[isp] = 0.5*(vtx1->dfdx[isp] + vtx2->dfdx[isp]);
                    dfdy[isp] = 0.5*(vtx1->dfdy[isp] + vtx2->dfdy[isp]);
                }
	        // Apply a diffusion model
	        double D_t = 0.0;
	        if ( G.turbulence_model != TM_NONE ) {
                    double Sc_t = G.turbulence_schmidt;
                    D_t = mu_t / (fs.gas->rho * Sc_t);
	        }
	        calculate_diffusion_fluxes(*(fs.gas),
				           D_t,
					   dfdx, dfdy, dfdz,
					   jx, jy, jz);
		// NOTE: now applying diffusion_factor to diffusive fluxes --> instead of viscous_factor!
		for ( size_t isp = 0; isp < nsp; ++isp ) {
		    jx[isp] *= diffusion_factor;
		    jy[isp] *= diffusion_factor;
		    jz[isp] *= diffusion_factor;
		}
	    }
	    
	    if ( G.axisymmetric ) {
		// Viscous stresses at the mid-point of the interface.
		// Axisymmetric terms no longer include the radial multiplier
		// as that has been absorbed into the interface area calculation.
                if (ybar > VERY_SMALL) {
                    tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy + fs.vel.y / ybar);
                    tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy + fs.vel.y / ybar);
                } else {
                    tau_xx = 0.0;
                    tau_yy = 0.0;
                }
                tau_xy = mu_eff * (dudy + dvdx);
	    } else {
		// 2-dimensional-planar stresses.
                tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy);
                tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy);
                tau_xy = mu_eff * (dudy + dvdx);
	    }  
	    
	    // Thermal conductivity
	    // NOTE: q[0] is total energy flux
	    qx[0] = k_eff[0] * dTdx[0];
	    qy[0] = k_eff[0] * dTdy[0];
	    for ( size_t itm=1; itm<ntm; ++itm ) {
		qx[itm] = k_eff[itm] * dTdx[itm];
		qy[itm] = k_eff[itm] * dTdy[itm];
		qx[0] += qx[itm];
		qy[0] += qy[itm];
	    }
	    
	    if ( G.diffusion ) {
		for ( size_t isp = 0; isp < nsp; ++isp ) {
		    double h = gmodel->enthalpy(*(fs.gas), isp);
		    qx[0] -= jx[isp] * h;
		    qy[0] -= jy[isp] * h;
		    for ( size_t itm=1; itm<ntm; ++itm ) {
			double hmode = gmodel->modal_enthalpy(*(fs.gas), isp, itm);
			qx[itm] -= jx[isp] * hmode;
			qy[itm] -= jy[isp] * hmode;
		    }
		}
	    }	    
	    
	    if ( G.turbulence_model == TM_K_OMEGA &&
		 !(G.axisymmetric && (ybar <= VERY_SMALL)) ) {
		// Turbulence contribution to the shear stresses.
		// 2016-05-30 PJ.  If we are doing an axisymmetric flow,
		// only have non-zero shear stresses if we are not near the axis.
		tau_xx -= 0.66667 * fs.gas->rho * fs.tke;
		tau_yy -= 0.66667 * fs.gas->rho * fs.tke;
		// Turbulence contribution to heat transfer.
		mu_effective = mu_lam + sigma_star * mu_t;
		dtkedx = 0.5*(vtx1->dtkedx + vtx2->dtkedx);
		dtkedy = 0.5*(vtx1->dtkedy + vtx2->dtkedy);
		qx[0] += mu_effective * dtkedx;
		qy[0] += mu_effective * dtkedy;
		// Turbulence transport of the turbulence properties themselves.
		tau_kx = mu_effective * dtkedx; 
		tau_ky = mu_effective * dtkedy;
		mu_effective = mu_lam + sigma * mu_t;
		domegadx = 0.5*(vtx1->domegadx + vtx2->domegadx);
		domegady = 0.5*(vtx1->domegady + vtx2->domegady);
		tau_wx = mu_effective * domegadx; 
		tau_wy = mu_effective * domegady; 
	    } else {
		tau_kx = 0.0;
		tau_ky = 0.0;
		tau_wx = 0.0;
		tau_wy = 0.0;
	    }

	    // Combine into fluxes: store as the dot product (F.n).
	    ConservedQuantities &F = *(IFace->F);
            nx = IFace->n.x;
            ny = IFace->n.y;
            // Mass flux -- NO CONTRIBUTION
	    if ( (G.wall_function && i == A->imin && A->bcp[WEST]->is_wall() && A->bcp[WEST]->type_code != SLIP_WALL)
	    || (G.wall_function && i == A->imax+1 && A->bcp[EAST]->is_wall() && A->bcp[EAST]->type_code != SLIP_WALL) ) {
	        tau_wall_x = IFace->tau_wall_x;
	        tau_wall_y = IFace->tau_wall_y;	    	  	        
                F.momentum.x -= tau_xx * nx + tau_wall_x;
                F.momentum.y -= tau_yy * ny + tau_wall_y;       
		F.total_energy -= tau_xx * fs.vel.x * nx + tau_yy * fs.vel.y * ny + 
		                  tau_wall_x * fs.vel.x + tau_wall_y * fs.vel.y +
		                  IFace->q_wall;
	    } else {
                F.momentum.x -= tau_xx * nx + tau_xy * ny;
                F.momentum.y -= tau_xy * nx + tau_yy * ny;
                F.total_energy -= (tau_xx * fs.vel.x + tau_xy * fs.vel.y + qx[0]) * nx
                    + (tau_xy * fs.vel.x + tau_yy * fs.vel.y + qy[0]) * ny;	    
	    }
	    // Viscous transport of k-omega turbulence quantities.
	    // Only built for 2D planar geometry at the moment.
	    if ( G.turbulence_model == TM_K_OMEGA ) {
		F.tke -= tau_kx * nx + tau_ky * ny;
		F.omega -= tau_wx * nx + tau_wy * ny;
	    }
            // Species mass flux
	    if( G.diffusion ) {
		if ( ( i == A->imin && A->bcp[WEST]->type_code == USER_DEFINED_MASS_FLUX ) ||
		     ( i == A->imax+1 && A->bcp[EAST]->type_code == USER_DEFINED_MASS_FLUX ) ) {
		      // Retain species mass flux set earlier
		    ; // Do nothing statement.
		}
		else {
		    for ( size_t isp = 0; isp < nsp; ++isp ) {
			F.massf[isp] += jx[isp] * nx + jy[isp] * ny;
		    }
		}
	    }
	    // Modal energy flux (skipping first mode as this is handled by total energy)
	    for ( size_t itm=1; itm<ntm; ++itm ) {
	    	F.energies[itm] -= qx[itm] * nx + qy[itm] * ny;
	    }
        } // j loop
    } // i loop
    /*
     * North-facing interfaces
     */
    for ( size_t i = A->imin; i <= A->imax; ++i ) {
        for ( size_t j = A->jmin; j <= A->jmax+1; ++j ) {
	    if ( (j == A->jmin && A->bcp[SOUTH]->sets_visc_flux()) ||
		 (j == A->jmax+1 && A->bcp[NORTH]->sets_visc_flux()) ) {
		// Retain the b.c. set flux by doing nothing, just continue
		continue;
	    }
            IFace = A->get_ifj(i,j);
	    FlowState &fs = *(IFace->fs);
            vtx1 = A->get_vtx(i,j);
            vtx2 = A->get_vtx(i+1,j);
	    // Determine some of the interface properties.
            ybar = IFace->Ybar;
            dudx = 0.5*(vtx1->dudx + vtx2->dudx);
            dudy = 0.5*(vtx1->dudy + vtx2->dudy);
            dvdx = 0.5*(vtx1->dvdx + vtx2->dvdx);
            dvdy = 0.5*(vtx1->dvdy + vtx2->dvdy);
            for ( size_t itm=0; itm<ntm; ++itm ) {
                dTdx[itm] = 0.5*(vtx1->dTdx[itm] + vtx2->dTdx[itm]);
                dTdy[itm] = 0.5*(vtx1->dTdy[itm] + vtx2->dTdy[itm]);
                k_eff[itm] = viscous_factor * fs.gas->k[itm];
            }
            mu_lam = viscous_factor * fs.gas->mu;
	    mu_t = viscous_factor * fs.mu_t;
	    k_t = viscous_factor * fs.k_t;
	    // CHECK-ME: Assume turbulent conductivity only applies to primary mode?
	    k_eff[0] += k_t;
	    mu_eff = mu_lam + mu_t;
	    lmbda = -2.0/3.0 * mu_eff;
            if ( G.diffusion ) {
                for ( size_t isp = 0; isp < nsp; ++isp ) {
                    dfdx[isp] = 0.5*(vtx1->dfdx[isp] + vtx2->dfdx[isp]);
                    dfdy[isp] = 0.5*(vtx1->dfdy[isp] + vtx2->dfdy[isp]);
                }
		// Apply a diffusion model
		double D_t = 0.0;
		if ( G.turbulence_model != TM_NONE ) {
                    double Sc_t = G.turbulence_schmidt;
                    D_t = mu_t / (fs.gas->rho * Sc_t);
		}
		calculate_diffusion_fluxes(*(fs.gas),
					   D_t,
					   dfdx, dfdy, dfdz,
					   jx, jy, jz);
		// NOTE: now applying diffusion_factor to diffusive fluxes --> instead of viscous_factor!
		for( size_t isp = 0; isp < nsp; ++isp ) {
		    jx[isp] *= diffusion_factor;
		    jy[isp] *= diffusion_factor;
		    jz[isp] *= diffusion_factor;
		}
	    }

	    if ( G.axisymmetric ) {
		// Viscous stresses at the mid-point of the interface.
		// Axisymmetric terms no longer include the radial multiplier
		// as that has been absorbed into the interface area calculation.
                if (ybar > VERY_SMALL) {
                    tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy + fs.vel.y / ybar);
                    tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy + fs.vel.y / ybar);
                } else {
                    tau_xx = 0.0;
                    tau_yy = 0.0;
                }
                tau_xy = mu_eff * (dudy + dvdx);
	    } else {
		// 2-dimensional-planar stresses.
                tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy);
                tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy);
                tau_xy = mu_eff * (dudy + dvdx);
	    }

	    // Thermal conductivity
	    // NOTE: q[0] is total energy flux
	    // [todo] 2013-12-04 check the comment above.
	    qx[0] = k_eff[0] * dTdx[0];
	    qy[0] = k_eff[0] * dTdy[0];
	    for ( size_t itm=1; itm<ntm; ++itm ) {
		qx[itm] = k_eff[itm] * dTdx[itm];
		qy[itm] = k_eff[itm] * dTdy[itm];
		qx[0] += qx[itm];
		qy[0] += qy[itm];
	    }
	    
	    if( G.diffusion ) {
		for( size_t isp = 0; isp < nsp; ++isp ) {
		    double h = gmodel->enthalpy(*(fs.gas), isp);
		    qx[0] -= jx[isp] * h;
		    qy[0] -= jy[isp] * h;
		    for ( size_t itm=1; itm<ntm; ++itm ) {
			double hmode = gmodel->modal_enthalpy(*(fs.gas), isp, itm);
			qx[itm] -= jx[isp] * hmode;
			qy[itm] -= jy[isp] * hmode;
		    }
		}
	    }
	    
	    if ( G.turbulence_model == TM_K_OMEGA &&
		 !(G.axisymmetric && (ybar <= VERY_SMALL)) ) {
		// Turbulence contribution to the shear stresses.
		// 2016-05-30 PJ.  If we are doing an axisymmetric flow,
		// only have non-zero shear stresses if we are not near the axis.
		tau_xx -= 0.66667 * fs.gas->rho * fs.tke;
		tau_yy -= 0.66667 * fs.gas->rho * fs.tke;
		// Turbulence contribution to heat transfer.
		mu_effective = mu_lam + sigma_star * mu_t;
		dtkedx = 0.5*(vtx1->dtkedx + vtx2->dtkedx);
		dtkedy = 0.5*(vtx1->dtkedy + vtx2->dtkedy);
		qx[0] += mu_effective * dtkedx;
		qy[0] += mu_effective * dtkedy;
		// Turbulence transport of the turbulence properties themselves.
		tau_kx = mu_effective * dtkedx; 
		tau_ky = mu_effective * dtkedy;
		mu_effective = mu_lam + sigma * mu_t;
		domegadx = 0.5*(vtx1->domegadx + vtx2->domegadx);
		domegady = 0.5*(vtx1->domegady + vtx2->domegady);
		tau_wx = mu_effective * domegadx; 
		tau_wy = mu_effective * domegady;
	    } else {
		tau_kx = 0.0;
		tau_ky = 0.0;
		tau_wx = 0.0;
		tau_wy = 0.0;
	    }

	    // Combine into fluxes: store as the dot product (F.n).
	    ConservedQuantities &F = *(IFace->F);
            nx = IFace->n.x;
            ny = IFace->n.y;
	    // Mass flux -- NO CONTRIBUTION            
	    if ( (G.wall_function && j == A->jmin && A->bcp[SOUTH]->is_wall() && A->bcp[SOUTH]->type_code != SLIP_WALL)
	    || (G.wall_function && j == A->jmax+1 && A->bcp[NORTH]->is_wall() && A->bcp[NORTH]->type_code != SLIP_WALL) ) {   	            
	        tau_wall_x = IFace->tau_wall_x;
	        tau_wall_y = IFace->tau_wall_y;	        
                F.momentum.x -= tau_xx * nx + tau_wall_x;
                F.momentum.y -= tau_yy * ny + tau_wall_y;	    
		F.total_energy -= tau_xx * fs.vel.x * nx + tau_yy * fs.vel.y * ny + 
		                  tau_wall_x * fs.vel.x + tau_wall_y * fs.vel.y +
		                  IFace->q_wall;
	    } else {
                F.momentum.x -= tau_xx * nx + tau_xy * ny;
                F.momentum.y -= tau_xy * nx + tau_yy * ny;
	        if ( j == A->jmax+1 && ( A->bcp[NORTH]->type_code == USER_DEFINED_ENERGY_FLUX ||
				       ( (A->bcp[NORTH]->type_code == CONJUGATE_HT) &&
				       ( (G.cht_coupling == TFS_QWS) || (G.cht_coupling == QFS_QWS) ) ) ) ) {
		    // Retain the flux set by the b.c. by doing nothing.
		    ; // Do nothing statement
	        }
	        else {
		    F.total_energy -= (tau_xx * fs.vel.x + tau_xy * fs.vel.y + qx[0]) * nx
		        + (tau_xy * fs.vel.x + tau_yy * fs.vel.y + qy[0]) * ny;	    
	        }
	    } 
	    	    
	    // Viscous transport of k-omega turbulence quantities.
	    // Only built for 2D planar geometry at the moment.
	    if ( G.turbulence_model == TM_K_OMEGA ) {
		F.tke -= tau_kx * nx + tau_ky * ny;
		F.omega -= tau_wx * nx + tau_wy * ny;
	    }
	    // Species mass flux
	    if( G.diffusion ) {
		if ( (j == A->jmin && A->bcp[SOUTH]->type_code == USER_DEFINED_MASS_FLUX ) ||
		     (j == A->jmax+1 && A->bcp[NORTH]->type_code == USER_DEFINED_MASS_FLUX ) ) {
		    // Retain the b.c. set species fluxes by doing nothing, just continue
		    ; // Do nothing statement
		}
		else {
		    for( size_t isp = 0; isp < nsp; ++isp ) {
			F.massf[isp] += jx[isp] * nx + jy[isp] * ny;
		    }
		}
	    }
	    // Modal energy flux (skipping first mode as this is handled by total energy)
	    for ( size_t itm=1; itm<ntm; ++itm ) {
	    	F.energies[itm] -= qx[itm] * nx + qy[itm] * ny;
	    }
        } // j loop
    } // i loop
    return SUCCESS;
} // end viscous_flux()


// --------------------
// Viscous Book keeping
// --------------------
#define APPLY_DIVERGENCE_THEOREM() \
            A->get_vtx(i,j)->dudx = 0.5 * area_inv * \
		((uB + uA) * (yB - yA) + (uC + uB) * (yC - yB) + \
		 (uD + uC) * (yD - yC) + (uA + uD) * (yA - yD)); \
            A->get_vtx(i,j)->dudy = -0.5 * area_inv * \
		((uB + uA) * (xB - xA) + (uC + uB) * (xC - xB) + \
		 (uD + uC) * (xD - xC) + (uA + uD) * (xA - xD)); \
            A->get_vtx(i,j)->dvdx = 0.5 * area_inv * \
		((vB + vA) * (yB - yA) + (vC + vB) * (yC - yB) + \
		 (vD + vC) * (yD - yC) + (vA + vD) * (yA - yD)); \
            A->get_vtx(i,j)->dvdy = -0.5 * area_inv * \
		((vB + vA) * (xB - xA) + (vC + vB) * (xC - xB) + \
		 (vD + vC) * (xD - xC) + (vA + vD) * (xA - xD)); \
	    for ( size_t itm=0; itm<ntm; ++itm ) { \
	        A->get_vtx(i,j)->dTdx[itm] = 0.5 * area_inv * \
		    ((TB[itm] + TA[itm]) * (yB - yA) + (TC[itm] + TB[itm]) * (yC - yB) + \
		     (TD[itm] + TC[itm]) * (yD - yC) + (TA[itm] + TD[itm]) * (yA - yD)); \
                A->get_vtx(i,j)->dTdy[itm] = -0.5 * area_inv * \
		    ((TB[itm] + TA[itm]) * (xB - xA) + (TC[itm] + TB[itm]) * (xC - xB) + \
		     (TD[itm] + TC[itm]) * (xD - xC) + (TA[itm] + TD[itm]) * (xA - xD)); \
            } \
            if( G.diffusion ) { \
                for( size_t isp = 0; isp < nsp; ++isp ) { \
                    A->get_vtx(i,j)->dfdx[isp] = 0.5 * area_inv * \
		        ((fB[isp] + fA[isp]) * (yB - yA) + (fC[isp] + fB[isp]) * (yC - yB) + \
		        (fD[isp] + fC[isp]) * (yD - yC) + (fA[isp] + fD[isp]) * (yA - yD)); \
                    A->get_vtx(i,j)->dfdy[isp] = -0.5 * area_inv * \
		        ((fB[isp] + fA[isp]) * (xB - xA) + (fC[isp] + fB[isp]) * (xC - xB) + \
		        (fD[isp] + fC[isp]) * (xD - xC) + (fA[isp] + fD[isp]) * (xA - xD)); \
	        } \
	    } \
	    if( G.electric_field_work ) { \
	        A->get_vtx(i,j)->dpedx = -0.5 * area_inv * \
		    ((peB + peA) * (yB - yA) + (peC + peB) * (yC - yB) + \
		    (peD + peC) * (yD - yC) + (peA + peD) * (yA - yD)); \
	        A->get_vtx(i,j)->dpedy = -0.5 * area_inv * \
		    ((peB + peA) * (xB - yA) + (peC + peB) * (xC - xB) + \
		    (peD + peC) * (xD - xC) + (peA + peD) * (xA - xD)); \
	    }

#define APPLY_DIVERGENCE_THEOREM_2() \
            A->get_vtx(i,j)->dtkedx = 0.5 * area_inv * \
		((tkeB + tkeA) * (yB - yA) + (tkeC + tkeB) * (yC - yB) + \
		 (tkeD + tkeC) * (yD - yC) + (tkeA + tkeD) * (yA - yD)); \
            A->get_vtx(i,j)->dtkedy = -0.5 * area_inv * \
		((tkeB + tkeA) * (xB - xA) + (tkeC + tkeB) * (xC - xB) + \
		 (tkeD + tkeC) * (xD - xC) + (tkeA + tkeD) * (xA - xD)); \
            A->get_vtx(i,j)->domegadx = 0.5 * area_inv * \
		((omegaB + omegaA) * (yB - yA) + (omegaC + omegaB) * (yC - yB) + \
		 (omegaD + omegaC) * (yD - yC) + (omegaA + omegaD) * (yA - yD)); \
            A->get_vtx(i,j)->domegady = -0.5 * area_inv * \
		((omegaB + omegaA) * (xB - xA) + (omegaC + omegaB) * (xC - xB) + \
		 (omegaD + omegaC) * (xD - xC) + (omegaA + omegaD) * (xA - xD));

#define APPLY_DIVERGENCE_THEOREM_3() \
            A->get_vtx(i,j)->dsigma_Tdx = 0.5 * area_inv *				\
		((sigma_TB + sigma_TA) * (yB - yA) + (sigma_TC + sigma_TB) * (yC - yB) + \
		 (sigma_TD + sigma_TC) * (yD - yC) + (sigma_TA + sigma_TD) * (yA - yD)); \
            A->get_vtx(i,j)->dsigma_Tdy = -0.5 * area_inv * \
		((sigma_TB + sigma_TA) * (xB - xA) + (sigma_TC + sigma_TB) * (xC - xB) + \
		 (sigma_TD + sigma_TC) * (xD - xC) + (sigma_TA + sigma_TD) * (xA - xD)); \
            A->get_vtx(i,j)->dsigma_cdx = 0.5 * area_inv * \
		((sigma_cB + sigma_cA) * (yB - yA) + (sigma_cC + sigma_cB) * (yC - yB) + \
		 (sigma_cD + sigma_cC) * (yD - yC) + (sigma_cA + sigma_cD) * (yA - yD)); \
            A->get_vtx(i,j)->dsigma_cdy = -0.5 * area_inv * \
		((sigma_cB + sigma_cA) * (xB - xA) + (sigma_cC + sigma_cB) * (xC - xB) + \
		 (sigma_cD + sigma_cC) * (xD - xC) + (sigma_cA + sigma_cD) * (xA - xD)); \
           for( size_t isp = 0; isp < nsp; ++isp ) { \
                    A->get_vtx(i,j)->dcdx[isp] = 0.5 * area_inv * \
		        ((cB[isp] + cA[isp]) * (yB - yA) + (cC[isp] + cB[isp]) * (yC - yB) + \
		        (cD[isp] + cC[isp]) * (yD - yC) + (cA[isp] + cD[isp]) * (yA - yD)); \
                    A->get_vtx(i,j)->dcdy[isp] = -0.5 * area_inv * \
		        ((cB[isp] + cA[isp]) * (xB - xA) + (cC[isp] + cB[isp]) * (xC - xB) + \
		        (cD[isp] + cC[isp]) * (xD - xC) + (cA[isp] + cD[isp]) * (xA - xD)); \
	   }

/// \brief Compute the derivatives of velocity, temperature and mass fractions 
///        at primary cell vertices.
///
/// This is done by applying the divergence theorem around the perimeters 
/// of SECONDARY cells.  The secondary cells are centred on the vertices 
/// of the primary cells and have primary cell centres as their corners.
///
int viscous_derivatives_2D(Block *A, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double uA, uB, uC, uD;
    double vA, vB, vC, vD;
    double tkeA, tkeB, tkeC, tkeD;
    double omegaA, omegaB, omegaC, omegaD;
    double peA = 0.0, peB = 0.0, peC = 0.0, peD = 0.0;
    double area_inv;
    
    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    if( fA.size() == 0 ) {
	fA.resize(nsp); 
	fB.resize(nsp);
	fC.resize(nsp);
	fD.resize(nsp);
    }
    
    size_t ntm = get_gas_model_ptr()->get_number_of_modes();
    if ( TA.size() == 0 ) {
    	TA.resize(ntm);
    	TB.resize(ntm);
    	TC.resize(ntm);
    	TD.resize(ntm);
    }

    // First, do all of the internal secondary cells.
    // i.e. Those not on a boundary.
    for ( size_t i = A->imin+1; i <= A->imax; ++i ) {
        for ( size_t j = A->jmin+1; j <= A->jmax; ++j ) {
            area_inv = 1.0 / A->get_vtx(i,j)->area;
	    // These are the corners of the secondary cell.
            xA = A->get_cell(i,j-1)->pos[gtl].x;
            yA = A->get_cell(i,j-1)->pos[gtl].y;
            xB = A->get_cell(i,j)->pos[gtl].x;
            yB = A->get_cell(i,j)->pos[gtl].y;
            xC = A->get_cell(i-1,j)->pos[gtl].x;
            yC = A->get_cell(i-1,j)->pos[gtl].y;
            xD = A->get_cell(i-1,j-1)->pos[gtl].x;
            yD = A->get_cell(i-1,j-1)->pos[gtl].y;
	    // These are the flow properties at the corners.
            uA = A->get_cell(i,j-1)->fs->vel.x;
            uB = A->get_cell(i,j)->fs->vel.x;
            uC = A->get_cell(i-1,j)->fs->vel.x;
            uD = A->get_cell(i-1,j-1)->fs->vel.x;
	    //
            vA = A->get_cell(i,j-1)->fs->vel.y;
            vB = A->get_cell(i,j)->fs->vel.y;
            vC = A->get_cell(i-1,j)->fs->vel.y;
            vD = A->get_cell(i-1,j-1)->fs->vel.y;
	    //
	    for ( size_t itm=0; itm<ntm; ++itm ) {
		TA[itm] = A->get_cell(i,j-1)->fs->gas->T[itm];
		TB[itm] = A->get_cell(i,j)->fs->gas->T[itm];
		TC[itm] = A->get_cell(i-1,j)->fs->gas->T[itm];
		TD[itm] = A->get_cell(i-1,j-1)->fs->gas->T[itm];
            }
	    //
	    if( G.diffusion ) {
		for( size_t isp = 0; isp < nsp; ++isp ) {
		    fA[isp] = A->get_cell(i,j-1)->fs->gas->massf[isp];
		    fB[isp] = A->get_cell(i,j)->fs->gas->massf[isp];
		    fC[isp] = A->get_cell(i-1,j)->fs->gas->massf[isp];
		    fD[isp] = A->get_cell(i-1,j-1)->fs->gas->massf[isp];
		}
	    }
	    //
	    if ( G.electric_field_work ) {
	        peA = A->get_cell(i,j-1)->fs->gas->p_e;
	        peB = A->get_cell(i,j)->fs->gas->p_e;
	        peC = A->get_cell(i-1,j)->fs->gas->p_e;
	        peD = A->get_cell(i-1,j-1)->fs->gas->p_e;
	    }
	    //
	    APPLY_DIVERGENCE_THEOREM()
	    //
            tkeA = A->get_cell(i,j-1)->fs->tke;
            tkeB = A->get_cell(i,j)->fs->tke;
            tkeC = A->get_cell(i-1,j)->fs->tke;
            tkeD = A->get_cell(i-1,j-1)->fs->tke;
	    //
            omegaA = A->get_cell(i,j-1)->fs->omega;
            omegaB = A->get_cell(i,j)->fs->omega;
            omegaC = A->get_cell(i-1,j)->fs->omega;
            omegaD = A->get_cell(i-1,j-1)->fs->omega;
	    //
	    APPLY_DIVERGENCE_THEOREM_2()
        } // j loop
    } // i loop

    // Now, do the boundaries as half cells.
    viscous_derivatives_edges(A, gtl);
    // ...and pick up (fudge) the corner values.
    viscous_derivatives_corners(A, gtl);
    return SUCCESS;
} // end viscous_derivatives()


/// \brief Compute the derivatives of velocity, temperature and mass fractions at
///        primary cell vertices along the block boundaries.
///
/// NOTE that the secondary cells along boundaries are HALF cells.
///
int viscous_derivatives_edges(Block *A, size_t gtl)
{
    global_data &G = *get_global_data_ptr();
    size_t i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double uA, uB, uC, uD;
    double vA, vB, vC, vD;
    double tkeA, tkeB, tkeC, tkeD;
    double omegaA, omegaB, omegaC, omegaD;
    double peA = 0.0, peB = 0.0, peC = 0.0, peD = 0.0;
    double area_inv;

    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    if( fA.size() == 0 ) {
	fA.resize(nsp); 
	fB.resize(nsp);
	fC.resize(nsp);
	fD.resize(nsp);
    }
    
    size_t ntm = get_gas_model_ptr()->get_number_of_modes();
    if ( TA.size() == 0 ) {
    	TA.resize(ntm);
    	TB.resize(ntm);
    	TC.resize(ntm);
    	TD.resize(ntm);
    }

    // EAST boundary
    i = A->imax+1;
    for (j = A->jmin+1; j <= A->jmax; ++j) {
        area_inv = 1.0 / A->get_vtx(i,j)->area;
	// These are the corners of the secondary cell.
	xA = A->get_ifi(i,j-1)->pos.x;
	yA = A->get_ifi(i,j-1)->pos.y;
	xB = A->get_ifi(i,j)->pos.x;
	yB = A->get_ifi(i,j)->pos.y;
        xC = A->get_cell(i-1,j)->pos[gtl].x;
        yC = A->get_cell(i-1,j)->pos[gtl].y;
        xD = A->get_cell(i-1,j-1)->pos[gtl].x;
        yD = A->get_cell(i-1,j-1)->pos[gtl].y;
	//
	// These are the flow properties at the corners of the secondary cell.
	uA = A->get_ifi(i,j-1)->fs->vel.x;
	uB = A->get_ifi(i,j)->fs->vel.x;
        uC = A->get_cell(i-1,j)->fs->vel.x;
        uD = A->get_cell(i-1,j-1)->fs->vel.x;
	//
	vA = A->get_ifi(i,j-1)->fs->vel.y;
	vB = A->get_ifi(i,j)->fs->vel.y;
        vC = A->get_cell(i-1,j)->fs->vel.y;
        vD = A->get_cell(i-1,j-1)->fs->vel.y;
	//
	for ( size_t itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_ifi(i,j-1)->fs->gas->T[itm];
	    TB[itm] = A->get_ifi(i,j)->fs->gas->T[itm];
	    TC[itm] = A->get_cell(i-1,j)->fs->gas->T[itm];
	    TD[itm] = A->get_cell(i-1,j-1)->fs->gas->T[itm];
	}
	//
	if( G.diffusion ) {
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_ifi(i,j-1)->fs->gas->massf[isp];
		fB[isp] = A->get_ifi(i,j)->fs->gas->massf[isp];
		fC[isp] = A->get_cell(i-1,j)->fs->gas->massf[isp];
		fD[isp] = A->get_cell(i-1,j-1)->fs->gas->massf[isp];
	    }
	}
	//
	if ( G.electric_field_work ) {
	    peA = A->get_ifi(i,j-1)->fs->gas->p_e;
	    peB = A->get_ifi(i,j)->fs->gas->p_e;
	    peC = A->get_ifi(i-1,j)->fs->gas->p_e;
	    peD = A->get_ifi(i-1,j-1)->fs->gas->p_e;
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
	tkeA = A->get_ifi(i,j-1)->fs->tke;
	tkeB = A->get_ifi(i,j)->fs->tke;
	tkeC = A->get_cell(i-1,j)->fs->tke;
	tkeD = A->get_cell(i-1,j-1)->fs->tke;
	//
	omegaA = A->get_ifi(i,j-1)->fs->omega;
	omegaB = A->get_ifi(i,j)->fs->omega;
	omegaC = A->get_cell(i-1,j)->fs->omega;
	omegaD = A->get_cell(i-1,j-1)->fs->omega;
	//
	APPLY_DIVERGENCE_THEOREM_2()
	//
    } // j loop

    // WEST boundary
    i = A->imin;
    for (j = A->jmin+1; j <= A->jmax; ++j) {
        area_inv = 1.0 / A->get_vtx(i,j)->area;
	// These are the corners of the secondary cell.
        xA = A->get_cell(i,j-1)->pos[gtl].x;
        yA = A->get_cell(i,j-1)->pos[gtl].y;
        xB = A->get_cell(i,j)->pos[gtl].x;
        yB = A->get_cell(i,j)->pos[gtl].y;
	xC = A->get_ifi(i,j)->pos.x;
	yC = A->get_ifi(i,j)->pos.y;
	xD = A->get_ifi(i,j-1)->pos.x;
	yD = A->get_ifi(i,j-1)->pos.y;
	// These are the flow properties at the corners of the secondary cell.
        uA = A->get_cell(i,j-1)->fs->vel.x;
        uB = A->get_cell(i,j)->fs->vel.x;
	uC = A->get_ifi(i,j)->fs->vel.x;
	uD = A->get_ifi(i,j-1)->fs->vel.x;
	//
        vA = A->get_cell(i,j-1)->fs->vel.y;
        vB = A->get_cell(i,j)->fs->vel.y;
	vC = A->get_ifi(i,j)->fs->vel.y;
	vD = A->get_ifi(i,j-1)->fs->vel.y;
	//
	for ( size_t itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_cell(i,j-1)->fs->gas->T[itm];
	    TB[itm] = A->get_cell(i,j)->fs->gas->T[itm];
	    TC[itm] = A->get_ifi(i,j)->fs->gas->T[itm];
	    TD[itm] = A->get_ifi(i,j-1)->fs->gas->T[itm];
	}
        //
	if( G.diffusion ) {
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_cell(i,j-1)->fs->gas->massf[isp];
		fB[isp] = A->get_cell(i,j)->fs->gas->massf[isp];
		fC[isp] = A->get_ifi(i,j)->fs->gas->massf[isp];
		fD[isp] = A->get_ifi(i,j-1)->fs->gas->massf[isp];
	    }
	}
	//
	if ( G.electric_field_work ) {
	    peA = A->get_cell(i,j-1)->fs->gas->p_e;
	    peB = A->get_cell(i,j)->fs->gas->p_e;
	    peC = A->get_ifi(i,j)->fs->gas->p_e;
	    peD = A->get_ifi(i,j-1)->fs->gas->p_e;
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
        tkeA = A->get_cell(i,j-1)->fs->tke;
        tkeB = A->get_cell(i,j)->fs->tke;
	tkeC = A->get_ifi(i,j)->fs->tke;
	tkeD = A->get_ifi(i,j-1)->fs->tke;
	//
        omegaA = A->get_cell(i,j-1)->fs->omega;
        omegaB = A->get_cell(i,j)->fs->omega;
	omegaC = A->get_ifi(i,j)->fs->omega;
	omegaD = A->get_ifi(i,j-1)->fs->omega;
	//
	APPLY_DIVERGENCE_THEOREM_2()
	//
    } // j loop

    // NORTH boundary
    j = A->jmax+1;
    for (i = A->imin+1; i <= A->imax; ++i) {
        area_inv = 1.0 / A->get_vtx(i,j)->area;
	// These are the corners of the secondary cell.
        xA = A->get_cell(i,j-1)->pos[gtl].x;
        yA = A->get_cell(i,j-1)->pos[gtl].y;
	xB = A->get_ifj(i,j)->pos.x;
	yB = A->get_ifj(i,j)->pos.y;
	xC = A->get_ifj(i-1,j)->pos.x;
	yC = A->get_ifj(i-1,j)->pos.y;
        xD = A->get_cell(i-1,j-1)->pos[gtl].x;
        yD = A->get_cell(i-1,j-1)->pos[gtl].y;
	// These are the flow properties at the corners of
	// the secondary cell.
        uA = A->get_cell(i,j-1)->fs->vel.x;
	uB = A->get_ifj(i,j)->fs->vel.x;
	uC = A->get_ifj(i-1,j)->fs->vel.x;
        uD = A->get_cell(i-1,j-1)->fs->vel.x;
	//
        vA = A->get_cell(i,j-1)->fs->vel.y;
	vB = A->get_ifj(i,j)->fs->vel.y;
	vC = A->get_ifj(i-1,j)->fs->vel.y;
        vD = A->get_cell(i-1,j-1)->fs->vel.y;
	//
	for ( size_t itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_cell(i,j-1)->fs->gas->T[itm];
	    TB[itm] = A->get_ifj(i,j)->fs->gas->T[itm];
	    TC[itm] = A->get_ifj(i-1,j)->fs->gas->T[itm];
	    TD[itm] = A->get_cell(i-1,j-1)->fs->gas->T[itm];
	}
	//
	if( G.diffusion ) {
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_cell(i,j-1)->fs->gas->massf[isp];
		fB[isp] = A->get_ifj(i,j)->fs->gas->massf[isp];
		fC[isp] = A->get_ifj(i-1,j)->fs->gas->massf[isp];
		fD[isp] = A->get_cell(i-1,j-1)->fs->gas->massf[isp];
	    }
	}
	//
	if ( G.electric_field_work ) {
	    peA = A->get_cell(i,j-1)->fs->gas->p_e;
	    peB = A->get_ifi(i,j)->fs->gas->p_e;
	    peC = A->get_ifi(i-1,j)->fs->gas->p_e;
	    peD = A->get_cell(i-1,j-1)->fs->gas->p_e;
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
        tkeA = A->get_cell(i,j-1)->fs->tke;
	tkeB = A->get_ifj(i,j)->fs->tke;
	tkeC = A->get_ifj(i-1,j)->fs->tke;
        tkeD = A->get_cell(i-1,j-1)->fs->tke;
	//
        omegaA = A->get_cell(i,j-1)->fs->omega;
	omegaB = A->get_ifj(i,j)->fs->omega;
	omegaC = A->get_ifj(i-1,j)->fs->omega;
        omegaD = A->get_cell(i-1,j-1)->fs->omega;
	//
	APPLY_DIVERGENCE_THEOREM_2()
	//
    } // i loop

    // SOUTH boundary
    j = A->jmin;
    for (i = A->imin+1; i <= A->imax; ++i) {
        area_inv = 1.0 / A->get_vtx(i,j)->area;
	// These are the corners of the secondary cell.
	xA = A->get_ifj(i,j)->pos.x;
	yA = A->get_ifj(i,j)->pos.y;
        xB = A->get_cell(i,j)->pos[gtl].x;
        yB = A->get_cell(i,j)->pos[gtl].y;
        xC = A->get_cell(i-1,j)->pos[gtl].x;
        yC = A->get_cell(i-1,j)->pos[gtl].y;
	xD = A->get_ifj(i-1,j)->pos.x;
	yD = A->get_ifj(i-1,j)->pos.y;
	// These are the flow properties at the corners of the secondary cell.
	uA = A->get_ifj(i,j)->fs->vel.x;
        uB = A->get_cell(i,j)->fs->vel.x;
        uC = A->get_cell(i-1,j)->fs->vel.x;
	uD = A->get_ifj(i-1,j)->fs->vel.x;
	//
	vA = A->get_ifj(i,j)->fs->vel.y;
        vB = A->get_cell(i,j)->fs->vel.y;
        vC = A->get_cell(i-1,j)->fs->vel.y;
	vD = A->get_ifj(i-1,j)->fs->vel.y;
	//
	for ( size_t itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_ifj(i,j)->fs->gas->T[itm];
	    TB[itm] = A->get_cell(i,j)->fs->gas->T[itm];
	    TC[itm] = A->get_cell(i-1,j)->fs->gas->T[itm];
	    TD[itm] = A->get_ifj(i-1,j)->fs->gas->T[itm];
	}
	//
	if( G.diffusion ) { 
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_ifj(i,j)->fs->gas->massf[isp];
		fB[isp] = A->get_cell(i,j)->fs->gas->massf[isp];
		fC[isp] = A->get_cell(i-1,j)->fs->gas->massf[isp];
		fD[isp] = A->get_ifj(i-1,j)->fs->gas->massf[isp];
	    }
	}
	//
	if ( G.electric_field_work ) {
	    peA = A->get_ifi(i,j)->fs->gas->p_e;
	    peB = A->get_cell(i,j)->fs->gas->p_e;
	    peC = A->get_cell(i-1,j)->fs->gas->p_e;
	    peD = A->get_ifi(i-1,j-1)->fs->gas->p_e;
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
	tkeA = A->get_ifj(i,j)->fs->tke;
        tkeB = A->get_cell(i,j)->fs->tke;
        tkeC = A->get_cell(i-1,j)->fs->tke;
	tkeD = A->get_ifj(i-1,j)->fs->tke;
	//
	omegaA = A->get_ifj(i,j)->fs->omega;
        omegaB = A->get_cell(i,j)->fs->omega;
        omegaC = A->get_cell(i-1,j)->fs->omega;
	omegaD = A->get_ifj(i-1,j)->fs->omega;
	//
	APPLY_DIVERGENCE_THEOREM_2()
	//
    } // i loop
    return SUCCESS;
} // end viscous_derivatives_edges()


/// \brief Compute the derivatives of velocity, temperature and mass fractions at
///        primary cell vertices at the corners of the block.
///
int viscous_derivatives_corners(Block *bdp, size_t gtl)
{
    size_t i, j;
    FV_Vertex *vtx;
    FV_Interface *a, *b;
    FV_Cell *c;
    double fa, xa, ya, fb, xb, yb, fc, xc, yc, denom;
    size_t nsp = get_gas_model_ptr()->get_number_of_species();
    size_t ntm = get_gas_model_ptr()->get_number_of_modes();
    // North-East corner
    i = bdp->imax;
    j = bdp->jmax;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i+1,j+1);
    a = bdp->get_ifj(i,j+1);
    b = bdp->get_ifi(i+1,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos[gtl].x; yc = c->pos[gtl].y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( size_t itm=0; itm<ntm; ++itm ) {
    	 fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm];
    	 vtx->dTdx[itm] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    	 vtx->dTdy[itm] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke;
    vtx->dtkedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dtkedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega;
    vtx->domegadx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->domegady = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for( size_t isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->gas->p_e; fb = b->fs->gas->p_e; fc = c->fs->gas->p_e;
    vtx->dpedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dpedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    // South-East corner
    i = bdp->imax;
    j = bdp->jmin;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i+1,j);
    a = bdp->get_ifj(i,j);
    b = bdp->get_ifi(i+1,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos[gtl].x; yc = c->pos[gtl].y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( size_t itm=0; itm<ntm; ++itm ) {
    	fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm];
    	vtx->dTdx[itm] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    	vtx->dTdy[itm] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke;
    vtx->dtkedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dtkedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega;
    vtx->domegadx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->domegady = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for( size_t isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->gas->p_e; fb = b->fs->gas->p_e; fc = c->fs->gas->p_e;
    vtx->dpedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dpedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    // South-West corner
    i = bdp->imin;
    j = bdp->jmin;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i,j);
    a = bdp->get_ifj(i,j);
    b = bdp->get_ifi(i,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos[gtl].x; yc = c->pos[gtl].y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( size_t itm=0; itm<ntm; ++itm ) {
    	fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm];
    	vtx->dTdx[itm] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    	vtx->dTdy[itm] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke;
    vtx->dtkedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dtkedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega;
    vtx->domegadx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->domegady = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for( size_t isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->gas->p_e; fb = b->fs->gas->p_e; fc = c->fs->gas->p_e;
    vtx->dpedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dpedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    // North-West corner
    i = bdp->imin;
    j = bdp->jmax;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i,j+1);
    a = bdp->get_ifj(i,j+1);
    b = bdp->get_ifi(i,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos[gtl].x; yc = c->pos[gtl].y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( size_t itm=0; itm<ntm; ++itm ) {
    	fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm];
    	vtx->dTdx[itm] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    	vtx->dTdy[itm] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke;
    vtx->dtkedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dtkedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega;
    vtx->domegadx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->domegady = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for( size_t isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    fa = a->fs->gas->p_e; fb = b->fs->gas->p_e; fc = c->fs->gas->p_e;
    vtx->dpedx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dpedy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    return SUCCESS;
} // end viscous_derivatives_corners()
