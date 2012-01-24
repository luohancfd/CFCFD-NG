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
 * 02-Mar-2008: Elmer3 port
 *
 */

/*-----------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "block.hh"
#include "kernel.hh"
#include "bc_defs.hh"
#include "bc.hh"
#include "visc.hh"
#include "diffusion.hh"
#include "baldwin_lomax.hh"

#define VERY_SMALL  1.0e-10

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
int estimate_turbulence_viscosity( global_data *gdp, Block *bdp )
{
    if ( get_turbulence_flag() == 0 ) {
	bdp->apply( &FV_Cell::turbulence_viscosity_zero, "turbulence_viscosity_zero" );
	return SUCCESS;
    }

    if ( get_k_omega_flag() == 1 ) {
	bdp->apply( &FV_Cell::turbulence_viscosity_k_omega, "turbulence_viscosity_k_omega" );
    } else if ( get_baldwin_lomax_flag() == 1 ) {
	baldwin_lomax_turbulence_model( *gdp, *bdp );
    } else {
	cout << "Turbulence model requested but not available." << endl;
	exit( NOT_IMPLEMENTED_ERROR );
    }
    bdp->apply( &FV_Cell::turbulence_viscosity_factor, gdp->transient_mu_t_factor, 
		"turbulence_viscosity_factor" );
    bdp->apply( &FV_Cell::turbulence_viscosity_limit, gdp->max_mu_t_factor, 
		"turbulence_viscosity_limit" );
    bdp->apply( &FV_Cell::turbulence_viscosity_zero_if_not_in_zone,
		"turbulence_viscosity_zero_if_not_in_zone" );
    return SUCCESS;
}

/// \brief Compute the viscous contribution to the cell interface fluxes.
///
/// This contribution is added to the flux variables so make sure that
/// the inviscid (Euler) contributions have already been computed and stored.
///
int viscous_flux_2D( Block *A )
{
    Gas_model *gmodel = get_gas_model_ptr();
    FV_Vertex *vtx1, *vtx2;
    FV_Interface *IFace;
    double nx, ny;
    int i, j, nsp, ntm;
    double dudx, dudy, dvdx, dvdy;
    double mu_lam; // molecular-scale transport properties
    double lmbda; // second-coefficient of viscosity
    double mu_t;  // turbulence viscosity
    double k_t;   // turbulence thermal conductivity
    double viscous_factor; // so that we can scale down the viscous effects
    double mu_eff; // combined laminar and turbulent viscosity after scaling
    double dtkedx, dtkedy, domegadx, domegady;
    double sigma = 0.5;
    double sigma_star = 0.6;
    double mu_effective;
    double tau_kx, tau_ky, tau_wx, tau_wy;
    double tau_xx, tau_yy, tau_xy;
    double ybar;

    nsp = gmodel->get_number_of_species();
    if( dfdx.size() == 0 ) {
	dfdx.resize(nsp); 
	dfdy.resize(nsp);
	dfdz.resize(nsp);
	jx.resize(nsp);
	jy.resize(nsp);
	jz.resize(nsp);
    }
    
    ntm = gmodel->get_number_of_modes();
    if( dTdx.size() == 0 ) {
	dTdx.resize(ntm);
	dTdy.resize(ntm);
	qx.resize(ntm);
	qy.resize(ntm);
	k_eff.resize(ntm);
    }
    
    viscous_factor = get_viscous_factor();

    // East-facing interfaces.
    for (i = A->imin; i <= A->imax+1; ++i) {
        for (j = A->jmin; j <= A->jmax; ++j) {
            IFace = A->get_ifi(i,j);
	    FlowState &fs = *(IFace->fs);
            vtx1 = A->get_vtx(i,j+1);
            vtx2 = A->get_vtx(i,j);
	    // Determine some of the interface properties.
            ybar = IFace->Ybar;
            dudx = 0.5 * (vtx1->dudx + vtx2->dudx);
            dudy = 0.5 * (vtx1->dudy + vtx2->dudy);
            dvdx = 0.5 * (vtx1->dvdx + vtx2->dvdx);
            dvdy = 0.5 * (vtx1->dvdy + vtx2->dvdy);
            for ( int itm=0; itm<ntm; ++itm ) {
            	dTdx[itm] = 0.5 * (vtx1->dTdx[itm] + vtx2->dTdx[itm]);
            	dTdy[itm] = 0.5 * (vtx1->dTdy[itm] + vtx2->dTdy[itm]);
            	k_eff[itm] = viscous_factor * fs.gas->k[itm];
            }
            mu_lam = viscous_factor * fs.gas->mu;
	    mu_t = viscous_factor * fs.mu_t;
	    k_t = viscous_factor * fs.k_t;
	    // CHECK-ME: Assume turbulent conductivity only applies to primary mode?
	    k_eff[0] += k_t;
	    mu_eff = mu_lam + mu_t;
	    lmbda = -2.0/3.0 * mu_eff;
	    if( get_diffusion_flag() == 1 ) {
		for( int isp = 0; isp < nsp; ++isp ) {
		    dfdx[isp] = 0.5 * (vtx1->dfdx[isp] + vtx2->dfdx[isp]);
		    dfdy[isp] = 0.5 * (vtx1->dfdy[isp] + vtx2->dfdy[isp]);
		}
		// Apply a diffusion model
		double D_t = 0.0;
		if ( get_k_omega_flag() == 1 ) {
                    double Sc_t = get_turbulence_schmidt_number();
                    D_t = mu_t / (fs.gas->rho * Sc_t);
		}
		calculate_diffusion_fluxes(*(fs.gas),
					   D_t,
					   dfdx, dfdy, dfdz,
					   jx, jy, jz);
		// NOTE: now applying viscous_factor to diffusive fluxes
		for( int isp = 0; isp < nsp; ++isp ) {
		    jx[isp] *= viscous_factor;
		    jy[isp] *= viscous_factor;
		    jz[isp] *= viscous_factor;
		}
	    }
	    
	    if ( get_axisymmetric_flag() == 1 ) {
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
	    for ( int itm=1; itm<ntm; ++itm ) {
		qx[itm] = k_eff[itm] * dTdx[itm];
		qy[itm] = k_eff[itm] * dTdy[itm];
		qx[0] += qx[itm];
		qy[0] += qy[itm];
	    }
	    
	    if( get_diffusion_flag() == 1 ) {
		for( int isp = 0; isp < nsp; ++isp ) {
		    double h = gmodel->enthalpy(*(fs.gas), isp);
		    qx[0] -= jx[isp] * h;
		    qy[0] -= jy[isp] * h;
		    for ( int itm=1; itm<ntm; ++itm ) {
			double hmode = gmodel->modal_enthalpy(fs.gas->T[itm], isp, itm);
			qx[itm] -= jx[isp] * hmode;
			qy[itm] -= jy[isp] * hmode;
		    }
		}
	    }	    
	    
	    if ( get_k_omega_flag() == 1 ) {
		// Turbulence contribution to the shear stresses.
		tau_xx -= 0.66667 * fs.gas->rho * fs.tke;
		tau_yy -= 0.66667 * fs.gas->rho * fs.tke;
		// Turbulence contribution to heat transfer.
		mu_effective = mu_lam + sigma_star * mu_t;
		dtkedx = 0.5 * (vtx1->dtkedx + vtx2->dtkedx);
		dtkedy = 0.5 * (vtx1->dtkedy + vtx2->dtkedy);
		qx[0] += mu_effective * dtkedx;
		qy[0] += mu_effective * dtkedy;
		// Turbulence transport of the turbulence properties themselves.
		tau_kx = mu_effective * dtkedx; 
		tau_ky = mu_effective * dtkedy;
		mu_effective = mu_lam + sigma * mu_t;
		domegadx = 0.5 * (vtx1->domegadx + vtx2->domegadx);
		domegady = 0.5 * (vtx1->domegady + vtx2->domegady);
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
            F.momentum.x -= tau_xx * nx + tau_xy * ny;
            F.momentum.y -= tau_xy * nx + tau_yy * ny;
            F.total_energy -= (tau_xx * fs.vel.x + tau_xy * fs.vel.y + qx[0]) * nx
                + (tau_xy * fs.vel.x + tau_yy * fs.vel.y + qy[0]) * ny;
	    // Viscous transport of k-omega turbulence quantities.
	    // Only built for 2D planar geometry at the moment.
	    if ( get_k_omega_flag() == 1 ) {
		F.tke -= tau_kx * nx + tau_ky * ny;
		F.omega -= tau_wx * nx + tau_wy * ny;
	    }
            // Species mass flux
	    if( get_diffusion_flag() == 1 ) {
		for( int isp = 0; isp < nsp; ++isp ) {
		    F.massf[isp] += jx[isp] * nx + jy[isp] * ny;
		}
	    }
	    // Modal energy flux (skipping first mode as this is handled by total energy)
	    for ( int itm=1; itm<ntm; ++itm ) {
	    	F.energies[itm] -= qx[itm] * nx + qy[itm] * ny;
	    }
        } // j loop
    } // i loop
    /*
     * North-facing interfaces
     */
    for (i = A->imin; i <= A->imax; ++i) {
        for (j = A->jmin; j <= A->jmax+1; ++j) {
            IFace = A->get_ifj(i,j);
	    FlowState &fs = *(IFace->fs);
            vtx1 = A->get_vtx(i,j);
            vtx2 = A->get_vtx(i+1,j);
	    // Determine some of the interface properties.
            ybar = IFace->Ybar;
            dudx = 0.5 * (vtx1->dudx + vtx2->dudx);
            dudy = 0.5 * (vtx1->dudy + vtx2->dudy);
            dvdx = 0.5 * (vtx1->dvdx + vtx2->dvdx);
            dvdy = 0.5 * (vtx1->dvdy + vtx2->dvdy);
            for ( int itm=0; itm<ntm; ++itm ) {
            	dTdx[itm] = 0.5 * (vtx1->dTdx[itm] + vtx2->dTdx[itm]);
            	dTdy[itm] = 0.5 * (vtx1->dTdy[itm] + vtx2->dTdy[itm]);
            	k_eff[itm] = viscous_factor * fs.gas->k[itm];
            }
            mu_lam = viscous_factor * fs.gas->mu;
	    mu_t = viscous_factor * fs.mu_t;
	    k_t = viscous_factor * fs.k_t;
	    // CHECK-ME: Assume turbulent conductivity only applies to primary mode?
	    k_eff[0] += k_t;
	    mu_eff = mu_lam + mu_t;
	    lmbda = -2.0/3.0 * mu_eff;
	    if( get_diffusion_flag() == 1 ) {
		for( int isp = 0; isp < nsp; ++isp ) {
		    dfdx[isp] = 0.5 * (vtx1->dfdx[isp] + vtx2->dfdx[isp]);
		    dfdy[isp] = 0.5 * (vtx1->dfdy[isp] + vtx2->dfdy[isp]);
		}
		// Apply a diffusion model
		double D_t = 0.0;
		if ( get_k_omega_flag() == 1 ) {
                    double Sc_t = get_turbulence_schmidt_number();
                    D_t = mu_t / (fs.gas->rho * Sc_t);
		}
		calculate_diffusion_fluxes(*(fs.gas),
					   D_t,
					   dfdx, dfdy, dfdz,
					   jx, jy, jz);
		// NOTE: now applying viscous_factor to diffusive fluxes
		for( int isp = 0; isp < nsp; ++isp ) {
		    jx[isp] *= viscous_factor;
		    jy[isp] *= viscous_factor;
		    jz[isp] *= viscous_factor;
		}
	    }

	    if ( get_axisymmetric_flag() == 1 ) {
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
	    for ( int itm=1; itm<ntm; ++itm ) {
		qx[itm] = k_eff[itm] * dTdx[itm];
		qy[itm] = k_eff[itm] * dTdy[itm];
		qx[0] += qx[itm];
		qy[0] += qy[itm];
	    }
	    
	    if( get_diffusion_flag() == 1 ) {
		for( int isp = 0; isp < nsp; ++isp ) {
		    double h = gmodel->enthalpy(*(fs.gas), isp);
		    qx[0] -= jx[isp] * h;
		    qy[0] -= jy[isp] * h;
		    for ( int itm=1; itm<ntm; ++itm ) {
			double hmode = gmodel->modal_enthalpy(fs.gas->T[itm], isp, itm);
			qx[itm] -= jx[isp] * hmode;
			qy[itm] -= jy[isp] * hmode;
		    }
		}
	    }
	    
	    if ( get_k_omega_flag() == 1 ) {
		// Turbulence contribution to the shear stresses.
		tau_xx -= 0.66667 * fs.gas->rho * fs.tke;
		tau_yy -= 0.66667 * fs.gas->rho * fs.tke;
		// Turbulence contribution to heat transfer.
		mu_effective = mu_lam + sigma_star * mu_t;
		dtkedx = 0.5 * (vtx1->dtkedx + vtx2->dtkedx);
		dtkedy = 0.5 * (vtx1->dtkedy + vtx2->dtkedy);
		qx[0] += mu_effective * dtkedx;
		qy[0] += mu_effective * dtkedy;
		// Turbulence transport of the turbulence properties themselves.
		tau_kx = mu_effective * dtkedx; 
		tau_ky = mu_effective * dtkedy;
		mu_effective = mu_lam + sigma * mu_t;
		domegadx = 0.5 * (vtx1->domegadx + vtx2->domegadx);
		domegady = 0.5 * (vtx1->domegady + vtx2->domegady);
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
            F.momentum.x -= tau_xx * nx + tau_xy * ny;
            F.momentum.y -= tau_xy * nx + tau_yy * ny;
            F.total_energy -= (tau_xx * fs.vel.x + tau_xy * fs.vel.y + qx[0]) * nx
                + (tau_xy * fs.vel.x + tau_yy * fs.vel.y + qy[0]) * ny;
	    // Viscous transport of k-omega turbulence quantities.
	    // Only built for 2D planar geometry at the moment.
	    if ( get_k_omega_flag() == 1 ) {
		F.tke -= tau_kx * nx + tau_ky * ny;
		F.omega -= tau_wx * nx + tau_wy * ny;
	    }
	    // Species mass flux
	    if( get_diffusion_flag() == 1 ) {
		for( int isp = 0; isp < nsp; ++isp ) {
		    F.massf[isp] += jx[isp] * nx + jy[isp] * ny;
		}
	    }
	    // Modal energy flux (skipping first mode as this is handled by total energy)
	    for ( int itm=1; itm<ntm; ++itm ) {
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
	    for ( int itm=0; itm<ntm; ++itm ) { \
	        A->get_vtx(i,j)->dTdx[itm] = 0.5 * area_inv * \
		    ((TB[itm] + TA[itm]) * (yB - yA) + (TC[itm] + TB[itm]) * (yC - yB) + \
		     (TD[itm] + TC[itm]) * (yD - yC) + (TA[itm] + TD[itm]) * (yA - yD)); \
                A->get_vtx(i,j)->dTdy[itm] = -0.5 * area_inv * \
		    ((TB[itm] + TA[itm]) * (xB - xA) + (TC[itm] + TB[itm]) * (xC - xB) + \
		     (TD[itm] + TC[itm]) * (xD - xC) + (TA[itm] + TD[itm]) * (xA - xD)); \
            } \
            if( get_diffusion_flag() == 1) { \
                for( int isp = 0; isp < nsp; ++isp ) { \
                    A->get_vtx(i,j)->dfdx[isp] = 0.5 * area_inv * \
		        ((fB[isp] + fA[isp]) * (yB - yA) + (fC[isp] + fB[isp]) * (yC - yB) + \
		        (fD[isp] + fC[isp]) * (yD - yC) + (fA[isp] + fD[isp]) * (yA - yD)); \
                    A->get_vtx(i,j)->dfdy[isp] = -0.5 * area_inv * \
		        ((fB[isp] + fA[isp]) * (xB - xA) + (fC[isp] + fB[isp]) * (xC - xB) + \
		        (fD[isp] + fC[isp]) * (xD - xC) + (fA[isp] + fD[isp]) * (xA - xD)); \
	        } \
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
           for( int isp = 0; isp < nsp; ++isp ) { \
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
int viscous_derivatives_2D( Block *A )
{
    int i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double uA, uB, uC, uD;
    double vA, vB, vC, vD;
    double tkeA, tkeB, tkeC, tkeD;
    double omegaA, omegaB, omegaC, omegaD;
    double area_inv;
    
    int nsp = get_gas_model_ptr()->get_number_of_species();
    if( fA.size() == 0 ) {
	fA.resize(nsp); 
	fB.resize(nsp);
	fC.resize(nsp);
	fD.resize(nsp);
    }
    
    int ntm = get_gas_model_ptr()->get_number_of_modes();
    if ( TA.size() == 0 ) {
    	TA.resize(ntm);
    	TB.resize(ntm);
    	TC.resize(ntm);
    	TD.resize(ntm);
    }

    // First, do all of the internal secondary cells.
    // i.e. Those not on a boundary.
    for (i = A->imin+1; i <= A->imax; ++i) {
        for (j = A->jmin+1; j <= A->jmax; ++j) {
            area_inv = 1.0 / A->get_vtx(i,j)->area;
	    // These are the corners of the secondary cell.
            xA = A->get_cell(IADSH,JADSH)->pos.x;
            yA = A->get_cell(IADSH,JADSH)->pos.y;
            xB = A->get_cell(IBDSH,JBDSH)->pos.x;
            yB = A->get_cell(IBDSH,JBDSH)->pos.y;
            xC = A->get_cell(ICDSH,JCDSH)->pos.x;
            yC = A->get_cell(ICDSH,JCDSH)->pos.y;
            xD = A->get_cell(IDDSH,JDDSH)->pos.x;
            yD = A->get_cell(IDDSH,JDDSH)->pos.y;
	    // These are the flow properties at the corners.
            uA = A->get_cell(IADSH,JADSH)->fs->vel.x;
            uB = A->get_cell(IBDSH,JBDSH)->fs->vel.x;
            uC = A->get_cell(ICDSH,JCDSH)->fs->vel.x;
            uD = A->get_cell(IDDSH,JDDSH)->fs->vel.x;
	    //
            vA = A->get_cell(IADSH,JADSH)->fs->vel.y;
            vB = A->get_cell(IBDSH,JBDSH)->fs->vel.y;
            vC = A->get_cell(ICDSH,JCDSH)->fs->vel.y;
            vD = A->get_cell(IDDSH,JDDSH)->fs->vel.y;
	    //
	    for ( int itm=0; itm<ntm; ++itm ) {
		TA[itm] = A->get_cell(IADSH,JADSH)->fs->gas->T[itm];
		TB[itm] = A->get_cell(IBDSH,JBDSH)->fs->gas->T[itm];
		TC[itm] = A->get_cell(ICDSH,JCDSH)->fs->gas->T[itm];
		TD[itm] = A->get_cell(IDDSH,JDDSH)->fs->gas->T[itm];
            }
	    //
	    if( get_diffusion_flag() == 1 ) {
		for( int isp = 0; isp < nsp; ++isp ) {
		    fA[isp] = A->get_cell(IADSH,JADSH)->fs->gas->massf[isp];
		    fB[isp] = A->get_cell(IBDSH,JBDSH)->fs->gas->massf[isp];
		    fC[isp] = A->get_cell(ICDSH,JCDSH)->fs->gas->massf[isp];
		    fD[isp] = A->get_cell(IDDSH,JDDSH)->fs->gas->massf[isp];
		}
	    }
	    //
	    APPLY_DIVERGENCE_THEOREM()
	    //
            tkeA = A->get_cell(IADSH,JADSH)->fs->tke;
            tkeB = A->get_cell(IBDSH,JBDSH)->fs->tke;
            tkeC = A->get_cell(ICDSH,JCDSH)->fs->tke;
            tkeD = A->get_cell(IDDSH,JDDSH)->fs->tke;
	    //
            omegaA = A->get_cell(IADSH,JADSH)->fs->omega;
            omegaB = A->get_cell(IBDSH,JBDSH)->fs->omega;
            omegaC = A->get_cell(ICDSH,JCDSH)->fs->omega;
            omegaD = A->get_cell(IDDSH,JDDSH)->fs->omega;
	    //
	    APPLY_DIVERGENCE_THEOREM_2()
        } // j loop
    } // i loop

    // Now, do the boundaries as half cells.
    viscous_derivatives_edges(A);
    // ...and pick up (fudge) the corner values.
    viscous_derivatives_corners(A);
    return SUCCESS;
} // end viscous_derivatives()


/// \brief Compute the derivatives of velocity, temperature and mass fractions at
///        primary cell vertices along the block boundaries.
///
/// NOTE that the secondary cells along boundaries are HALF cells.
///
int viscous_derivatives_edges( Block *A )
{
    int i, j;
    double xA, yA, xB, yB, xC, yC, xD, yD;
    double uA, uB, uC, uD;
    double vA, vB, vC, vD;
    double tkeA, tkeB, tkeC, tkeD;
    double omegaA, omegaB, omegaC, omegaD;
    double area_inv;

    int nsp = get_gas_model_ptr()->get_number_of_species();
    if( fA.size() == 0 ) {
	fA.resize(nsp); 
	fB.resize(nsp);
	fC.resize(nsp);
	fD.resize(nsp);
    }
    
    int ntm = get_gas_model_ptr()->get_number_of_modes();
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
        xC = A->get_cell(ICDSH,JCDSH)->pos.x;
        yC = A->get_cell(ICDSH,JCDSH)->pos.y;
        xD = A->get_cell(IDDSH,JDDSH)->pos.x;
        yD = A->get_cell(IDDSH,JDDSH)->pos.y;
	//
	// These are the flow properties at the corners of the secondary cell.
	uA = A->get_ifi(i,j-1)->fs->vel.x;
	uB = A->get_ifi(i,j)->fs->vel.x;
        uC = A->get_cell(ICDSH,JCDSH)->fs->vel.x;
        uD = A->get_cell(IDDSH,JDDSH)->fs->vel.x;
	//
	vA = A->get_ifi(i,j-1)->fs->vel.y;
	vB = A->get_ifi(i,j)->fs->vel.y;
        vC = A->get_cell(ICDSH,JCDSH)->fs->vel.y;
        vD = A->get_cell(IDDSH,JDDSH)->fs->vel.y;
	//
	for ( int itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_ifi(i,j-1)->fs->gas->T[itm];
	    TB[itm] = A->get_ifi(i,j)->fs->gas->T[itm];
	    TC[itm] = A->get_cell(ICDSH,JCDSH)->fs->gas->T[itm];
	    TD[itm] = A->get_cell(IDDSH,JDDSH)->fs->gas->T[itm];
	}
	//
	if( get_diffusion_flag() == 1) {
	    for( int isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_ifi(i,j-1)->fs->gas->massf[isp];
		fB[isp] = A->get_ifi(i,j)->fs->gas->massf[isp];
		fC[isp] = A->get_cell(ICDSH,JCDSH)->fs->gas->massf[isp];
		fD[isp] = A->get_cell(IDDSH,JDDSH)->fs->gas->massf[isp];
	    }
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
	tkeA = A->get_ifi(i,j-1)->fs->tke;
	tkeB = A->get_ifi(i,j)->fs->tke;
	tkeC = A->get_cell(ICDSH,JCDSH)->fs->tke;
	tkeD = A->get_cell(IDDSH,JDDSH)->fs->tke;
	//
	omegaA = A->get_ifi(i,j-1)->fs->omega;
	omegaB = A->get_ifi(i,j)->fs->omega;
	omegaC = A->get_cell(ICDSH,JCDSH)->fs->omega;
	omegaD = A->get_cell(IDDSH,JDDSH)->fs->omega;
	//
	APPLY_DIVERGENCE_THEOREM_2()
	//
    } // j loop

    // WEST boundary
    i = A->imin;
    for (j = A->jmin+1; j <= A->jmax; ++j) {
        area_inv = 1.0 / A->get_vtx(i,j)->area;
	// These are the corners of the secondary cell.
        xA = A->get_cell(IADSH,JADSH)->pos.x;
        yA = A->get_cell(IADSH,JADSH)->pos.y;
        xB = A->get_cell(IBDSH,JBDSH)->pos.x;
        yB = A->get_cell(IBDSH,JBDSH)->pos.y;
	xC = A->get_ifi(i,j)->pos.x;
	yC = A->get_ifi(i,j)->pos.y;
	xD = A->get_ifi(i,j-1)->pos.x;
	yD = A->get_ifi(i,j-1)->pos.y;
	// These are the flow properties at the corners of the secondary cell.
        uA = A->get_cell(IADSH,JADSH)->fs->vel.x;
        uB = A->get_cell(IBDSH,JBDSH)->fs->vel.x;
	uC = A->get_ifi(i,j)->fs->vel.x;
	uD = A->get_ifi(i,j-1)->fs->vel.x;
	//
        vA = A->get_cell(IADSH,JADSH)->fs->vel.y;
        vB = A->get_cell(IBDSH,JBDSH)->fs->vel.y;
	vC = A->get_ifi(i,j)->fs->vel.y;
	vD = A->get_ifi(i,j-1)->fs->vel.y;
	//
	for ( int itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_cell(IADSH,JADSH)->fs->gas->T[itm];
	    TB[itm] = A->get_cell(IBDSH,JBDSH)->fs->gas->T[itm];
	    TC[itm] = A->get_ifi(i,j)->fs->gas->T[itm];
	    TD[itm] = A->get_ifi(i,j-1)->fs->gas->T[itm];
	}
        //
	if( get_diffusion_flag() == 1) {
	    for( int isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_cell(IADSH,JADSH)->fs->gas->massf[isp];
		fB[isp] = A->get_cell(IBDSH,JBDSH)->fs->gas->massf[isp];
		fC[isp] = A->get_ifi(i,j)->fs->gas->massf[isp];
		fD[isp] = A->get_ifi(i,j-1)->fs->gas->massf[isp];
	    }
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
        tkeA = A->get_cell(IADSH,JADSH)->fs->tke;
        tkeB = A->get_cell(IBDSH,JBDSH)->fs->tke;
	tkeC = A->get_ifi(i,j)->fs->tke;
	tkeD = A->get_ifi(i,j-1)->fs->tke;
	//
        omegaA = A->get_cell(IADSH,JADSH)->fs->omega;
        omegaB = A->get_cell(IBDSH,JBDSH)->fs->omega;
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
        xA = A->get_cell(IADSH,JADSH)->pos.x;
        yA = A->get_cell(IADSH,JADSH)->pos.y;
	xB = A->get_ifj(i,j)->pos.x;
	yB = A->get_ifj(i,j)->pos.y;
	xC = A->get_ifj(i-1,j)->pos.x;
	yC = A->get_ifj(i-1,j)->pos.y;
        xD = A->get_cell(IDDSH,JDDSH)->pos.x;
        yD = A->get_cell(IDDSH,JDDSH)->pos.y;
	// These are the flow properties at the corners of
	// the secondary cell.
        uA = A->get_cell(IADSH,JADSH)->fs->vel.x;
	uB = A->get_ifj(i,j)->fs->vel.x;
	uC = A->get_ifj(i-1,j)->fs->vel.x;
        uD = A->get_cell(IDDSH,JDDSH)->fs->vel.x;
	//
        vA = A->get_cell(IADSH,JADSH)->fs->vel.y;
	vB = A->get_ifj(i,j)->fs->vel.y;
	vC = A->get_ifj(i-1,j)->fs->vel.y;
        vD = A->get_cell(IDDSH,JDDSH)->fs->vel.y;
	//
	for ( int itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_cell(IADSH,JADSH)->fs->gas->T[itm];
	    TB[itm] = A->get_ifj(i,j)->fs->gas->T[itm];
	    TC[itm] = A->get_ifj(i-1,j)->fs->gas->T[itm];
	    TD[itm] = A->get_cell(IDDSH,JDDSH)->fs->gas->T[itm];
	}
	//
	if( get_diffusion_flag() == 1) {
	    for( int isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_cell(IADSH,JADSH)->fs->gas->massf[isp];
		fB[isp] = A->get_ifj(i,j)->fs->gas->massf[isp];
		fC[isp] = A->get_ifj(i-1,j)->fs->gas->massf[isp];
		fD[isp] = A->get_cell(IDDSH,JDDSH)->fs->gas->massf[isp];
	    }
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
        tkeA = A->get_cell(IADSH,JADSH)->fs->tke;
	tkeB = A->get_ifj(i,j)->fs->tke;
	tkeC = A->get_ifj(i-1,j)->fs->tke;
        tkeD = A->get_cell(IDDSH,JDDSH)->fs->tke;
	//
        omegaA = A->get_cell(IADSH,JADSH)->fs->omega;
	omegaB = A->get_ifj(i,j)->fs->omega;
	omegaC = A->get_ifj(i-1,j)->fs->omega;
        omegaD = A->get_cell(IDDSH,JDDSH)->fs->omega;
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
        xB = A->get_cell(IBDSH,JBDSH)->pos.x;
        yB = A->get_cell(IBDSH,JBDSH)->pos.y;
        xC = A->get_cell(ICDSH,JCDSH)->pos.x;
        yC = A->get_cell(ICDSH,JCDSH)->pos.y;
	xD = A->get_ifj(i-1,j)->pos.x;
	yD = A->get_ifj(i-1,j)->pos.y;
	// These are the flow properties at the corners of the secondary cell.
	uA = A->get_ifj(i,j)->fs->vel.x;
        uB = A->get_cell(IBDSH,JBDSH)->fs->vel.x;
        uC = A->get_cell(ICDSH,JCDSH)->fs->vel.x;
	uD = A->get_ifj(i-1,j)->fs->vel.x;
	//
	vA = A->get_ifj(i,j)->fs->vel.y;
        vB = A->get_cell(IBDSH,JBDSH)->fs->vel.y;
        vC = A->get_cell(ICDSH,JCDSH)->fs->vel.y;
	vD = A->get_ifj(i-1,j)->fs->vel.y;
	//
	for ( int itm=0; itm<ntm; ++itm ) {
	    TA[itm] = A->get_ifj(i,j)->fs->gas->T[itm];
	    TB[itm] = A->get_cell(IBDSH,JBDSH)->fs->gas->T[itm];
	    TC[itm] = A->get_cell(ICDSH,JCDSH)->fs->gas->T[itm];
	    TD[itm] = A->get_ifj(i-1,j)->fs->gas->T[itm];
	}
	//
	if( get_diffusion_flag() == 1) { 
	    for( int isp = 0; isp < nsp; ++isp ) {
		fA[isp] = A->get_ifj(i,j)->fs->gas->massf[isp];
		fB[isp] = A->get_cell(IBDSH,JBDSH)->fs->gas->massf[isp];
		fC[isp] = A->get_cell(ICDSH,JCDSH)->fs->gas->massf[isp];
		fD[isp] = A->get_ifj(i-1,j)->fs->gas->massf[isp];
	    }
	}
	//
	APPLY_DIVERGENCE_THEOREM()
	//
	tkeA = A->get_ifj(i,j)->fs->tke;
        tkeB = A->get_cell(IBDSH,JBDSH)->fs->tke;
        tkeC = A->get_cell(ICDSH,JCDSH)->fs->tke;
	tkeD = A->get_ifj(i-1,j)->fs->tke;
	//
	omegaA = A->get_ifj(i,j)->fs->omega;
        omegaB = A->get_cell(IBDSH,JBDSH)->fs->omega;
        omegaC = A->get_cell(ICDSH,JCDSH)->fs->omega;
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
int viscous_derivatives_corners( Block *bdp )
{
    int i, j;
    FV_Vertex *vtx;
    FV_Interface *a, *b;
    FV_Cell *c;
    double fa, xa, ya, fb, xb, yb, fc, xc, yc, denom;
    int nsp = get_gas_model_ptr()->get_number_of_species();
    int ntm = get_gas_model_ptr()->get_number_of_modes();
    // North-East corner
    i = bdp->imax;
    j = bdp->jmax;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i+1,j+1);
    a = bdp->get_ifj(i,j+1);
    b = bdp->get_ifi(i+1,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos.x; yc = c->pos.y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( int itm=0; itm<ntm; ++itm ) {
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
    for( int isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    // South-East corner
    i = bdp->imax;
    j = bdp->jmin;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i+1,j);
    a = bdp->get_ifj(i,j);
    b = bdp->get_ifi(i+1,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos.x; yc = c->pos.y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( int itm=0; itm<ntm; ++itm ) {
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
    for( int isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    // South-West corner
    i = bdp->imin;
    j = bdp->jmin;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i,j);
    a = bdp->get_ifj(i,j);
    b = bdp->get_ifi(i,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos.x; yc = c->pos.y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( int itm=0; itm<ntm; ++itm ) {
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
    for( int isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    // North-West corner
    i = bdp->imin;
    j = bdp->jmax;
    c = bdp->get_cell(i,j);
    vtx = bdp->get_vtx(i,j+1);
    a = bdp->get_ifj(i,j+1);
    b = bdp->get_ifi(i,j);
    xa = a->pos.x; ya = a->pos.y;
    xb = b->pos.x; yb = b->pos.y;
    xc = c->pos.x; yc = c->pos.y;
    denom = (xc-xa)*(yb-ya) - (xb-xa)*(yc-ya);
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x;
    vtx->dudx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dudy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y;
    vtx->dvdx = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
    vtx->dvdy = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    for ( int itm=0; itm<ntm; ++itm ) {
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
    for( int isp = 0; isp < nsp; ++isp ) {
	fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp];
	vtx->dfdx[isp] = ((fc-fa)*(yb-ya) - (fb-fa)*(yc-ya))/denom;
	vtx->dfdy[isp] = ((fb-fa)*(xc-xa) - (fc-fa)*(xb-xa))/denom;
    }
    return SUCCESS;
} // end viscous_derivatives_corners()
