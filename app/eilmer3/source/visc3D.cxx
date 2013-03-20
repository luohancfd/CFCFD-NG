/** -*-m4-*- 
 * \file visc3D.cxx.m4
 * \ingroup eilmer3
 * \brief Functions to compute viscous fluxes for Eilmer3.
 *
 * \author Andrew Denman and (more recently) PJ
 * \version August 2004 bring code over from mb_cns.
 * \version November 2008 port from Elmer2 to Eilmer3 (via m4 macroprocessor), PJ
 * \version January 2010 more complete, with better edge and corner calculations, PJ
 *
 * Note that visc3D.cxx.m4 is the manually edited file.
 * The C++ file, visc3D.cxx, is machine generated.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "block.hh"
#include "kernel.hh"
#include "bc_defs.hh"
#include "bc.hh"
#include "visc3D.hh"
#include "diffusion.hh"
#include "baldwin_lomax.hh"
extern "C" {
#   include "../../../lib/util/source/useful.h"
}

// The following macro is used in the derivative functions.
// It applies macro $1 for each of the parameter sets.
// Within the set, first parameter is for cell center,
// second for cell interface and third is the result label.



// Working arrays for species derivatives
static std::vector<double> dfdx, dfdy, dfdz, jx, jy, jz;
// Working arrays for thermal derivatives
static std::vector<double> dTdx, dTdy, dTdz, qx, qy, qz, k_eff;

/** \brief Compute the viscous contribution to the cell interface fluxes.
 *
 * \param A      : pointer to the single-block data structure
 *
 * Compute the viscous contribution to the cell interface fluxes.
 * This contribution is added to the flux variables so make sure
 * that the inviscid (Euler) contributions have already been
 * computed and stored in the flux variables.
 */

int viscous_flux_3D(Block *A)
{
    FV_Vertex *Vtx1, *Vtx2, *Vtx3, *Vtx4;
    FV_Interface *IFace;
    double nx, ny, nz;
    size_t i, j, k;
    double dudx, dudy, dudz;
    double dvdx, dvdy, dvdz;
    double dwdx, dwdy, dwdz;
    double dtkedx, dtkedy, dtkedz;
    double domegadx, domegady, domegadz;
    double tau_xx, tau_yy, tau_xy;
    double tau_xz, tau_yz, tau_zz;
    double mu_eff, lmbda;
    double sigma = 0.5;
    double sigma_star = 0.6;
    double mu_effective;
    double tau_kx, tau_ky, tau_kz, tau_wx, tau_wy, tau_wz;

    double viscous_factor = get_viscous_factor();
    Gas_model *gmodel = get_gas_model_ptr();

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
	dTdz.resize(ntm);
	qx.resize(ntm);
	qy.resize(ntm);
	qz.resize(ntm);
	k_eff.resize(ntm);
    }

    // i-direction interfaces are treated as East interfaces.
    // t1 vector aligned with j-index direction
    // t2 vector aligned with k-index direction
    // The cycle [Vtx1,Vtx2,Vtx3,Vtx4] progresses counter-clockwise around 
    // the periphery of the face when the normal unit vector is pointing toward you.
    for (i = A->imin - 1; i <= A->imax; ++i) {
        for (j = A->jmin; j <= A->jmax; ++j) {
	    for (k = A->kmin; k <= A->kmax; ++k) {
		IFace = A->get_ifi(i+1,j,k);
		FlowState &fs = *(IFace->fs);
		Vtx1 = A->get_vtx(i+1,j,k);
		Vtx2 = A->get_vtx(i+1,j+1,k);
		Vtx3 = A->get_vtx(i+1,j+1,k+1);
		Vtx4 = A->get_vtx(i+1,j,k+1);
		// Determine some of the interface properties.
                if ( get_viscous_upwinding_flag() == 1 ) {
                    // Select one corner, based on the wind direction.
	            // When getting the velocity for upwinding, use the interface value
	            // unless we are at one of the block boundaries. 
	            // Use the interior cell value for boundary faces because we want to 
	            // know which direction is upwind, even for no-slip boundaries.
	            double vt1dp = 0.0;
	            double vt2dp = 0.0;
	            if ( i == A->imin-1 ) {
		        vt1dp = dot(A->get_cell(i+1,j,k)->fs->vel, IFace->t1);
		        vt2dp = dot(A->get_cell(i+1,j,k)->fs->vel, IFace->t2);
	            } else if ( i == A->imax ) {
		        vt1dp = dot(A->get_cell(i,j,k)->fs->vel, IFace->t1);
		        vt2dp = dot(A->get_cell(i,j,k)->fs->vel, IFace->t2);
	            } else {
		        vt1dp = dot(fs.vel,IFace->t1);
		        vt2dp = dot(fs.vel,IFace->t2);
	            }
		    		    		    if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudx = Vtx1->dudx; } else { dudx = Vtx4->dudx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudx = Vtx2->dudx; } else { dudx = Vtx3->dudx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudy = Vtx1->dudy; } else { dudy = Vtx4->dudy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudy = Vtx2->dudy; } else { dudy = Vtx3->dudy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudz = Vtx1->dudz; } else { dudz = Vtx4->dudz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudz = Vtx2->dudz; } else { dudz = Vtx3->dudz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdx = Vtx1->dvdx; } else { dvdx = Vtx4->dvdx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdx = Vtx2->dvdx; } else { dvdx = Vtx3->dvdx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdy = Vtx1->dvdy; } else { dvdy = Vtx4->dvdy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdy = Vtx2->dvdy; } else { dvdy = Vtx3->dvdy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdz = Vtx1->dvdz; } else { dvdz = Vtx4->dvdz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdz = Vtx2->dvdz; } else { dvdz = Vtx3->dvdz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdx = Vtx1->dwdx; } else { dwdx = Vtx4->dwdx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdx = Vtx2->dwdx; } else { dwdx = Vtx3->dwdx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdy = Vtx1->dwdy; } else { dwdy = Vtx4->dwdy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdy = Vtx2->dwdy; } else { dwdy = Vtx3->dwdy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdz = Vtx1->dwdz; } else { dwdz = Vtx4->dwdz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdz = Vtx2->dwdz; } else { dwdz = Vtx3->dwdz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedx = Vtx1->dtkedx; } else { dtkedx = Vtx4->dtkedx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedx = Vtx2->dtkedx; } else { dtkedx = Vtx3->dtkedx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedy = Vtx1->dtkedy; } else { dtkedy = Vtx4->dtkedy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedy = Vtx2->dtkedy; } else { dtkedy = Vtx3->dtkedy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedz = Vtx1->dtkedz; } else { dtkedz = Vtx4->dtkedz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedz = Vtx2->dtkedz; } else { dtkedz = Vtx3->dtkedz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegadx = Vtx1->domegadx; } else { domegadx = Vtx4->domegadx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegadx = Vtx2->domegadx; } else { domegadx = Vtx3->domegadx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegady = Vtx1->domegady; } else { domegady = Vtx4->domegady; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegady = Vtx2->domegady; } else { domegady = Vtx3->domegady; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegadz = Vtx1->domegadz; } else { domegadz = Vtx4->domegadz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegadz = Vtx2->domegadz; } else { domegadz = Vtx3->domegadz; }
                           }

		    for ( size_t itm=0; itm<ntm; ++itm ) {
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdx[itm] = Vtx1->dTdx[itm]; } else { dTdx[itm] = Vtx4->dTdx[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdx[itm] = Vtx2->dTdx[itm]; } else { dTdx[itm] = Vtx3->dTdx[itm]; }
                           }
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdy[itm] = Vtx1->dTdy[itm]; } else { dTdy[itm] = Vtx4->dTdy[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdy[itm] = Vtx2->dTdy[itm]; } else { dTdy[itm] = Vtx3->dTdy[itm]; }
                           }
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdz[itm] = Vtx1->dTdz[itm]; } else { dTdz[itm] = Vtx4->dTdz[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdz[itm] = Vtx2->dTdz[itm]; } else { dTdz[itm] = Vtx3->dTdz[itm]; }
                           }
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( size_t isp = 0; isp < nsp; ++isp ) {
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdx[isp] = Vtx1->dfdx[isp]; } else { dfdx[isp] = Vtx4->dfdx[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdx[isp] = Vtx2->dfdx[isp]; } else { dfdx[isp] = Vtx3->dfdx[isp]; }
                           }
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdy[isp] = Vtx1->dfdy[isp]; } else { dfdy[isp] = Vtx4->dfdy[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdy[isp] = Vtx2->dfdy[isp]; } else { dfdy[isp] = Vtx3->dfdy[isp]; }
                           }
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdz[isp] = Vtx1->dfdz[isp]; } else { dfdz[isp] = Vtx4->dfdz[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdz[isp] = Vtx2->dfdz[isp]; } else { dfdz[isp] = Vtx3->dfdz[isp]; }
                           }
		        }
                    }
                } else {
                    // Symmetric average.
		    		    		    dudx = 0.25*(Vtx1->dudx+Vtx2->dudx+Vtx3->dudx+Vtx4->dudx);
       dudy = 0.25*(Vtx1->dudy+Vtx2->dudy+Vtx3->dudy+Vtx4->dudy);
       dudz = 0.25*(Vtx1->dudz+Vtx2->dudz+Vtx3->dudz+Vtx4->dudz);
       dvdx = 0.25*(Vtx1->dvdx+Vtx2->dvdx+Vtx3->dvdx+Vtx4->dvdx);
       dvdy = 0.25*(Vtx1->dvdy+Vtx2->dvdy+Vtx3->dvdy+Vtx4->dvdy);
       dvdz = 0.25*(Vtx1->dvdz+Vtx2->dvdz+Vtx3->dvdz+Vtx4->dvdz);
       dwdx = 0.25*(Vtx1->dwdx+Vtx2->dwdx+Vtx3->dwdx+Vtx4->dwdx);
       dwdy = 0.25*(Vtx1->dwdy+Vtx2->dwdy+Vtx3->dwdy+Vtx4->dwdy);
       dwdz = 0.25*(Vtx1->dwdz+Vtx2->dwdz+Vtx3->dwdz+Vtx4->dwdz);
       dtkedx = 0.25*(Vtx1->dtkedx+Vtx2->dtkedx+Vtx3->dtkedx+Vtx4->dtkedx);
       dtkedy = 0.25*(Vtx1->dtkedy+Vtx2->dtkedy+Vtx3->dtkedy+Vtx4->dtkedy);
       dtkedz = 0.25*(Vtx1->dtkedz+Vtx2->dtkedz+Vtx3->dtkedz+Vtx4->dtkedz);
       domegadx = 0.25*(Vtx1->domegadx+Vtx2->domegadx+Vtx3->domegadx+Vtx4->domegadx);
       domegady = 0.25*(Vtx1->domegady+Vtx2->domegady+Vtx3->domegady+Vtx4->domegady);
       domegadz = 0.25*(Vtx1->domegadz+Vtx2->domegadz+Vtx3->domegadz+Vtx4->domegadz);

		    for ( size_t itm=0; itm<ntm; ++itm ) {
                        dTdx[itm] = 0.25*(Vtx1->dTdx[itm]+Vtx2->dTdx[itm]+Vtx3->dTdx[itm]+Vtx4->dTdx[itm]);
                        dTdy[itm] = 0.25*(Vtx1->dTdy[itm]+Vtx2->dTdy[itm]+Vtx3->dTdy[itm]+Vtx4->dTdy[itm]);
                        dTdz[itm] = 0.25*(Vtx1->dTdz[itm]+Vtx2->dTdz[itm]+Vtx3->dTdz[itm]+Vtx4->dTdz[itm]);
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( size_t isp = 0; isp < nsp; ++isp ) {
                            dfdx[isp] = 0.25*(Vtx1->dfdx[isp]+Vtx2->dfdx[isp]+Vtx3->dfdx[isp]+Vtx4->dfdx[isp]);
                            dfdy[isp] = 0.25*(Vtx1->dfdy[isp]+Vtx2->dfdy[isp]+Vtx3->dfdy[isp]+Vtx4->dfdy[isp]);
                            dfdz[isp] = 0.25*(Vtx1->dfdz[isp]+Vtx2->dfdz[isp]+Vtx3->dfdz[isp]+Vtx4->dfdz[isp]);
		        }
                    }
                }		    
                k_eff[0] = viscous_factor * (fs.gas->k[0] + fs.k_t);
		for ( size_t itm=1; itm<ntm; ++itm ) {
		    k_eff[itm] = viscous_factor * fs.gas->k[itm];
                }
		mu_eff =  viscous_factor * (fs.gas->mu + fs.mu_t);
		lmbda = -2.0/3.0 * mu_eff;
	        if( get_diffusion_flag() == 1 ) {
		    // Apply a diffusion model
		    double D_t = 0.0;
		    if ( get_k_omega_flag() == 1 ) {
                        double Sc_t = get_turbulence_schmidt_number();
                        D_t = fs.mu_t / (fs.gas->rho * Sc_t);
		    }
		    calculate_diffusion_fluxes(*(fs.gas),
					       D_t,
					       dfdx, dfdy, dfdz,
					       jx, jy, jz);
		    // NOTE: now applying viscous_factor to diffusive fluxes
		    for( size_t isp = 0; isp < nsp; ++isp ) {
		        jx[isp] *= viscous_factor;
		        jy[isp] *= viscous_factor;
		        jz[isp] *= viscous_factor;
		    }
	        } // end if get_diffusion_flag
		// 3-dimensional planar stresses.
		tau_xx = 2.0 * mu_eff * dudx + lmbda * (dudx + dvdy + dwdz);
		tau_yy = 2.0 * mu_eff * dvdy + lmbda * (dudx + dvdy + dwdz);
		tau_zz = 2.0 * mu_eff * dwdz + lmbda * (dudx + dvdy + dwdz);
		tau_xy = mu_eff * (dudy + dvdx);
		tau_xz = mu_eff * (dudz + dwdx);
		tau_yz = mu_eff * (dvdz + dwdy);
	        // Thermal conductivity
	        // NOTE: q[0] is total energy flux
	        qx[0] = k_eff[0] * dTdx[0];
	        qy[0] = k_eff[0] * dTdy[0];
	        qz[0] = k_eff[0] * dTdz[0];
	        for ( size_t itm=1; itm<ntm; ++itm ) {
		    qx[itm] = k_eff[itm] * dTdx[itm];
		    qy[itm] = k_eff[itm] * dTdy[itm];
		    qz[itm] = k_eff[itm] * dTdz[itm];
		    qx[0] += qx[itm];
		    qy[0] += qy[itm];
		    qz[0] += qz[itm];
	        }
	        if( get_diffusion_flag() == 1 ) {
		    for( size_t isp = 0; isp < nsp; ++isp ) {
		    	double h = gmodel->enthalpy(*(fs.gas), isp);
		        qx[0] -= jx[isp] * h;
		        qy[0] -= jy[isp] * h;
		        qz[0] -= jz[isp] * h;
		        for ( size_t itm=1; itm<ntm; ++itm ) {
                            double hmode = gmodel->modal_enthalpy(*(fs.gas), isp, itm);
			    qx[itm] -= jx[isp] * hmode;
			    qy[itm] -= jy[isp] * hmode;
			    qz[itm] -= jz[isp] * hmode;
		        }
		    }
	        }	    
	        if ( get_k_omega_flag() == 1 ) {
		    // Turbulence contribution to the shear stresses.
		    tau_xx -= 0.66667 * fs.gas->rho * fs.tke;
		    tau_yy -= 0.66667 * fs.gas->rho * fs.tke;
		    tau_zz -= 0.66667 * fs.gas->rho * fs.tke;
		    // Turbulence contribution to heat transfer.
		    mu_effective = fs.gas->mu + sigma_star * fs.mu_t;
		    qx[0] += mu_effective * dtkedx;
		    qy[0] += mu_effective * dtkedy;
		    qz[0] += mu_effective * dtkedz;
		    // Turbulence transport of the turbulence properties themselves.
		    tau_kx = mu_effective * dtkedx; 
		    tau_ky = mu_effective * dtkedy;
		    tau_kz = mu_effective * dtkedz;
		    mu_effective = fs.gas->mu + sigma * fs.mu_t;
		    tau_wx = mu_effective * domegadx; 
		    tau_wy = mu_effective * domegady; 
		    tau_wz = mu_effective * domegadz; 
	        } else {
		    tau_kx = 0.0;
		    tau_ky = 0.0;
		    tau_kz = 0.0;
		    tau_wx = 0.0;
		    tau_wy = 0.0;
		    tau_wz = 0.0;
	        }
		// Combine into fluxes: store as the dot product (F.n).
		ConservedQuantities &F = *(IFace->F);
		nx = IFace->n.x;
		ny = IFace->n.y;
		nz = IFace->n.z;
		// Mass flux -- NO CONTRIBUTION
		F.momentum.x -= tau_xx*nx + tau_xy*ny + tau_xz*nz;
		F.momentum.y -= tau_xy*nx + tau_yy*ny + tau_yz*nz;
		F.momentum.z -= tau_xz*nx + tau_yz*ny + tau_zz*nz;
		F.total_energy -=
		    (tau_xx*fs.vel.x + tau_xy*fs.vel.y + tau_xz*fs.vel.z + qx[0])*nx +
		    (tau_xy*fs.vel.x + tau_yy*fs.vel.y + tau_yz*fs.vel.z + qy[0])*ny +
		    (tau_xz*fs.vel.x + tau_yz*fs.vel.y + tau_zz*fs.vel.z + qz[0])*nz;
	        if ( get_k_omega_flag() == 1 ) {
		    F.tke -= tau_kx * nx + tau_ky * ny + tau_kz * nz;
		    F.omega -= tau_wx * nx + tau_wy * ny + tau_wz * nz;
	        }
                // Species mass flux
	        if( get_diffusion_flag() == 1 ) {
	  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		        F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
		    }
	        }
	        // Modal energy flux (skipping first mode as this is handled by total energy)
	        for ( size_t itm=1; itm<ntm; ++itm ) {
	    	    F.energies[itm] -= qx[itm]*nx + qy[itm]*ny + qz[itm]*nz;
	        }
	    } // k loop
	} // j loop
    } // i loop

    // j-interfaces are treated as North interfaces
    // t1 vector aligned with k-index direction
    // t2 vector aligned with i-index direction
    for (i = A->imin; i <= A->imax; ++i) {
        for (j = A->jmin - 1; j <= A->jmax; ++j) {
	    for (k = A->kmin; k <= A->kmax; ++k) {
		IFace = A->get_ifj(i,j+1,k);
		FlowState &fs = *(IFace->fs);
		Vtx1 = A->get_vtx(i,j+1,k);
		Vtx2 = A->get_vtx(i,j+1,k+1);
		Vtx3 = A->get_vtx(i+1,j+1,k+1);
		Vtx4 = A->get_vtx(i+1,j+1,k);
		// Determine some of the interface properties.
                if ( get_viscous_upwinding_flag() == 1 ) {
                    // Select one corner, based on the wind direction.
	            // When getting the velocity for upwinding, use the interface value
	            // unless we are at one of the block boundaries. 
	            // Use the interior cell value for boundary faces because we want to 
	            // know which direction is upwind, even for no-slip boundaries.
	            double vt1dp = 0.0;
	            double vt2dp = 0.0;
	            if ( j == A->jmin-1 ) {
		        vt1dp = dot(A->get_cell(i,j+1,k)->fs->vel, IFace->t1);
		        vt2dp = dot(A->get_cell(i,j+1,k)->fs->vel, IFace->t2);
	            } else if ( j == A->jmax ) {
		        vt1dp = dot(A->get_cell(i,j,k)->fs->vel, IFace->t1);
		        vt2dp = dot(A->get_cell(i,j,k)->fs->vel, IFace->t2);
	            } else {
		        vt1dp = dot(fs.vel,IFace->t1);
		        vt2dp = dot(fs.vel,IFace->t2);
	            }
		    if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudx = Vtx1->dudx; } else { dudx = Vtx4->dudx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudx = Vtx2->dudx; } else { dudx = Vtx3->dudx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudy = Vtx1->dudy; } else { dudy = Vtx4->dudy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudy = Vtx2->dudy; } else { dudy = Vtx3->dudy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudz = Vtx1->dudz; } else { dudz = Vtx4->dudz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudz = Vtx2->dudz; } else { dudz = Vtx3->dudz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdx = Vtx1->dvdx; } else { dvdx = Vtx4->dvdx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdx = Vtx2->dvdx; } else { dvdx = Vtx3->dvdx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdy = Vtx1->dvdy; } else { dvdy = Vtx4->dvdy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdy = Vtx2->dvdy; } else { dvdy = Vtx3->dvdy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdz = Vtx1->dvdz; } else { dvdz = Vtx4->dvdz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdz = Vtx2->dvdz; } else { dvdz = Vtx3->dvdz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdx = Vtx1->dwdx; } else { dwdx = Vtx4->dwdx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdx = Vtx2->dwdx; } else { dwdx = Vtx3->dwdx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdy = Vtx1->dwdy; } else { dwdy = Vtx4->dwdy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdy = Vtx2->dwdy; } else { dwdy = Vtx3->dwdy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdz = Vtx1->dwdz; } else { dwdz = Vtx4->dwdz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdz = Vtx2->dwdz; } else { dwdz = Vtx3->dwdz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedx = Vtx1->dtkedx; } else { dtkedx = Vtx4->dtkedx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedx = Vtx2->dtkedx; } else { dtkedx = Vtx3->dtkedx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedy = Vtx1->dtkedy; } else { dtkedy = Vtx4->dtkedy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedy = Vtx2->dtkedy; } else { dtkedy = Vtx3->dtkedy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedz = Vtx1->dtkedz; } else { dtkedz = Vtx4->dtkedz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedz = Vtx2->dtkedz; } else { dtkedz = Vtx3->dtkedz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegadx = Vtx1->domegadx; } else { domegadx = Vtx4->domegadx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegadx = Vtx2->domegadx; } else { domegadx = Vtx3->domegadx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegady = Vtx1->domegady; } else { domegady = Vtx4->domegady; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegady = Vtx2->domegady; } else { domegady = Vtx3->domegady; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegadz = Vtx1->domegadz; } else { domegadz = Vtx4->domegadz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegadz = Vtx2->domegadz; } else { domegadz = Vtx3->domegadz; }
                           }

		    for ( size_t itm=0; itm<ntm; ++itm ) {
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdx[itm] = Vtx1->dTdx[itm]; } else { dTdx[itm] = Vtx4->dTdx[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdx[itm] = Vtx2->dTdx[itm]; } else { dTdx[itm] = Vtx3->dTdx[itm]; }
                           }
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdy[itm] = Vtx1->dTdy[itm]; } else { dTdy[itm] = Vtx4->dTdy[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdy[itm] = Vtx2->dTdy[itm]; } else { dTdy[itm] = Vtx3->dTdy[itm]; }
                           }
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdz[itm] = Vtx1->dTdz[itm]; } else { dTdz[itm] = Vtx4->dTdz[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdz[itm] = Vtx2->dTdz[itm]; } else { dTdz[itm] = Vtx3->dTdz[itm]; }
                           }
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( size_t isp = 0; isp < nsp; ++isp ) {
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdx[isp] = Vtx1->dfdx[isp]; } else { dfdx[isp] = Vtx4->dfdx[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdx[isp] = Vtx2->dfdx[isp]; } else { dfdx[isp] = Vtx3->dfdx[isp]; }
                           }
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdy[isp] = Vtx1->dfdy[isp]; } else { dfdy[isp] = Vtx4->dfdy[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdy[isp] = Vtx2->dfdy[isp]; } else { dfdy[isp] = Vtx3->dfdy[isp]; }
                           }
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdz[isp] = Vtx1->dfdz[isp]; } else { dfdz[isp] = Vtx4->dfdz[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdz[isp] = Vtx2->dfdz[isp]; } else { dfdz[isp] = Vtx3->dfdz[isp]; }
                           }
		        }
                    }
                } else {
                    // Symmetric average.
 		    dudx = 0.25*(Vtx1->dudx+Vtx2->dudx+Vtx3->dudx+Vtx4->dudx);
       dudy = 0.25*(Vtx1->dudy+Vtx2->dudy+Vtx3->dudy+Vtx4->dudy);
       dudz = 0.25*(Vtx1->dudz+Vtx2->dudz+Vtx3->dudz+Vtx4->dudz);
       dvdx = 0.25*(Vtx1->dvdx+Vtx2->dvdx+Vtx3->dvdx+Vtx4->dvdx);
       dvdy = 0.25*(Vtx1->dvdy+Vtx2->dvdy+Vtx3->dvdy+Vtx4->dvdy);
       dvdz = 0.25*(Vtx1->dvdz+Vtx2->dvdz+Vtx3->dvdz+Vtx4->dvdz);
       dwdx = 0.25*(Vtx1->dwdx+Vtx2->dwdx+Vtx3->dwdx+Vtx4->dwdx);
       dwdy = 0.25*(Vtx1->dwdy+Vtx2->dwdy+Vtx3->dwdy+Vtx4->dwdy);
       dwdz = 0.25*(Vtx1->dwdz+Vtx2->dwdz+Vtx3->dwdz+Vtx4->dwdz);
       dtkedx = 0.25*(Vtx1->dtkedx+Vtx2->dtkedx+Vtx3->dtkedx+Vtx4->dtkedx);
       dtkedy = 0.25*(Vtx1->dtkedy+Vtx2->dtkedy+Vtx3->dtkedy+Vtx4->dtkedy);
       dtkedz = 0.25*(Vtx1->dtkedz+Vtx2->dtkedz+Vtx3->dtkedz+Vtx4->dtkedz);
       domegadx = 0.25*(Vtx1->domegadx+Vtx2->domegadx+Vtx3->domegadx+Vtx4->domegadx);
       domegady = 0.25*(Vtx1->domegady+Vtx2->domegady+Vtx3->domegady+Vtx4->domegady);
       domegadz = 0.25*(Vtx1->domegadz+Vtx2->domegadz+Vtx3->domegadz+Vtx4->domegadz);

		    for ( size_t itm=0; itm<ntm; ++itm ) {
                        dTdx[itm] = 0.25*(Vtx1->dTdx[itm]+Vtx2->dTdx[itm]+Vtx3->dTdx[itm]+Vtx4->dTdx[itm]);
                        dTdy[itm] = 0.25*(Vtx1->dTdy[itm]+Vtx2->dTdy[itm]+Vtx3->dTdy[itm]+Vtx4->dTdy[itm]);
                        dTdz[itm] = 0.25*(Vtx1->dTdz[itm]+Vtx2->dTdz[itm]+Vtx3->dTdz[itm]+Vtx4->dTdz[itm]);
                    }
 	            if( get_diffusion_flag() == 1 ) {
                        // derivatives needed for diffusion model, below
		        for( size_t isp = 0; isp < nsp; ++isp ) {
                            dfdx[isp] = 0.25*(Vtx1->dfdx[isp]+Vtx2->dfdx[isp]+Vtx3->dfdx[isp]+Vtx4->dfdx[isp]);
                            dfdy[isp] = 0.25*(Vtx1->dfdy[isp]+Vtx2->dfdy[isp]+Vtx3->dfdy[isp]+Vtx4->dfdy[isp]);
                            dfdz[isp] = 0.25*(Vtx1->dfdz[isp]+Vtx2->dfdz[isp]+Vtx3->dfdz[isp]+Vtx4->dfdz[isp]);
		        }
                    }
                }
                k_eff[0] = viscous_factor * (fs.gas->k[0] + fs.k_t);
		for ( size_t itm=1; itm<ntm; ++itm ) {
		    k_eff[itm] = viscous_factor * fs.gas->k[itm];
                }
		mu_eff =  viscous_factor * (fs.gas->mu + fs.mu_t);
		lmbda = -2.0/3.0 * mu_eff;
	        if( get_diffusion_flag() == 1 ) {
		    // Apply a diffusion model
		    double D_t = 0.0;
		    if ( get_k_omega_flag() == 1 ) {
                        double Sc_t = get_turbulence_schmidt_number();
                        D_t = fs.mu_t / (fs.gas->rho * Sc_t);
		    }
		    calculate_diffusion_fluxes(*(fs.gas),
					       D_t,
					       dfdx, dfdy, dfdz,
					       jx, jy, jz);
		    // NOTE: now applying viscous_factor to diffusive fluxes
		    for( size_t isp = 0; isp < nsp; ++isp ) {
		        jx[isp] *= viscous_factor;
		        jy[isp] *= viscous_factor;
		        jz[isp] *= viscous_factor;
		    }
	        } // end if get_diffusion_flag
		// 3-dimensional planar stresses.
		tau_xx = 2.0*mu_eff*dudx + lmbda*(dudx + dvdy + dwdz);
		tau_yy = 2.0*mu_eff*dvdy + lmbda*(dudx + dvdy + dwdz);
		tau_zz = 2.0*mu_eff*dwdz + lmbda*(dudx + dvdy + dwdz);
		tau_xy = mu_eff * (dudy + dvdx);
		tau_xz = mu_eff * (dudz + dwdx);
		tau_yz = mu_eff * (dvdz + dwdy);
	        // Thermal conductivity
	        // NOTE: q[0] is total energy flux
	        qx[0] = k_eff[0] * dTdx[0];
	        qy[0] = k_eff[0] * dTdy[0];
	        qz[0] = k_eff[0] * dTdz[0];
	        for ( size_t itm=1; itm<ntm; ++itm ) {
		    qx[itm] = k_eff[itm] * dTdx[itm];
		    qy[itm] = k_eff[itm] * dTdy[itm];
		    qz[itm] = k_eff[itm] * dTdz[itm];
		    qx[0] += qx[itm];
		    qy[0] += qy[itm];
		    qz[0] += qz[itm];
	        }
	        if( get_diffusion_flag() == 1 ) {
		    for( size_t isp = 0; isp < nsp; ++isp ) {
		    	double h = gmodel->enthalpy(*(fs.gas), isp);
		        qx[0] -= jx[isp] * h;
		        qy[0] -= jy[isp] * h;
		        qz[0] -= jz[isp] * h;
		        for ( size_t itm=1; itm<ntm; ++itm ) {
			    double hmode = gmodel->modal_enthalpy(*(fs.gas), isp, itm);
			    qx[itm] -= jx[isp] * hmode;
			    qy[itm] -= jy[isp] * hmode;
			    qz[itm] -= jz[isp] * hmode;
		        }
		    }
	        }	    
	        if ( get_k_omega_flag() == 1 ) {
		    // Turbulence contribution to the shear stresses.
		    tau_xx -= 0.66667 * fs.gas->rho * fs.tke;
		    tau_yy -= 0.66667 * fs.gas->rho * fs.tke;
		    tau_zz -= 0.66667 * fs.gas->rho * fs.tke;
		    // Turbulence contribution to heat transfer.
		    mu_effective = fs.gas->mu + sigma_star * fs.mu_t;
		    qx[0] += mu_effective * dtkedx;
		    qy[0] += mu_effective * dtkedy;
		    qz[0] += mu_effective * dtkedz;
		    // Turbulence transport of the turbulence properties themselves.
		    tau_kx = mu_effective * dtkedx; 
		    tau_ky = mu_effective * dtkedy;
		    tau_kz = mu_effective * dtkedz;
		    mu_effective = fs.gas->mu + sigma * fs.mu_t;
		    tau_wx = mu_effective * domegadx; 
		    tau_wy = mu_effective * domegady; 
		    tau_wz = mu_effective * domegadz; 
	        } else {
		    tau_kx = 0.0;
		    tau_ky = 0.0;
		    tau_kz = 0.0;
		    tau_wx = 0.0;
		    tau_wy = 0.0;
		    tau_wz = 0.0;
	        }
		// Combine into fluxes: store as the dot product (F.n).
		ConservedQuantities &F = *(IFace->F);
		nx = IFace->n.x;
		ny = IFace->n.y;
		nz = IFace->n.z;
		// Mass flux -- NO CONTRIBUTION
		F.momentum.x -= tau_xx*nx + tau_xy*ny + tau_xz*nz;
		F.momentum.y -= tau_xy*nx + tau_yy*ny + tau_yz*nz;
		F.momentum.z -= tau_xz*nx + tau_yz*ny + tau_zz*nz;
		F.total_energy -=
		    (tau_xx*fs.vel.x + tau_xy*fs.vel.y + tau_xz*fs.vel.z + qx[0])*nx +
		    (tau_xy*fs.vel.x + tau_yy*fs.vel.y + tau_yz*fs.vel.z + qy[0])*ny +
		    (tau_xz*fs.vel.x + tau_yz*fs.vel.y + tau_zz*fs.vel.z + qz[0])*nz;
	        if ( get_k_omega_flag() == 1 ) {
		    F.tke -= tau_kx * nx + tau_ky * ny + tau_kz * nz;
		    F.omega -= tau_wx * nx + tau_wy * ny + tau_wz * nz;
	        }
                // Species mass flux
	        if( get_diffusion_flag() == 1 ) {
	  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		        F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
		    }
	        }
	        // Modal energy flux (skipping first mode as this is handled by total energy)
	        for ( size_t itm=1; itm<ntm; ++itm ) {
	    	    F.energies[itm] -= qx[itm]*nx + qy[itm]*ny + qz[itm]*nz;
	        }
	    } // k loop
	} // j loop
    } // i loop

    // k-interfaces are treated as Top interfaces
    // t1 vector aligned with i-index direction
    // t2 vector aligned with j-index direction
    for (i = A->imin; i <= A->imax; ++i) {
        for (j = A->jmin; j <= A->jmax; ++j) {
	    for (k = A->kmin - 1; k <= A->kmax; ++k) {
		IFace = A->get_ifk(i,j,k+1);
		FlowState &fs = *(IFace->fs);
		Vtx1 = A->get_vtx(i,j,k+1);
		Vtx2 = A->get_vtx(i+1,j,k+1);
		Vtx3 = A->get_vtx(i+1,j+1,k+1);
		Vtx4 = A->get_vtx(i,j+1,k+1);
		// Determine some of the interface properties.
                if ( get_viscous_upwinding_flag() == 1 ) {
                    // Select one corner, based on the wind direction.
	            // When getting the velocity for upwinding, use the interface value
	            // unless we are at one of the block boundaries. 
	            // Use the interior cell value for boundary faces because we want to 
	            // know which direction is upwind, even for no-slip boundaries.
	            double vt1dp = 0.0;
	            double vt2dp = 0.0;
	            if ( k == A->kmin-1 ) {
		        vt1dp = dot(A->get_cell(i,j,k+1)->fs->vel, IFace->t1);
		        vt2dp = dot(A->get_cell(i,j,k+1)->fs->vel, IFace->t2);
	            } else if ( k == A->kmax ) {
		        vt1dp = dot(A->get_cell(i,j,k)->fs->vel, IFace->t1);
		        vt2dp = dot(A->get_cell(i,j,k)->fs->vel, IFace->t2);
	            } else {
		        vt1dp = dot(fs.vel,IFace->t1);
		        vt2dp = dot(fs.vel,IFace->t2);
	            }
		    if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudx = Vtx1->dudx; } else { dudx = Vtx4->dudx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudx = Vtx2->dudx; } else { dudx = Vtx3->dudx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudy = Vtx1->dudy; } else { dudy = Vtx4->dudy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudy = Vtx2->dudy; } else { dudy = Vtx3->dudy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dudz = Vtx1->dudz; } else { dudz = Vtx4->dudz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dudz = Vtx2->dudz; } else { dudz = Vtx3->dudz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdx = Vtx1->dvdx; } else { dvdx = Vtx4->dvdx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdx = Vtx2->dvdx; } else { dvdx = Vtx3->dvdx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdy = Vtx1->dvdy; } else { dvdy = Vtx4->dvdy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdy = Vtx2->dvdy; } else { dvdy = Vtx3->dvdy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dvdz = Vtx1->dvdz; } else { dvdz = Vtx4->dvdz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dvdz = Vtx2->dvdz; } else { dvdz = Vtx3->dvdz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdx = Vtx1->dwdx; } else { dwdx = Vtx4->dwdx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdx = Vtx2->dwdx; } else { dwdx = Vtx3->dwdx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdy = Vtx1->dwdy; } else { dwdy = Vtx4->dwdy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdy = Vtx2->dwdy; } else { dwdy = Vtx3->dwdy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dwdz = Vtx1->dwdz; } else { dwdz = Vtx4->dwdz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dwdz = Vtx2->dwdz; } else { dwdz = Vtx3->dwdz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedx = Vtx1->dtkedx; } else { dtkedx = Vtx4->dtkedx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedx = Vtx2->dtkedx; } else { dtkedx = Vtx3->dtkedx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedy = Vtx1->dtkedy; } else { dtkedy = Vtx4->dtkedy; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedy = Vtx2->dtkedy; } else { dtkedy = Vtx3->dtkedy; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dtkedz = Vtx1->dtkedz; } else { dtkedz = Vtx4->dtkedz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dtkedz = Vtx2->dtkedz; } else { dtkedz = Vtx3->dtkedz; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegadx = Vtx1->domegadx; } else { domegadx = Vtx4->domegadx; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegadx = Vtx2->domegadx; } else { domegadx = Vtx3->domegadx; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegady = Vtx1->domegady; } else { domegady = Vtx4->domegady; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegady = Vtx2->domegady; } else { domegady = Vtx3->domegady; }
                           }
       if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { domegadz = Vtx1->domegadz; } else { domegadz = Vtx4->domegadz; }
                           } else {
		               if ( vt2dp >= 0.0 ) { domegadz = Vtx2->domegadz; } else { domegadz = Vtx3->domegadz; }
                           }

		    for ( size_t itm=0; itm<ntm; ++itm ) {
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdx[itm] = Vtx1->dTdx[itm]; } else { dTdx[itm] = Vtx4->dTdx[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdx[itm] = Vtx2->dTdx[itm]; } else { dTdx[itm] = Vtx3->dTdx[itm]; }
                           }
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdy[itm] = Vtx1->dTdy[itm]; } else { dTdy[itm] = Vtx4->dTdy[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdy[itm] = Vtx2->dTdy[itm]; } else { dTdy[itm] = Vtx3->dTdy[itm]; }
                           }
                        if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dTdz[itm] = Vtx1->dTdz[itm]; } else { dTdz[itm] = Vtx4->dTdz[itm]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dTdz[itm] = Vtx2->dTdz[itm]; } else { dTdz[itm] = Vtx3->dTdz[itm]; }
                           }
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( size_t isp = 0; isp < nsp; ++isp ) {
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdx[isp] = Vtx1->dfdx[isp]; } else { dfdx[isp] = Vtx4->dfdx[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdx[isp] = Vtx2->dfdx[isp]; } else { dfdx[isp] = Vtx3->dfdx[isp]; }
                           }
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdy[isp] = Vtx1->dfdy[isp]; } else { dfdy[isp] = Vtx4->dfdy[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdy[isp] = Vtx2->dfdy[isp]; } else { dfdy[isp] = Vtx3->dfdy[isp]; }
                           }
                            if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { dfdz[isp] = Vtx1->dfdz[isp]; } else { dfdz[isp] = Vtx4->dfdz[isp]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { dfdz[isp] = Vtx2->dfdz[isp]; } else { dfdz[isp] = Vtx3->dfdz[isp]; }
                           }
		        }
                    }
                } else {
                    // Symmetric average.
		    dudx = 0.25*(Vtx1->dudx+Vtx2->dudx+Vtx3->dudx+Vtx4->dudx);
       dudy = 0.25*(Vtx1->dudy+Vtx2->dudy+Vtx3->dudy+Vtx4->dudy);
       dudz = 0.25*(Vtx1->dudz+Vtx2->dudz+Vtx3->dudz+Vtx4->dudz);
       dvdx = 0.25*(Vtx1->dvdx+Vtx2->dvdx+Vtx3->dvdx+Vtx4->dvdx);
       dvdy = 0.25*(Vtx1->dvdy+Vtx2->dvdy+Vtx3->dvdy+Vtx4->dvdy);
       dvdz = 0.25*(Vtx1->dvdz+Vtx2->dvdz+Vtx3->dvdz+Vtx4->dvdz);
       dwdx = 0.25*(Vtx1->dwdx+Vtx2->dwdx+Vtx3->dwdx+Vtx4->dwdx);
       dwdy = 0.25*(Vtx1->dwdy+Vtx2->dwdy+Vtx3->dwdy+Vtx4->dwdy);
       dwdz = 0.25*(Vtx1->dwdz+Vtx2->dwdz+Vtx3->dwdz+Vtx4->dwdz);
       dtkedx = 0.25*(Vtx1->dtkedx+Vtx2->dtkedx+Vtx3->dtkedx+Vtx4->dtkedx);
       dtkedy = 0.25*(Vtx1->dtkedy+Vtx2->dtkedy+Vtx3->dtkedy+Vtx4->dtkedy);
       dtkedz = 0.25*(Vtx1->dtkedz+Vtx2->dtkedz+Vtx3->dtkedz+Vtx4->dtkedz);
       domegadx = 0.25*(Vtx1->domegadx+Vtx2->domegadx+Vtx3->domegadx+Vtx4->domegadx);
       domegady = 0.25*(Vtx1->domegady+Vtx2->domegady+Vtx3->domegady+Vtx4->domegady);
       domegadz = 0.25*(Vtx1->domegadz+Vtx2->domegadz+Vtx3->domegadz+Vtx4->domegadz);

		    for ( size_t itm=0; itm<ntm; ++itm ) {
                        dTdx[itm] = 0.25*(Vtx1->dTdx[itm]+Vtx2->dTdx[itm]+Vtx3->dTdx[itm]+Vtx4->dTdx[itm]);
                        dTdy[itm] = 0.25*(Vtx1->dTdy[itm]+Vtx2->dTdy[itm]+Vtx3->dTdy[itm]+Vtx4->dTdy[itm]);
                        dTdz[itm] = 0.25*(Vtx1->dTdz[itm]+Vtx2->dTdz[itm]+Vtx3->dTdz[itm]+Vtx4->dTdz[itm]);
                    }
 	            if( get_diffusion_flag() == 1 ) {
                        // derivatives needed for diffusion model, below
		        for( size_t isp = 0; isp < nsp; ++isp ) {
                            dfdx[isp] = 0.25*(Vtx1->dfdx[isp]+Vtx2->dfdx[isp]+Vtx3->dfdx[isp]+Vtx4->dfdx[isp]);
                            dfdy[isp] = 0.25*(Vtx1->dfdy[isp]+Vtx2->dfdy[isp]+Vtx3->dfdy[isp]+Vtx4->dfdy[isp]);
                            dfdz[isp] = 0.25*(Vtx1->dfdz[isp]+Vtx2->dfdz[isp]+Vtx3->dfdz[isp]+Vtx4->dfdz[isp]);
		        }
                    }
                }
                k_eff[0] = viscous_factor * (fs.gas->k[0] + fs.k_t);
		for ( size_t itm=1; itm<ntm; ++itm ) {
		    k_eff[itm] = viscous_factor * fs.gas->k[itm];
                }
		mu_eff =  viscous_factor * (fs.gas->mu + fs.mu_t);
		lmbda = -2.0/3.0 * mu_eff;
 	        if( get_diffusion_flag() == 1 ) {
		    // Apply a diffusion model
		    double D_t = 0.0;
		    if ( get_k_omega_flag() == 1 ) {
                        double Sc_t = get_turbulence_schmidt_number();
                        D_t = fs.mu_t / (fs.gas->rho * Sc_t);
		    }
		    calculate_diffusion_fluxes(*(fs.gas),
					       D_t,
					       dfdx, dfdy, dfdz,
					       jx, jy, jz);
		    // NOTE: now applying viscous_factor to diffusive fluxes
		    for( size_t isp = 0; isp < nsp; ++isp ) {
		        jx[isp] *= viscous_factor;
		        jy[isp] *= viscous_factor;
		        jz[isp] *= viscous_factor;
		    }
	        } // end if get_diffusion_flag
		// 3-dimensional planar stresses.
		tau_xx = 2.0*mu_eff*dudx + lmbda*(dudx + dvdy + dwdz);
		tau_yy = 2.0*mu_eff*dvdy + lmbda*(dudx + dvdy + dwdz);
		tau_zz = 2.0*mu_eff*dwdz + lmbda*(dudx + dvdy + dwdz);
		tau_xy = mu_eff * (dudy + dvdx);
		tau_xz = mu_eff * (dudz + dwdx);
		tau_yz = mu_eff * (dvdz + dwdy);
	        // Thermal conductivity
	        // NOTE: q[0] is total energy flux
	        qx[0] = k_eff[0] * dTdx[0];
	        qy[0] = k_eff[0] * dTdy[0];
	        qz[0] = k_eff[0] * dTdz[0];
	        for ( size_t itm=1; itm<ntm; ++itm ) {
		    qx[itm] = k_eff[itm] * dTdx[itm];
		    qy[itm] = k_eff[itm] * dTdy[itm];
		    qz[itm] = k_eff[itm] * dTdz[itm];
		    qx[0] += qx[itm];
		    qy[0] += qy[itm];
		    qz[0] += qz[itm];
	        }
	        if( get_diffusion_flag() == 1 ) {
		    for( size_t isp = 0; isp < nsp; ++isp ) {
		    	double h = gmodel->enthalpy(*(fs.gas), isp);
		        qx[0] -= jx[isp] * h;
		        qy[0] -= jy[isp] * h;
		        qz[0] -= jz[isp] * h;
		        for ( size_t itm=1; itm<ntm; ++itm ) {
                            double hmode = gmodel->modal_enthalpy(*(fs.gas), isp, itm);
			    qx[itm] -= jx[isp] * hmode;
			    qy[itm] -= jy[isp] * hmode;
			    qz[itm] -= jz[isp] * hmode;
		        }
		    }
	        }	    
	        if ( get_k_omega_flag() == 1 ) {
		    // Turbulence contribution to the shear stresses.
		    tau_xx -= 0.66667 * fs.gas->rho * fs.tke;
		    tau_yy -= 0.66667 * fs.gas->rho * fs.tke;
		    tau_zz -= 0.66667 * fs.gas->rho * fs.tke;
		    // Turbulence contribution to heat transfer.
		    mu_effective = fs.gas->mu + sigma_star * fs.mu_t;
		    qx[0] += mu_effective * dtkedx;
		    qy[0] += mu_effective * dtkedy;
		    qz[0] += mu_effective * dtkedz;
		    // Turbulence transport of the turbulence properties themselves.
		    tau_kx = mu_effective * dtkedx; 
		    tau_ky = mu_effective * dtkedy;
		    tau_kz = mu_effective * dtkedz;
		    mu_effective = fs.gas->mu + sigma * fs.mu_t;
		    tau_wx = mu_effective * domegadx; 
		    tau_wy = mu_effective * domegady; 
		    tau_wz = mu_effective * domegadz; 
	        } else {
		    tau_kx = 0.0;
		    tau_ky = 0.0;
		    tau_kz = 0.0;
		    tau_wx = 0.0;
		    tau_wy = 0.0;
		    tau_wz = 0.0;
	        }
		// Combine into fluxes: store as the dot product (F.n).
		ConservedQuantities &F = *(IFace->F);
		nx = IFace->n.x;
		ny = IFace->n.y;
		nz = IFace->n.z;
		// Mass flux -- NO CONTRIBUTION
		F.momentum.x -= tau_xx*nx + tau_xy*ny + tau_xz*nz;
		F.momentum.y -= tau_xy*nx + tau_yy*ny + tau_yz*nz;
		F.momentum.z -= tau_xz*nx + tau_yz*ny + tau_zz*nz;
		F.total_energy -=
		    (tau_xx*fs.vel.x + tau_xy*fs.vel.y + tau_xz*fs.vel.z + qx[0])*nx +
		    (tau_xy*fs.vel.x + tau_yy*fs.vel.y + tau_yz*fs.vel.z + qy[0])*ny +
		    (tau_xz*fs.vel.x + tau_yz*fs.vel.y + tau_zz*fs.vel.z + qz[0])*nz;
	        if ( get_k_omega_flag() == 1 ) {
		    F.tke -= tau_kx * nx + tau_ky * ny + tau_kz * nz;
		    F.omega -= tau_wx * nx + tau_wy * ny + tau_wz * nz;
	        }
                // Species mass flux
	        if( get_diffusion_flag() == 1 ) {
	  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		        F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
		    }
	        }
	        // Modal energy flux (skipping first mode as this is handled by total energy)
	        for ( size_t itm=1; itm<ntm; ++itm ) {
	    	    F.energies[itm] -= qx[itm]*nx + qy[itm]*ny + qz[itm]*nz;
	        }
	    } // k loop
	} // j loop
    } // i loop
    return SUCCESS;
} // end viscous_flux_3D()

//-------------------------------------------------------------------------------------

/** \brief Compute the derivatives needed for viscous terms.
 *
 * These derivative values are associated with the primary cell vertices.
 * and are obtained by applying the divergence theorem, integrating over
 * the surfaces of the secondary cells.
 * The secondary cells are centred on the vertices of the primary cells
 * and have primary cell centres as their corners.
 */
int viscous_derivatives_3D(Block *A)
{
    size_t i, j, k;
    double q_e, q_w, q_n, q_s, q_top, q_bottom;
    double vol_inv;
    FV_Vertex *sec_ctr;
    // The ABC... notation is from Ian Johnston, and was used by Andrew Denman.
    // We'll keep using it to minimize the change from Elmer2 to Elmer3.
    FV_Cell *cA, *cB, *cC, *cD, *cE, *cF, *cG, *cH;
    FV_Interface *fA, *fB, *fC, *fD, *fE, *fF, *fG, *fH;
    FV_Interface *north, *east, *south, *west, *bottom, *top;

    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t ntm = gmodel->get_number_of_modes();

    // First, do all of the internal secondary cells.
    // i.e. Those not on a boundary.
    for (i = A->imin; i <= A->imax - 1; ++i) {
        for (j = A->jmin; j <= A->jmax - 1; ++j) {
            for (k = A->kmin; k <= A->kmax - 1; ++k) {
		sec_ctr = A->get_vtx(i+1,j+1,k+1);
		vol_inv = 1.0 / sec_ctr->volume;
		// These are the corners of the secondary cell.
		cA = A->get_cell(i,j,k+1);
		cB = A->get_cell(i+1,j,k+1);
		cC = A->get_cell(i+1,j,k);
		cD = A->get_cell(i,j,k);
		cE = A->get_cell(i,j+1,k+1);
		cF = A->get_cell(i+1,j+1,k+1);
		cG = A->get_cell(i+1,j+1,k);
		cH = A->get_cell(i,j+1,k);
		
		// Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);
                		
		// Average property value on each face.
		       q_e = 0.25 * (cB->fs->vel.x + cC->fs->vel.x + cF->fs->vel.x + cG->fs->vel.x);
		       q_w = 0.25 * (cA->fs->vel.x + cD->fs->vel.x + cH->fs->vel.x + cE->fs->vel.x);
		       q_n = 0.25 * (cE->fs->vel.x + cF->fs->vel.x + cG->fs->vel.x + cH->fs->vel.x);
		       q_s = 0.25 * (cA->fs->vel.x + cB->fs->vel.x + cC->fs->vel.x + cD->fs->vel.x);
		       q_top = 0.25 * (cA->fs->vel.x + cB->fs->vel.x + cE->fs->vel.x + cF->fs->vel.x);
		       q_bottom = 0.25 * (cD->fs->vel.x + cC->fs->vel.x + cG->fs->vel.x + cH->fs->vel.x);
		       // Apply the divergence theorem.
		       sec_ctr->dudx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
		       sec_ctr->dudy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
		       sec_ctr->dudz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
		       q_e = 0.25 * (cB->fs->vel.y + cC->fs->vel.y + cF->fs->vel.y + cG->fs->vel.y);
		       q_w = 0.25 * (cA->fs->vel.y + cD->fs->vel.y + cH->fs->vel.y + cE->fs->vel.y);
		       q_n = 0.25 * (cE->fs->vel.y + cF->fs->vel.y + cG->fs->vel.y + cH->fs->vel.y);
		       q_s = 0.25 * (cA->fs->vel.y + cB->fs->vel.y + cC->fs->vel.y + cD->fs->vel.y);
		       q_top = 0.25 * (cA->fs->vel.y + cB->fs->vel.y + cE->fs->vel.y + cF->fs->vel.y);
		       q_bottom = 0.25 * (cD->fs->vel.y + cC->fs->vel.y + cG->fs->vel.y + cH->fs->vel.y);
		       // Apply the divergence theorem.
		       sec_ctr->dvdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
		       sec_ctr->dvdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
		       sec_ctr->dvdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
		       q_e = 0.25 * (cB->fs->vel.z + cC->fs->vel.z + cF->fs->vel.z + cG->fs->vel.z);
		       q_w = 0.25 * (cA->fs->vel.z + cD->fs->vel.z + cH->fs->vel.z + cE->fs->vel.z);
		       q_n = 0.25 * (cE->fs->vel.z + cF->fs->vel.z + cG->fs->vel.z + cH->fs->vel.z);
		       q_s = 0.25 * (cA->fs->vel.z + cB->fs->vel.z + cC->fs->vel.z + cD->fs->vel.z);
		       q_top = 0.25 * (cA->fs->vel.z + cB->fs->vel.z + cE->fs->vel.z + cF->fs->vel.z);
		       q_bottom = 0.25 * (cD->fs->vel.z + cC->fs->vel.z + cG->fs->vel.z + cH->fs->vel.z);
		       // Apply the divergence theorem.
		       sec_ctr->dwdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
		       sec_ctr->dwdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
		       sec_ctr->dwdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
		       q_e = 0.25 * (cB->fs->tke + cC->fs->tke + cF->fs->tke + cG->fs->tke);
		       q_w = 0.25 * (cA->fs->tke + cD->fs->tke + cH->fs->tke + cE->fs->tke);
		       q_n = 0.25 * (cE->fs->tke + cF->fs->tke + cG->fs->tke + cH->fs->tke);
		       q_s = 0.25 * (cA->fs->tke + cB->fs->tke + cC->fs->tke + cD->fs->tke);
		       q_top = 0.25 * (cA->fs->tke + cB->fs->tke + cE->fs->tke + cF->fs->tke);
		       q_bottom = 0.25 * (cD->fs->tke + cC->fs->tke + cG->fs->tke + cH->fs->tke);
		       // Apply the divergence theorem.
		       sec_ctr->dtkedx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
		       sec_ctr->dtkedy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
		       sec_ctr->dtkedz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
		       q_e = 0.25 * (cB->fs->omega + cC->fs->omega + cF->fs->omega + cG->fs->omega);
		       q_w = 0.25 * (cA->fs->omega + cD->fs->omega + cH->fs->omega + cE->fs->omega);
		       q_n = 0.25 * (cE->fs->omega + cF->fs->omega + cG->fs->omega + cH->fs->omega);
		       q_s = 0.25 * (cA->fs->omega + cB->fs->omega + cC->fs->omega + cD->fs->omega);
		       q_top = 0.25 * (cA->fs->omega + cB->fs->omega + cE->fs->omega + cF->fs->omega);
		       q_bottom = 0.25 * (cD->fs->omega + cC->fs->omega + cG->fs->omega + cH->fs->omega);
		       // Apply the divergence theorem.
		       sec_ctr->domegadx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
		       sec_ctr->domegady = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
		       sec_ctr->domegadz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);

		
	        // Apply the divergence theorem to the primary temperature derivatives only.
                for ( size_t itm=0; itm<ntm; ++itm ) {
		    // Average property value on each face.
		       q_e = 0.25 * (cB->fs->gas->T[itm] + cC->fs->gas->T[itm] + cF->fs->gas->T[itm] + cG->fs->gas->T[itm]);
		       q_w = 0.25 * (cA->fs->gas->T[itm] + cD->fs->gas->T[itm] + cH->fs->gas->T[itm] + cE->fs->gas->T[itm]);
		       q_n = 0.25 * (cE->fs->gas->T[itm] + cF->fs->gas->T[itm] + cG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
		       q_s = 0.25 * (cA->fs->gas->T[itm] + cB->fs->gas->T[itm] + cC->fs->gas->T[itm] + cD->fs->gas->T[itm]);
		       q_top = 0.25 * (cA->fs->gas->T[itm] + cB->fs->gas->T[itm] + cE->fs->gas->T[itm] + cF->fs->gas->T[itm]);
		       q_bottom = 0.25 * (cD->fs->gas->T[itm] + cC->fs->gas->T[itm] + cG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
		       // Apply the divergence theorem.
		       sec_ctr->dTdx[itm] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
		       sec_ctr->dTdy[itm] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
		       sec_ctr->dTdz[itm] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
                }
  	        for( size_t isp = 0; isp < nsp; ++isp ) {
		    // Average property value on each face.
		       q_e = 0.25 * (cB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + cG->fs->gas->massf[isp]);
		       q_w = 0.25 * (cA->fs->gas->massf[isp] + cD->fs->gas->massf[isp] + cH->fs->gas->massf[isp] + cE->fs->gas->massf[isp]);
		       q_n = 0.25 * (cE->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
		       q_s = 0.25 * (cA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cD->fs->gas->massf[isp]);
		       q_top = 0.25 * (cA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + cE->fs->gas->massf[isp] + cF->fs->gas->massf[isp]);
		       q_bottom = 0.25 * (cD->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
		       // Apply the divergence theorem.
		       sec_ctr->dfdx[isp] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
		       sec_ctr->dfdy[isp] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
		       sec_ctr->dfdz[isp] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
                }
	    } // k loop
        } // j loop
    } // i loop

    // Now, do the boundaries as half cells.
    // East boundary surface.
    i = A->imax;
    for (j = A->jmin; j <= A->jmax - 1; ++j) {
	for (k = A->kmin; k <= A->kmax - 1; ++k) {
	    sec_ctr = A->get_vtx(i+1,j+1,k+1);
	    vol_inv = 1.0 / sec_ctr->volume;
	    // The following are the corners of the secondary cell.
	    cA = A->get_cell(i,j,k+1);
	    fB = A->get_ifi(i+1,j,k+1);
	    fC = A->get_ifi(i+1,j,k);
	    cD = A->get_cell(i,j,k);
	    cE = A->get_cell(i,j+1,k+1);
	    fF = A->get_ifi(i+1,j+1,k+1);
	    fG = A->get_ifi(i+1,j+1,k);
	    cH = A->get_cell(i,j+1,k);
	    // Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);
	    
	    // Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.x + fC->fs->vel.x + fF->fs->vel.x + fG->fs->vel.x);
	           q_w = 0.25 * (cA->fs->vel.x + cD->fs->vel.x + cH->fs->vel.x + cE->fs->vel.x);
	           q_n = 0.25 * (cE->fs->vel.x + fF->fs->vel.x + fG->fs->vel.x + cH->fs->vel.x);
	           q_s = 0.25 * (cA->fs->vel.x + fB->fs->vel.x + fC->fs->vel.x + cD->fs->vel.x);
	           q_top = 0.25 * (cA->fs->vel.x + fB->fs->vel.x + cE->fs->vel.x + fF->fs->vel.x);
	           q_bottom = 0.25 * (cD->fs->vel.x + fC->fs->vel.x + fG->fs->vel.x + cH->fs->vel.x);
	           // Apply the divergence theorem.
	           sec_ctr->dudx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dudy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dudz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.y + fC->fs->vel.y + fF->fs->vel.y + fG->fs->vel.y);
	           q_w = 0.25 * (cA->fs->vel.y + cD->fs->vel.y + cH->fs->vel.y + cE->fs->vel.y);
	           q_n = 0.25 * (cE->fs->vel.y + fF->fs->vel.y + fG->fs->vel.y + cH->fs->vel.y);
	           q_s = 0.25 * (cA->fs->vel.y + fB->fs->vel.y + fC->fs->vel.y + cD->fs->vel.y);
	           q_top = 0.25 * (cA->fs->vel.y + fB->fs->vel.y + cE->fs->vel.y + fF->fs->vel.y);
	           q_bottom = 0.25 * (cD->fs->vel.y + fC->fs->vel.y + fG->fs->vel.y + cH->fs->vel.y);
	           // Apply the divergence theorem.
	           sec_ctr->dvdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dvdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dvdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.z + fC->fs->vel.z + fF->fs->vel.z + fG->fs->vel.z);
	           q_w = 0.25 * (cA->fs->vel.z + cD->fs->vel.z + cH->fs->vel.z + cE->fs->vel.z);
	           q_n = 0.25 * (cE->fs->vel.z + fF->fs->vel.z + fG->fs->vel.z + cH->fs->vel.z);
	           q_s = 0.25 * (cA->fs->vel.z + fB->fs->vel.z + fC->fs->vel.z + cD->fs->vel.z);
	           q_top = 0.25 * (cA->fs->vel.z + fB->fs->vel.z + cE->fs->vel.z + fF->fs->vel.z);
	           q_bottom = 0.25 * (cD->fs->vel.z + fC->fs->vel.z + fG->fs->vel.z + cH->fs->vel.z);
	           // Apply the divergence theorem.
	           sec_ctr->dwdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dwdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dwdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->tke + fC->fs->tke + fF->fs->tke + fG->fs->tke);
	           q_w = 0.25 * (cA->fs->tke + cD->fs->tke + cH->fs->tke + cE->fs->tke);
	           q_n = 0.25 * (cE->fs->tke + fF->fs->tke + fG->fs->tke + cH->fs->tke);
	           q_s = 0.25 * (cA->fs->tke + fB->fs->tke + fC->fs->tke + cD->fs->tke);
	           q_top = 0.25 * (cA->fs->tke + fB->fs->tke + cE->fs->tke + fF->fs->tke);
	           q_bottom = 0.25 * (cD->fs->tke + fC->fs->tke + fG->fs->tke + cH->fs->tke);
	           // Apply the divergence theorem.
	           sec_ctr->dtkedx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dtkedy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dtkedz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->omega + fC->fs->omega + fF->fs->omega + fG->fs->omega);
	           q_w = 0.25 * (cA->fs->omega + cD->fs->omega + cH->fs->omega + cE->fs->omega);
	           q_n = 0.25 * (cE->fs->omega + fF->fs->omega + fG->fs->omega + cH->fs->omega);
	           q_s = 0.25 * (cA->fs->omega + fB->fs->omega + fC->fs->omega + cD->fs->omega);
	           q_top = 0.25 * (cA->fs->omega + fB->fs->omega + cE->fs->omega + fF->fs->omega);
	           q_bottom = 0.25 * (cD->fs->omega + fC->fs->omega + fG->fs->omega + cH->fs->omega);
	           // Apply the divergence theorem.
	           sec_ctr->domegadx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->domegady = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->domegadz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);

	    
            for ( size_t itm=0; itm<ntm; ++itm ) {
                // Average property value on each face.
	           q_e = 0.25 * (fB->fs->gas->T[itm] + fC->fs->gas->T[itm] + fF->fs->gas->T[itm] + fG->fs->gas->T[itm]);
	           q_w = 0.25 * (cA->fs->gas->T[itm] + cD->fs->gas->T[itm] + cH->fs->gas->T[itm] + cE->fs->gas->T[itm]);
	           q_n = 0.25 * (cE->fs->gas->T[itm] + fF->fs->gas->T[itm] + fG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
	           q_s = 0.25 * (cA->fs->gas->T[itm] + fB->fs->gas->T[itm] + fC->fs->gas->T[itm] + cD->fs->gas->T[itm]);
	           q_top = 0.25 * (cA->fs->gas->T[itm] + fB->fs->gas->T[itm] + cE->fs->gas->T[itm] + fF->fs->gas->T[itm]);
	           q_bottom = 0.25 * (cD->fs->gas->T[itm] + fC->fs->gas->T[itm] + fG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
	           // Apply the divergence theorem.
	           sec_ctr->dTdx[itm] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dTdy[itm] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dTdz[itm] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		// Average property value on each face.
	           q_e = 0.25 * (fB->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + fF->fs->gas->massf[isp] + fG->fs->gas->massf[isp]);
	           q_w = 0.25 * (cA->fs->gas->massf[isp] + cD->fs->gas->massf[isp] + cH->fs->gas->massf[isp] + cE->fs->gas->massf[isp]);
	           q_n = 0.25 * (cE->fs->gas->massf[isp] + fF->fs->gas->massf[isp] + fG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
	           q_s = 0.25 * (cA->fs->gas->massf[isp] + fB->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + cD->fs->gas->massf[isp]);
	           q_top = 0.25 * (cA->fs->gas->massf[isp] + fB->fs->gas->massf[isp] + cE->fs->gas->massf[isp] + fF->fs->gas->massf[isp]);
	           q_bottom = 0.25 * (cD->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + fG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
	           // Apply the divergence theorem.
	           sec_ctr->dfdx[isp] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dfdy[isp] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dfdz[isp] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
	} // k loop
    } // j loop

    // West boundary surface.
    i = A->imin - 1;
    for (j = A->jmin; j <= A->jmax - 1; ++j) {
	for (k = A->kmin; k <= A->kmax - 1; ++k) {
	    sec_ctr = A->get_vtx(i+1,j+1,k+1);
	    vol_inv = 1.0 / sec_ctr->volume;
	    // The following are the corners of the secondary cell.
	    fA = A->get_ifi(i+1,j,k+1);
	    cB = A->get_cell(i+1,j,k+1);
	    cC = A->get_cell(i+1,j,k);
	    fD = A->get_ifi(i+1,j,k);
	    fE = A->get_ifi(i+1,j+1,k+1);
	    cF = A->get_cell(i+1,j+1,k+1);
	    cG = A->get_cell(i+1,j+1,k);
	    fH = A->get_ifi(i+1,j+1,k);
	    // Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);
	    
	    // Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.x + cC->fs->vel.x + cF->fs->vel.x + cG->fs->vel.x);
	           q_w = 0.25 * (fA->fs->vel.x + fD->fs->vel.x + fH->fs->vel.x + fE->fs->vel.x);
	           q_n = 0.25 * (fE->fs->vel.x + cF->fs->vel.x + cG->fs->vel.x + fH->fs->vel.x);
	           q_s = 0.25 * (fA->fs->vel.x + cB->fs->vel.x + cC->fs->vel.x + fD->fs->vel.x);
	           q_top = 0.25 * (fA->fs->vel.x + cB->fs->vel.x + fE->fs->vel.x + cF->fs->vel.x);
	           q_bottom = 0.25 * (fD->fs->vel.x + cC->fs->vel.x + cG->fs->vel.x + fH->fs->vel.x);
	           // Apply the divergence theorem.
	           sec_ctr->dudx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dudy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dudz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.y + cC->fs->vel.y + cF->fs->vel.y + cG->fs->vel.y);
	           q_w = 0.25 * (fA->fs->vel.y + fD->fs->vel.y + fH->fs->vel.y + fE->fs->vel.y);
	           q_n = 0.25 * (fE->fs->vel.y + cF->fs->vel.y + cG->fs->vel.y + fH->fs->vel.y);
	           q_s = 0.25 * (fA->fs->vel.y + cB->fs->vel.y + cC->fs->vel.y + fD->fs->vel.y);
	           q_top = 0.25 * (fA->fs->vel.y + cB->fs->vel.y + fE->fs->vel.y + cF->fs->vel.y);
	           q_bottom = 0.25 * (fD->fs->vel.y + cC->fs->vel.y + cG->fs->vel.y + fH->fs->vel.y);
	           // Apply the divergence theorem.
	           sec_ctr->dvdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dvdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dvdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.z + cC->fs->vel.z + cF->fs->vel.z + cG->fs->vel.z);
	           q_w = 0.25 * (fA->fs->vel.z + fD->fs->vel.z + fH->fs->vel.z + fE->fs->vel.z);
	           q_n = 0.25 * (fE->fs->vel.z + cF->fs->vel.z + cG->fs->vel.z + fH->fs->vel.z);
	           q_s = 0.25 * (fA->fs->vel.z + cB->fs->vel.z + cC->fs->vel.z + fD->fs->vel.z);
	           q_top = 0.25 * (fA->fs->vel.z + cB->fs->vel.z + fE->fs->vel.z + cF->fs->vel.z);
	           q_bottom = 0.25 * (fD->fs->vel.z + cC->fs->vel.z + cG->fs->vel.z + fH->fs->vel.z);
	           // Apply the divergence theorem.
	           sec_ctr->dwdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dwdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dwdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->tke + cC->fs->tke + cF->fs->tke + cG->fs->tke);
	           q_w = 0.25 * (fA->fs->tke + fD->fs->tke + fH->fs->tke + fE->fs->tke);
	           q_n = 0.25 * (fE->fs->tke + cF->fs->tke + cG->fs->tke + fH->fs->tke);
	           q_s = 0.25 * (fA->fs->tke + cB->fs->tke + cC->fs->tke + fD->fs->tke);
	           q_top = 0.25 * (fA->fs->tke + cB->fs->tke + fE->fs->tke + cF->fs->tke);
	           q_bottom = 0.25 * (fD->fs->tke + cC->fs->tke + cG->fs->tke + fH->fs->tke);
	           // Apply the divergence theorem.
	           sec_ctr->dtkedx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dtkedy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dtkedz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->omega + cC->fs->omega + cF->fs->omega + cG->fs->omega);
	           q_w = 0.25 * (fA->fs->omega + fD->fs->omega + fH->fs->omega + fE->fs->omega);
	           q_n = 0.25 * (fE->fs->omega + cF->fs->omega + cG->fs->omega + fH->fs->omega);
	           q_s = 0.25 * (fA->fs->omega + cB->fs->omega + cC->fs->omega + fD->fs->omega);
	           q_top = 0.25 * (fA->fs->omega + cB->fs->omega + fE->fs->omega + cF->fs->omega);
	           q_bottom = 0.25 * (fD->fs->omega + cC->fs->omega + cG->fs->omega + fH->fs->omega);
	           // Apply the divergence theorem.
	           sec_ctr->domegadx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->domegady = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->domegadz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);

	    
            for ( size_t itm=0; itm<ntm; ++itm ) {
	        // Average property value on each face.
	           q_e = 0.25 * (cB->fs->gas->T[itm] + cC->fs->gas->T[itm] + cF->fs->gas->T[itm] + cG->fs->gas->T[itm]);
	           q_w = 0.25 * (fA->fs->gas->T[itm] + fD->fs->gas->T[itm] + fH->fs->gas->T[itm] + fE->fs->gas->T[itm]);
	           q_n = 0.25 * (fE->fs->gas->T[itm] + cF->fs->gas->T[itm] + cG->fs->gas->T[itm] + fH->fs->gas->T[itm]);
	           q_s = 0.25 * (fA->fs->gas->T[itm] + cB->fs->gas->T[itm] + cC->fs->gas->T[itm] + fD->fs->gas->T[itm]);
	           q_top = 0.25 * (fA->fs->gas->T[itm] + cB->fs->gas->T[itm] + fE->fs->gas->T[itm] + cF->fs->gas->T[itm]);
	           q_bottom = 0.25 * (fD->fs->gas->T[itm] + cC->fs->gas->T[itm] + cG->fs->gas->T[itm] + fH->fs->gas->T[itm]);
	           // Apply the divergence theorem.
	           sec_ctr->dTdx[itm] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dTdy[itm] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dTdz[itm] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		// Average property value on each face.
	           q_e = 0.25 * (cB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + cG->fs->gas->massf[isp]);
	           q_w = 0.25 * (fA->fs->gas->massf[isp] + fD->fs->gas->massf[isp] + fH->fs->gas->massf[isp] + fE->fs->gas->massf[isp]);
	           q_n = 0.25 * (fE->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + fH->fs->gas->massf[isp]);
	           q_s = 0.25 * (fA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + fD->fs->gas->massf[isp]);
	           q_top = 0.25 * (fA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + fE->fs->gas->massf[isp] + cF->fs->gas->massf[isp]);
	           q_bottom = 0.25 * (fD->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + fH->fs->gas->massf[isp]);
	           // Apply the divergence theorem.
	           sec_ctr->dfdx[isp] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dfdy[isp] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dfdz[isp] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
	} // k loop
    } // j loop

    // North boundary surface
    j = A->jmax;
    for (i = A->imin; i <= A->imax - 1; ++i) {
	for (k = A->kmin; k <= A->kmax - 1; ++k) {
	    sec_ctr = A->get_vtx(i+1,j+1,k+1);
	    vol_inv = 1.0 / sec_ctr->volume;
	    // The following are the corners of the secondary cell.
	    cA = A->get_cell(i,j,k+1);
	    cB = A->get_cell(i+1,j,k+1);
	    cC = A->get_cell(i+1,j,k);
	    cD = A->get_cell(i,j,k);
	    fE = A->get_ifj(i,j+1,k+1);
	    fF = A->get_ifj(i+1,j+1,k+1);
	    fG = A->get_ifj(i+1,j+1,k);
	    fH = A->get_ifj(i,j+1,k);
	    // Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);
	    
	    // Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.x + cC->fs->vel.x + fF->fs->vel.x + fG->fs->vel.x);
	           q_w = 0.25 * (cA->fs->vel.x + cD->fs->vel.x + fH->fs->vel.x + fE->fs->vel.x);
	           q_n = 0.25 * (fE->fs->vel.x + fF->fs->vel.x + fG->fs->vel.x + fH->fs->vel.x);
	           q_s = 0.25 * (cA->fs->vel.x + cB->fs->vel.x + cC->fs->vel.x + cD->fs->vel.x);
	           q_top = 0.25 * (cA->fs->vel.x + cB->fs->vel.x + fE->fs->vel.x + fF->fs->vel.x);
	           q_bottom = 0.25 * (cD->fs->vel.x + cC->fs->vel.x + fG->fs->vel.x + fH->fs->vel.x);
	           // Apply the divergence theorem.
	           sec_ctr->dudx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dudy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dudz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.y + cC->fs->vel.y + fF->fs->vel.y + fG->fs->vel.y);
	           q_w = 0.25 * (cA->fs->vel.y + cD->fs->vel.y + fH->fs->vel.y + fE->fs->vel.y);
	           q_n = 0.25 * (fE->fs->vel.y + fF->fs->vel.y + fG->fs->vel.y + fH->fs->vel.y);
	           q_s = 0.25 * (cA->fs->vel.y + cB->fs->vel.y + cC->fs->vel.y + cD->fs->vel.y);
	           q_top = 0.25 * (cA->fs->vel.y + cB->fs->vel.y + fE->fs->vel.y + fF->fs->vel.y);
	           q_bottom = 0.25 * (cD->fs->vel.y + cC->fs->vel.y + fG->fs->vel.y + fH->fs->vel.y);
	           // Apply the divergence theorem.
	           sec_ctr->dvdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dvdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dvdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.z + cC->fs->vel.z + fF->fs->vel.z + fG->fs->vel.z);
	           q_w = 0.25 * (cA->fs->vel.z + cD->fs->vel.z + fH->fs->vel.z + fE->fs->vel.z);
	           q_n = 0.25 * (fE->fs->vel.z + fF->fs->vel.z + fG->fs->vel.z + fH->fs->vel.z);
	           q_s = 0.25 * (cA->fs->vel.z + cB->fs->vel.z + cC->fs->vel.z + cD->fs->vel.z);
	           q_top = 0.25 * (cA->fs->vel.z + cB->fs->vel.z + fE->fs->vel.z + fF->fs->vel.z);
	           q_bottom = 0.25 * (cD->fs->vel.z + cC->fs->vel.z + fG->fs->vel.z + fH->fs->vel.z);
	           // Apply the divergence theorem.
	           sec_ctr->dwdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dwdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dwdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->tke + cC->fs->tke + fF->fs->tke + fG->fs->tke);
	           q_w = 0.25 * (cA->fs->tke + cD->fs->tke + fH->fs->tke + fE->fs->tke);
	           q_n = 0.25 * (fE->fs->tke + fF->fs->tke + fG->fs->tke + fH->fs->tke);
	           q_s = 0.25 * (cA->fs->tke + cB->fs->tke + cC->fs->tke + cD->fs->tke);
	           q_top = 0.25 * (cA->fs->tke + cB->fs->tke + fE->fs->tke + fF->fs->tke);
	           q_bottom = 0.25 * (cD->fs->tke + cC->fs->tke + fG->fs->tke + fH->fs->tke);
	           // Apply the divergence theorem.
	           sec_ctr->dtkedx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dtkedy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dtkedz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->omega + cC->fs->omega + fF->fs->omega + fG->fs->omega);
	           q_w = 0.25 * (cA->fs->omega + cD->fs->omega + fH->fs->omega + fE->fs->omega);
	           q_n = 0.25 * (fE->fs->omega + fF->fs->omega + fG->fs->omega + fH->fs->omega);
	           q_s = 0.25 * (cA->fs->omega + cB->fs->omega + cC->fs->omega + cD->fs->omega);
	           q_top = 0.25 * (cA->fs->omega + cB->fs->omega + fE->fs->omega + fF->fs->omega);
	           q_bottom = 0.25 * (cD->fs->omega + cC->fs->omega + fG->fs->omega + fH->fs->omega);
	           // Apply the divergence theorem.
	           sec_ctr->domegadx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->domegady = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->domegadz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);

	    
            for ( size_t itm=0; itm<ntm; ++itm ) {
	        // Average property value on each face.
	           q_e = 0.25 * (cB->fs->gas->T[itm] + cC->fs->gas->T[itm] + fF->fs->gas->T[itm] + fG->fs->gas->T[itm]);
	           q_w = 0.25 * (cA->fs->gas->T[itm] + cD->fs->gas->T[itm] + fH->fs->gas->T[itm] + fE->fs->gas->T[itm]);
	           q_n = 0.25 * (fE->fs->gas->T[itm] + fF->fs->gas->T[itm] + fG->fs->gas->T[itm] + fH->fs->gas->T[itm]);
	           q_s = 0.25 * (cA->fs->gas->T[itm] + cB->fs->gas->T[itm] + cC->fs->gas->T[itm] + cD->fs->gas->T[itm]);
	           q_top = 0.25 * (cA->fs->gas->T[itm] + cB->fs->gas->T[itm] + fE->fs->gas->T[itm] + fF->fs->gas->T[itm]);
	           q_bottom = 0.25 * (cD->fs->gas->T[itm] + cC->fs->gas->T[itm] + fG->fs->gas->T[itm] + fH->fs->gas->T[itm]);
	           // Apply the divergence theorem.
	           sec_ctr->dTdx[itm] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dTdy[itm] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dTdz[itm] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		// Average property value on each face.
	           q_e = 0.25 * (cB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + fF->fs->gas->massf[isp] + fG->fs->gas->massf[isp]);
	           q_w = 0.25 * (cA->fs->gas->massf[isp] + cD->fs->gas->massf[isp] + fH->fs->gas->massf[isp] + fE->fs->gas->massf[isp]);
	           q_n = 0.25 * (fE->fs->gas->massf[isp] + fF->fs->gas->massf[isp] + fG->fs->gas->massf[isp] + fH->fs->gas->massf[isp]);
	           q_s = 0.25 * (cA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cD->fs->gas->massf[isp]);
	           q_top = 0.25 * (cA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + fE->fs->gas->massf[isp] + fF->fs->gas->massf[isp]);
	           q_bottom = 0.25 * (cD->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + fG->fs->gas->massf[isp] + fH->fs->gas->massf[isp]);
	           // Apply the divergence theorem.
	           sec_ctr->dfdx[isp] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dfdy[isp] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dfdz[isp] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
        } // k loop
    } // i loop

    // South boundary surface
    j = A->jmin - 1;
    for (i = A->imin; i <= A->imax - 1; ++i) {
	for (k = A->kmin; k <= A->kmax - 1; ++k) {
	    sec_ctr = A->get_vtx(i+1,j+1,k+1);
	    vol_inv = 1.0 / sec_ctr->volume;
	    // The following are the corners of the secondary cell.
	    fA = A->get_ifj(i,j+1,k+1);
	    fB = A->get_ifj(i+1,j+1,k+1);
	    fC = A->get_ifj(i+1,j+1,k);
	    fD = A->get_ifj(i,j+1,k);
	    cE = A->get_cell(i,j+1,k+1);
	    cF = A->get_cell(i+1,j+1,k+1);
	    cG = A->get_cell(i+1,j+1,k);
	    cH = A->get_cell(i,j+1,k);
	    // Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);
	    
	    // Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.x + fC->fs->vel.x + cF->fs->vel.x + cG->fs->vel.x);
	           q_w = 0.25 * (fA->fs->vel.x + fD->fs->vel.x + cH->fs->vel.x + cE->fs->vel.x);
	           q_n = 0.25 * (cE->fs->vel.x + cF->fs->vel.x + cG->fs->vel.x + cH->fs->vel.x);
	           q_s = 0.25 * (fA->fs->vel.x + fB->fs->vel.x + fC->fs->vel.x + fD->fs->vel.x);
	           q_top = 0.25 * (fA->fs->vel.x + fB->fs->vel.x + cE->fs->vel.x + cF->fs->vel.x);
	           q_bottom = 0.25 * (fD->fs->vel.x + fC->fs->vel.x + cG->fs->vel.x + cH->fs->vel.x);
	           // Apply the divergence theorem.
	           sec_ctr->dudx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dudy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dudz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.y + fC->fs->vel.y + cF->fs->vel.y + cG->fs->vel.y);
	           q_w = 0.25 * (fA->fs->vel.y + fD->fs->vel.y + cH->fs->vel.y + cE->fs->vel.y);
	           q_n = 0.25 * (cE->fs->vel.y + cF->fs->vel.y + cG->fs->vel.y + cH->fs->vel.y);
	           q_s = 0.25 * (fA->fs->vel.y + fB->fs->vel.y + fC->fs->vel.y + fD->fs->vel.y);
	           q_top = 0.25 * (fA->fs->vel.y + fB->fs->vel.y + cE->fs->vel.y + cF->fs->vel.y);
	           q_bottom = 0.25 * (fD->fs->vel.y + fC->fs->vel.y + cG->fs->vel.y + cH->fs->vel.y);
	           // Apply the divergence theorem.
	           sec_ctr->dvdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dvdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dvdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.z + fC->fs->vel.z + cF->fs->vel.z + cG->fs->vel.z);
	           q_w = 0.25 * (fA->fs->vel.z + fD->fs->vel.z + cH->fs->vel.z + cE->fs->vel.z);
	           q_n = 0.25 * (cE->fs->vel.z + cF->fs->vel.z + cG->fs->vel.z + cH->fs->vel.z);
	           q_s = 0.25 * (fA->fs->vel.z + fB->fs->vel.z + fC->fs->vel.z + fD->fs->vel.z);
	           q_top = 0.25 * (fA->fs->vel.z + fB->fs->vel.z + cE->fs->vel.z + cF->fs->vel.z);
	           q_bottom = 0.25 * (fD->fs->vel.z + fC->fs->vel.z + cG->fs->vel.z + cH->fs->vel.z);
	           // Apply the divergence theorem.
	           sec_ctr->dwdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dwdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dwdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->tke + fC->fs->tke + cF->fs->tke + cG->fs->tke);
	           q_w = 0.25 * (fA->fs->tke + fD->fs->tke + cH->fs->tke + cE->fs->tke);
	           q_n = 0.25 * (cE->fs->tke + cF->fs->tke + cG->fs->tke + cH->fs->tke);
	           q_s = 0.25 * (fA->fs->tke + fB->fs->tke + fC->fs->tke + fD->fs->tke);
	           q_top = 0.25 * (fA->fs->tke + fB->fs->tke + cE->fs->tke + cF->fs->tke);
	           q_bottom = 0.25 * (fD->fs->tke + fC->fs->tke + cG->fs->tke + cH->fs->tke);
	           // Apply the divergence theorem.
	           sec_ctr->dtkedx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dtkedy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dtkedz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->omega + fC->fs->omega + cF->fs->omega + cG->fs->omega);
	           q_w = 0.25 * (fA->fs->omega + fD->fs->omega + cH->fs->omega + cE->fs->omega);
	           q_n = 0.25 * (cE->fs->omega + cF->fs->omega + cG->fs->omega + cH->fs->omega);
	           q_s = 0.25 * (fA->fs->omega + fB->fs->omega + fC->fs->omega + fD->fs->omega);
	           q_top = 0.25 * (fA->fs->omega + fB->fs->omega + cE->fs->omega + cF->fs->omega);
	           q_bottom = 0.25 * (fD->fs->omega + fC->fs->omega + cG->fs->omega + cH->fs->omega);
	           // Apply the divergence theorem.
	           sec_ctr->domegadx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->domegady = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->domegadz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);

	    
            for ( size_t itm=0; itm<ntm; ++itm ) {
	        // Average property value on each face.
	           q_e = 0.25 * (fB->fs->gas->T[itm] + fC->fs->gas->T[itm] + cF->fs->gas->T[itm] + cG->fs->gas->T[itm]);
	           q_w = 0.25 * (fA->fs->gas->T[itm] + fD->fs->gas->T[itm] + cH->fs->gas->T[itm] + cE->fs->gas->T[itm]);
	           q_n = 0.25 * (cE->fs->gas->T[itm] + cF->fs->gas->T[itm] + cG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
	           q_s = 0.25 * (fA->fs->gas->T[itm] + fB->fs->gas->T[itm] + fC->fs->gas->T[itm] + fD->fs->gas->T[itm]);
	           q_top = 0.25 * (fA->fs->gas->T[itm] + fB->fs->gas->T[itm] + cE->fs->gas->T[itm] + cF->fs->gas->T[itm]);
	           q_bottom = 0.25 * (fD->fs->gas->T[itm] + fC->fs->gas->T[itm] + cG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
	           // Apply the divergence theorem.
	           sec_ctr->dTdx[itm] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dTdy[itm] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dTdz[itm] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		// Average property value on each face.
	           q_e = 0.25 * (fB->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + cG->fs->gas->massf[isp]);
	           q_w = 0.25 * (fA->fs->gas->massf[isp] + fD->fs->gas->massf[isp] + cH->fs->gas->massf[isp] + cE->fs->gas->massf[isp]);
	           q_n = 0.25 * (cE->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
	           q_s = 0.25 * (fA->fs->gas->massf[isp] + fB->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + fD->fs->gas->massf[isp]);
	           q_top = 0.25 * (fA->fs->gas->massf[isp] + fB->fs->gas->massf[isp] + cE->fs->gas->massf[isp] + cF->fs->gas->massf[isp]);
	           q_bottom = 0.25 * (fD->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
	           // Apply the divergence theorem.
	           sec_ctr->dfdx[isp] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dfdy[isp] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dfdz[isp] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
	} // k loop
    } // i loop

    // Top boundary surface
    k = A->kmax;
    for (i = A->imin; i <= A->imax - 1; ++i) {
	for (j = A->jmin; j <= A->jmax - 1; ++j) {
	    sec_ctr = A->get_vtx(i+1,j+1,k+1);
	    vol_inv = 1.0 / sec_ctr->volume;
	    // The following are the corners of the secondary cell.
	    fA = A->get_ifk(i,j,k+1);
	    fB = A->get_ifk(i+1,j,k+1);
	    cC = A->get_cell(i+1,j,k);
	    cD = A->get_cell(i,j,k);
	    fE = A->get_ifk(i,j+1,k+1);
	    fF = A->get_ifk(i+1,j+1,k+1);
	    cG = A->get_cell(i+1,j+1,k);
	    cH = A->get_cell(i,j+1,k);
	    // Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);
	    
	    // Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.x + cC->fs->vel.x + fF->fs->vel.x + cG->fs->vel.x);
	           q_w = 0.25 * (fA->fs->vel.x + cD->fs->vel.x + cH->fs->vel.x + fE->fs->vel.x);
	           q_n = 0.25 * (fE->fs->vel.x + fF->fs->vel.x + cG->fs->vel.x + cH->fs->vel.x);
	           q_s = 0.25 * (fA->fs->vel.x + fB->fs->vel.x + cC->fs->vel.x + cD->fs->vel.x);
	           q_top = 0.25 * (fA->fs->vel.x + fB->fs->vel.x + fE->fs->vel.x + fF->fs->vel.x);
	           q_bottom = 0.25 * (cD->fs->vel.x + cC->fs->vel.x + cG->fs->vel.x + cH->fs->vel.x);
	           // Apply the divergence theorem.
	           sec_ctr->dudx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dudy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dudz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.y + cC->fs->vel.y + fF->fs->vel.y + cG->fs->vel.y);
	           q_w = 0.25 * (fA->fs->vel.y + cD->fs->vel.y + cH->fs->vel.y + fE->fs->vel.y);
	           q_n = 0.25 * (fE->fs->vel.y + fF->fs->vel.y + cG->fs->vel.y + cH->fs->vel.y);
	           q_s = 0.25 * (fA->fs->vel.y + fB->fs->vel.y + cC->fs->vel.y + cD->fs->vel.y);
	           q_top = 0.25 * (fA->fs->vel.y + fB->fs->vel.y + fE->fs->vel.y + fF->fs->vel.y);
	           q_bottom = 0.25 * (cD->fs->vel.y + cC->fs->vel.y + cG->fs->vel.y + cH->fs->vel.y);
	           // Apply the divergence theorem.
	           sec_ctr->dvdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dvdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dvdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->vel.z + cC->fs->vel.z + fF->fs->vel.z + cG->fs->vel.z);
	           q_w = 0.25 * (fA->fs->vel.z + cD->fs->vel.z + cH->fs->vel.z + fE->fs->vel.z);
	           q_n = 0.25 * (fE->fs->vel.z + fF->fs->vel.z + cG->fs->vel.z + cH->fs->vel.z);
	           q_s = 0.25 * (fA->fs->vel.z + fB->fs->vel.z + cC->fs->vel.z + cD->fs->vel.z);
	           q_top = 0.25 * (fA->fs->vel.z + fB->fs->vel.z + fE->fs->vel.z + fF->fs->vel.z);
	           q_bottom = 0.25 * (cD->fs->vel.z + cC->fs->vel.z + cG->fs->vel.z + cH->fs->vel.z);
	           // Apply the divergence theorem.
	           sec_ctr->dwdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dwdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dwdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->tke + cC->fs->tke + fF->fs->tke + cG->fs->tke);
	           q_w = 0.25 * (fA->fs->tke + cD->fs->tke + cH->fs->tke + fE->fs->tke);
	           q_n = 0.25 * (fE->fs->tke + fF->fs->tke + cG->fs->tke + cH->fs->tke);
	           q_s = 0.25 * (fA->fs->tke + fB->fs->tke + cC->fs->tke + cD->fs->tke);
	           q_top = 0.25 * (fA->fs->tke + fB->fs->tke + fE->fs->tke + fF->fs->tke);
	           q_bottom = 0.25 * (cD->fs->tke + cC->fs->tke + cG->fs->tke + cH->fs->tke);
	           // Apply the divergence theorem.
	           sec_ctr->dtkedx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dtkedy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dtkedz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (fB->fs->omega + cC->fs->omega + fF->fs->omega + cG->fs->omega);
	           q_w = 0.25 * (fA->fs->omega + cD->fs->omega + cH->fs->omega + fE->fs->omega);
	           q_n = 0.25 * (fE->fs->omega + fF->fs->omega + cG->fs->omega + cH->fs->omega);
	           q_s = 0.25 * (fA->fs->omega + fB->fs->omega + cC->fs->omega + cD->fs->omega);
	           q_top = 0.25 * (fA->fs->omega + fB->fs->omega + fE->fs->omega + fF->fs->omega);
	           q_bottom = 0.25 * (cD->fs->omega + cC->fs->omega + cG->fs->omega + cH->fs->omega);
	           // Apply the divergence theorem.
	           sec_ctr->domegadx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->domegady = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->domegadz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);

	    
            for ( size_t itm=0; itm<ntm; ++itm ) {
	        // Average property value on each face.
	           q_e = 0.25 * (fB->fs->gas->T[itm] + cC->fs->gas->T[itm] + fF->fs->gas->T[itm] + cG->fs->gas->T[itm]);
	           q_w = 0.25 * (fA->fs->gas->T[itm] + cD->fs->gas->T[itm] + cH->fs->gas->T[itm] + fE->fs->gas->T[itm]);
	           q_n = 0.25 * (fE->fs->gas->T[itm] + fF->fs->gas->T[itm] + cG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
	           q_s = 0.25 * (fA->fs->gas->T[itm] + fB->fs->gas->T[itm] + cC->fs->gas->T[itm] + cD->fs->gas->T[itm]);
	           q_top = 0.25 * (fA->fs->gas->T[itm] + fB->fs->gas->T[itm] + fE->fs->gas->T[itm] + fF->fs->gas->T[itm]);
	           q_bottom = 0.25 * (cD->fs->gas->T[itm] + cC->fs->gas->T[itm] + cG->fs->gas->T[itm] + cH->fs->gas->T[itm]);
	           // Apply the divergence theorem.
	           sec_ctr->dTdx[itm] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dTdy[itm] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dTdz[itm] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		// Average property value on each face.
	           q_e = 0.25 * (fB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + fF->fs->gas->massf[isp] + cG->fs->gas->massf[isp]);
	           q_w = 0.25 * (fA->fs->gas->massf[isp] + cD->fs->gas->massf[isp] + cH->fs->gas->massf[isp] + fE->fs->gas->massf[isp]);
	           q_n = 0.25 * (fE->fs->gas->massf[isp] + fF->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
	           q_s = 0.25 * (fA->fs->gas->massf[isp] + fB->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cD->fs->gas->massf[isp]);
	           q_top = 0.25 * (fA->fs->gas->massf[isp] + fB->fs->gas->massf[isp] + fE->fs->gas->massf[isp] + fF->fs->gas->massf[isp]);
	           q_bottom = 0.25 * (cD->fs->gas->massf[isp] + cC->fs->gas->massf[isp] + cG->fs->gas->massf[isp] + cH->fs->gas->massf[isp]);
	           // Apply the divergence theorem.
	           sec_ctr->dfdx[isp] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dfdy[isp] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dfdz[isp] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
        } // j loop
    } // i loop

    // Bottom boundary surface
    k = A->kmin - 1;
    for (i = A->imin; i <= A->imax - 1; ++i) {
        for (j = A->jmin; j <= A->jmax - 1; ++j) {
	    sec_ctr = A->get_vtx(i+1,j+1,k+1);
	    vol_inv = 1.0 / sec_ctr->volume;
	    // The following are the corners of the secondary cell.
	    cA = A->get_cell(i,j,k+1);
	    cB = A->get_cell(i+1,j,k+1);
	    fC = A->get_ifk(i+1,j,k+1);
	    fD = A->get_ifk(i,j,k+1);
	    cE = A->get_cell(i,j+1,k+1);
	    cF = A->get_cell(i+1,j+1,k+1);
	    fG = A->get_ifk(i+1,j+1,k+1);
	    fH = A->get_ifk(i,j+1,k+1);
	    // Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);
	    
	    // Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.x + fC->fs->vel.x + cF->fs->vel.x + fG->fs->vel.x);
	           q_w = 0.25 * (cA->fs->vel.x + fD->fs->vel.x + fH->fs->vel.x + cE->fs->vel.x);
	           q_n = 0.25 * (cE->fs->vel.x + cF->fs->vel.x + fG->fs->vel.x + fH->fs->vel.x);
	           q_s = 0.25 * (cA->fs->vel.x + cB->fs->vel.x + fC->fs->vel.x + fD->fs->vel.x);
	           q_top = 0.25 * (cA->fs->vel.x + cB->fs->vel.x + cE->fs->vel.x + cF->fs->vel.x);
	           q_bottom = 0.25 * (fD->fs->vel.x + fC->fs->vel.x + fG->fs->vel.x + fH->fs->vel.x);
	           // Apply the divergence theorem.
	           sec_ctr->dudx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dudy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dudz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.y + fC->fs->vel.y + cF->fs->vel.y + fG->fs->vel.y);
	           q_w = 0.25 * (cA->fs->vel.y + fD->fs->vel.y + fH->fs->vel.y + cE->fs->vel.y);
	           q_n = 0.25 * (cE->fs->vel.y + cF->fs->vel.y + fG->fs->vel.y + fH->fs->vel.y);
	           q_s = 0.25 * (cA->fs->vel.y + cB->fs->vel.y + fC->fs->vel.y + fD->fs->vel.y);
	           q_top = 0.25 * (cA->fs->vel.y + cB->fs->vel.y + cE->fs->vel.y + cF->fs->vel.y);
	           q_bottom = 0.25 * (fD->fs->vel.y + fC->fs->vel.y + fG->fs->vel.y + fH->fs->vel.y);
	           // Apply the divergence theorem.
	           sec_ctr->dvdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dvdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dvdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->vel.z + fC->fs->vel.z + cF->fs->vel.z + fG->fs->vel.z);
	           q_w = 0.25 * (cA->fs->vel.z + fD->fs->vel.z + fH->fs->vel.z + cE->fs->vel.z);
	           q_n = 0.25 * (cE->fs->vel.z + cF->fs->vel.z + fG->fs->vel.z + fH->fs->vel.z);
	           q_s = 0.25 * (cA->fs->vel.z + cB->fs->vel.z + fC->fs->vel.z + fD->fs->vel.z);
	           q_top = 0.25 * (cA->fs->vel.z + cB->fs->vel.z + cE->fs->vel.z + cF->fs->vel.z);
	           q_bottom = 0.25 * (fD->fs->vel.z + fC->fs->vel.z + fG->fs->vel.z + fH->fs->vel.z);
	           // Apply the divergence theorem.
	           sec_ctr->dwdx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dwdy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dwdz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->tke + fC->fs->tke + cF->fs->tke + fG->fs->tke);
	           q_w = 0.25 * (cA->fs->tke + fD->fs->tke + fH->fs->tke + cE->fs->tke);
	           q_n = 0.25 * (cE->fs->tke + cF->fs->tke + fG->fs->tke + fH->fs->tke);
	           q_s = 0.25 * (cA->fs->tke + cB->fs->tke + fC->fs->tke + fD->fs->tke);
	           q_top = 0.25 * (cA->fs->tke + cB->fs->tke + cE->fs->tke + cF->fs->tke);
	           q_bottom = 0.25 * (fD->fs->tke + fC->fs->tke + fG->fs->tke + fH->fs->tke);
	           // Apply the divergence theorem.
	           sec_ctr->dtkedx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dtkedy = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dtkedz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
	// Average property value on each face.
	           q_e = 0.25 * (cB->fs->omega + fC->fs->omega + cF->fs->omega + fG->fs->omega);
	           q_w = 0.25 * (cA->fs->omega + fD->fs->omega + fH->fs->omega + cE->fs->omega);
	           q_n = 0.25 * (cE->fs->omega + cF->fs->omega + fG->fs->omega + fH->fs->omega);
	           q_s = 0.25 * (cA->fs->omega + cB->fs->omega + fC->fs->omega + fD->fs->omega);
	           q_top = 0.25 * (cA->fs->omega + cB->fs->omega + cE->fs->omega + cF->fs->omega);
	           q_bottom = 0.25 * (fD->fs->omega + fC->fs->omega + fG->fs->omega + fH->fs->omega);
	           // Apply the divergence theorem.
	           sec_ctr->domegadx = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->domegady = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->domegadz = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);

	    
            for ( size_t itm=0; itm<ntm; ++itm ) {
	        // Average property value on each face.
	           q_e = 0.25 * (cB->fs->gas->T[itm] + fC->fs->gas->T[itm] + cF->fs->gas->T[itm] + fG->fs->gas->T[itm]);
	           q_w = 0.25 * (cA->fs->gas->T[itm] + fD->fs->gas->T[itm] + fH->fs->gas->T[itm] + cE->fs->gas->T[itm]);
	           q_n = 0.25 * (cE->fs->gas->T[itm] + cF->fs->gas->T[itm] + fG->fs->gas->T[itm] + fH->fs->gas->T[itm]);
	           q_s = 0.25 * (cA->fs->gas->T[itm] + cB->fs->gas->T[itm] + fC->fs->gas->T[itm] + fD->fs->gas->T[itm]);
	           q_top = 0.25 * (cA->fs->gas->T[itm] + cB->fs->gas->T[itm] + cE->fs->gas->T[itm] + cF->fs->gas->T[itm]);
	           q_bottom = 0.25 * (fD->fs->gas->T[itm] + fC->fs->gas->T[itm] + fG->fs->gas->T[itm] + fH->fs->gas->T[itm]);
	           // Apply the divergence theorem.
	           sec_ctr->dTdx[itm] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dTdy[itm] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dTdz[itm] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
  	    for( size_t isp = 0; isp < nsp; ++isp ) {
		// Average property value on each face.
	           q_e = 0.25 * (cB->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + fG->fs->gas->massf[isp]);
	           q_w = 0.25 * (cA->fs->gas->massf[isp] + fD->fs->gas->massf[isp] + fH->fs->gas->massf[isp] + cE->fs->gas->massf[isp]);
	           q_n = 0.25 * (cE->fs->gas->massf[isp] + cF->fs->gas->massf[isp] + fG->fs->gas->massf[isp] + fH->fs->gas->massf[isp]);
	           q_s = 0.25 * (cA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + fD->fs->gas->massf[isp]);
	           q_top = 0.25 * (cA->fs->gas->massf[isp] + cB->fs->gas->massf[isp] + cE->fs->gas->massf[isp] + cF->fs->gas->massf[isp]);
	           q_bottom = 0.25 * (fD->fs->gas->massf[isp] + fC->fs->gas->massf[isp] + fG->fs->gas->massf[isp] + fH->fs->gas->massf[isp]);
	           // Apply the divergence theorem.
	           sec_ctr->dfdx[isp] = vol_inv * (q_e*east->area*east->n.x - q_w*west->area*west->n.x +
		        q_n*north->area*north->n.x - q_s*south->area*south->n.x +
		        q_top*top->area*top->n.x - q_bottom*bottom->area*bottom->n.x);
	           sec_ctr->dfdy[isp] = vol_inv * (q_e*east->area*east->n.y - q_w*west->area*west->n.y +
		        q_n*north->area*north->n.y - q_s*south->area*south->n.y +
		        q_top*top->area*top->n.y - q_bottom*bottom->area*bottom->n.y);
	           sec_ctr->dfdz[isp] = vol_inv * (q_e*east->area*east->n.z - q_w*west->area*west->n.z +
		        q_n*north->area*north->n.z - q_s*south->area*south->n.z +
		        q_top*top->area*top->n.z - q_bottom*bottom->area*bottom->n.z);
            }
 	} // j loop
    } // i loop

    /*
     * ...and process the edge and corner values separately.
     */
    viscous_derivatives_edge_3D(A);
    viscous_derivatives_corners_3D(A);
    return SUCCESS;
} // end viscous_derivatives_3D()


/** \brief Corner derivatives from averages of near-by surface vertices.
 *
 * Fit a linear function to the nearest data points as per notebook Jan-2009.
 */
int viscous_derivatives_corners_3D(Block *bdp)
{
    size_t i, j, k;
    FV_Vertex *vtx;
    FV_Interface *a, *b, *d;
    FV_Cell *c;
    double xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd;
    double fa, fb, fc, fd, denom;

    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t ntm = gmodel->get_number_of_modes();

    // South-West-Bottom corner [0]
    i = bdp->imin; j = bdp->jmin; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j,k);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k);
        xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
                    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

        for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    // South-East-Bottom corner [1]
    i = bdp->imax; j = bdp->jmin; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j,k);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k);
    xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

    for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    // North-East-Bottom corner [2]
    i = bdp->imax; j = bdp->jmax; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j+1,k);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k);
    xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

    for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    // North-West-Bottom corner [3]
    i = bdp->imin; j = bdp->jmax; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j+1,k);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k);
    xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

    for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    // South-West-Top corner [4]
    i = bdp->imin; j = bdp->jmin; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j,k+1);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k+1);
    xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

    for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    // South-East-Top corner [5]
    i = bdp->imax; j = bdp->jmin; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j,k+1);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k+1);
    xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

    for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    // North-East-Top corner [6]
    i = bdp->imax; j = bdp->jmax; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j+1,k+1);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k+1);
    xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

    for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    // North-West-Top corner [7]
    i = bdp->imin; j = bdp->jmax; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j+1,k+1);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k+1);
    xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;
    fa = a->fs->vel.x; fb = b->fs->vel.x; fc = c->fs->vel.x; fd = d->fs->vel.x;
           vtx->dudx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dudy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dudz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.y; fb = b->fs->vel.y; fc = c->fs->vel.y; fd = d->fs->vel.y;
           vtx->dvdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dvdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dvdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->vel.z; fb = b->fs->vel.z; fc = c->fs->vel.z; fd = d->fs->vel.z;
           vtx->dwdx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dwdy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dwdz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->tke; fb = b->fs->tke; fc = c->fs->tke; fd = d->fs->tke;
           vtx->dtkedx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dtkedy = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dtkedz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
	fa = a->fs->omega; fb = b->fs->omega; fc = c->fs->omega; fd = d->fs->omega;
           vtx->domegadx = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->domegady = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->domegadz = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;

    for ( size_t itm=0; itm<ntm; ++itm ) {
        fa = a->fs->gas->T[itm]; fb = b->fs->gas->T[itm]; fc = c->fs->gas->T[itm]; fd = d->fs->gas->T[itm];
           vtx->dTdx[itm] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dTdy[itm] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dTdz[itm] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }
    for( size_t isp = 0; isp < nsp; ++isp ) {
        fa = a->fs->gas->massf[isp]; fb = b->fs->gas->massf[isp]; fc = c->fs->gas->massf[isp]; fd = d->fs->gas->massf[isp];
           vtx->dfdx[isp] = (fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za) / denom;
           vtx->dfdy[isp] = -(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za) / denom;
           vtx->dfdz[isp] = (fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya) / denom;
    }

    return SUCCESS;
} // end viscous_derivatives_corners_3D()


/** \brief Edge derivatives from a least-squares fit of near-by data.
 *
 * See workbook pages dated 07-Jan-2010.
 */
int viscous_derivatives_edge_3D(Block *bdp)
{
    size_t i, j, k, imin, jmin, kmin, imax, jmax, kmax;
    FV_Vertex *vtx;
    FV_Cell *ca, *cb;
    FV_Interface *fa, *fb, *fc, *fd;
    double ca_x, ca_y, ca_z, ca_f;
    double cb_x, cb_y, cb_z, cb_f;
    double fa_x, fa_y, fa_z, fa_f;
    double fb_x, fb_y, fb_z, fb_f;
    double fc_x, fc_y, fc_z, fc_f;
    double fd_x, fd_y, fd_z, fd_f;
    double S1, Sx, Sy, Sz, Sxx, Syx, Syy, Szx, Szy, Szz, Sf, Sfx, Sfy, Sfz;

    Gas_model *gmodel = get_gas_model_ptr();
    size_t nsp = gmodel->get_number_of_species();
    size_t ntm = gmodel->get_number_of_modes();

    Valmatrix A1( 4, 4);
    valarray<double> B1(4);
    valarray<double> x1(4);

    imin = bdp->imin; imax = bdp->imax;
    jmin = bdp->jmin; jmax = bdp->jmax;
    kmin = bdp->kmin; kmax = bdp->kmax;

    // Bottom-South edge [0]-->[1]
    j = jmin; k = kmin;    	
    for ( i = imin+1; i <= imax; ++i ) {
	vtx = bdp->get_vtx(i,j,k);
        ca = bdp->get_cell(i-1,j,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifj(i-1,j,k);
        fb = bdp->get_ifk(i-1,j,k);
        fc = bdp->get_ifj(i,j,k);
        fd = bdp->get_ifk(i,j,k);
	        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
                	                                        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // Bottom-North edge [3]-->[2]
    j = jmax; k = kmin;
    for ( i = imin+1; i <= imax; ++i ) {
	vtx = bdp->get_vtx(i,j+1,k);
        ca = bdp->get_cell(i-1,j,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifj(i-1,j+1,k);
        fb = bdp->get_ifk(i-1,j,k);
        fc = bdp->get_ifj(i,j+1,k);
        fd = bdp->get_ifk(i,j,k);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // Bottom-West edge [0]-->[3]
    i = imin; k = kmin;
    for ( j = jmin+1; j <= jmax; ++j ) {
	vtx = bdp->get_vtx(i,j,k);
        ca = bdp->get_cell(i,j-1,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i,j-1,k);
        fb = bdp->get_ifk(i,j-1,k);
        fc = bdp->get_ifi(i,j,k);
        fd = bdp->get_ifk(i,j,k);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // Bottom-East edge [1]-->[2]
    i = imax; k = kmin;
    for ( j = jmin+1; j <= jmax; ++j ) {
	vtx = bdp->get_vtx(i+1,j,k);
        ca = bdp->get_cell(i,j-1,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i+1,j-1,k);
        fb = bdp->get_ifk(i,j-1,k);
        fc = bdp->get_ifi(i+1,j,k);
        fd = bdp->get_ifk(i,j,k);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // Top-South edge [4]-->[5]
    j = jmin; k = kmax;
    for ( i = imin+1; i <= imax; ++i ) {
	vtx = bdp->get_vtx(i,j,k+1);
        ca = bdp->get_cell(i-1,j,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifj(i-1,j,k);
        fb = bdp->get_ifk(i-1,j,k+1);
        fc = bdp->get_ifj(i,j,k);
        fd = bdp->get_ifk(i,j,k+1);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // Top-North edge [7]-->[6]
    j = jmax; k = kmax;
    for ( i = imin+1; i <= imax; ++i ) {
	vtx = bdp->get_vtx(i,j+1,k+1);
        ca = bdp->get_cell(i-1,j,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifj(i-1,j+1,k);
        fb = bdp->get_ifk(i-1,j,k+1);
        fc = bdp->get_ifj(i,j+1,k);
        fd = bdp->get_ifk(i,j,k+1);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // Top-West edge [4]-->[7]
    i = imin; k = kmax;
    for ( j = jmin+1; j <= jmax; ++j ) {
	vtx = bdp->get_vtx(i,j,k+1);
        ca = bdp->get_cell(i,j-1,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i,j-1,k);
        fb = bdp->get_ifk(i,j-1,k+1);
        fc = bdp->get_ifi(i,j,k);
        fd = bdp->get_ifk(i,j,k+1);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // Top-East edge [5]-->[6]
    i = imax; k = kmax;
    for ( j = jmin+1; j <= jmax; ++j ) {
	vtx = bdp->get_vtx(i+1,j,k+1);
        ca = bdp->get_cell(i,j-1,k);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i+1,j-1,k);
        fb = bdp->get_ifk(i,j-1,k+1);
        fc = bdp->get_ifi(i+1,j,k);
        fd = bdp->get_ifk(i,j,k+1);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // South-West edge [0]-->[4]
    i = imin; j = jmin;
    for ( k = kmin+1; k <= kmax; ++k ) {
	vtx = bdp->get_vtx(i,j,k);
        ca = bdp->get_cell(i,j,k-1);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i,j,k-1);
        fb = bdp->get_ifj(i,j,k-1);
        fc = bdp->get_ifi(i,j,k);
        fd = bdp->get_ifj(i,j,k);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // South-East edge [1]-->[5]
    i = imax; j = jmin;
    for ( k = kmin+1; k <= kmax; ++k ) {
	vtx = bdp->get_vtx(i+1,j,k);
        ca = bdp->get_cell(i,j,k-1);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i+1,j,k-1);
        fb = bdp->get_ifj(i,j,k-1);
        fc = bdp->get_ifi(i+1,j,k);
        fd = bdp->get_ifj(i,j,k);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // North-East edge [2]-->[6]
    i = bdp->imax; j = bdp->jmax;
    for ( k = kmin+1; k <= kmax; ++k ) {
	vtx = bdp->get_vtx(i+1,j+1,k);
        ca = bdp->get_cell(i,j,k-1);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i+1,j,k-1);
        fb = bdp->get_ifj(i,j+1,k-1);
        fc = bdp->get_ifi(i+1,j,k);
        fd = bdp->get_ifj(i,j+1,k);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }

    // North-West edge [3]-->[7]
    i = imin; j = jmax;
    for ( k = kmin+1; k <= kmax; ++k ) {
	vtx = bdp->get_vtx(i,j+1,k);
        ca = bdp->get_cell(i,j,k-1);
        cb = bdp->get_cell(i,j,k);
        fa = bdp->get_ifi(i,j,k-1);
        fb = bdp->get_ifj(i,j+1,k-1);
        fc = bdp->get_ifi(i,j,k);
        fd = bdp->get_ifj(i,j+1,k);
        ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
                cb_x = cb->pos.x; cb_y = cb->pos.y; cb_z = cb->pos.z;
                fa_x = fa->pos.x; fa_y = fa->pos.y; fa_z = fa->pos.z;
                fb_x = fb->pos.x; fb_y = fb->pos.y; fb_z = fb->pos.z;
                fc_x = fc->pos.x; fc_y = fc->pos.y; fc_z = fc->pos.z;
                fd_x = fd->pos.x; fd_y = fd->pos.y; fd_z = fd->pos.z;
                S1 = 6;
	        Sx = ca_x + cb_x + fa_x + fb_x + fc_x + fd_x;
	        Sy = ca_y + cb_y + fa_y + fb_y + fc_y + fd_y;
	        Sz = ca_z + cb_z + fa_z + fb_z + fc_z + fd_z;
                Sxx = ca_x*ca_x + cb_x*cb_x + fa_x*fa_x + fb_x*fb_x + fc_x*fc_x + fd_x*fd_x;
                Syx = ca_y*ca_x + cb_y*cb_x + fa_y*fa_x + fb_y*fb_x + fc_y*fc_x + fd_y*fd_x;
                Syy = ca_y*ca_y + cb_y*cb_y + fa_y*fa_y + fb_y*fb_y + fc_y*fc_y + fd_y*fd_y;
                Szx = ca_z*ca_x + cb_z*cb_x + fa_z*fa_x + fb_z*fb_x + fc_z*fc_x + fd_z*fd_x;
                Szy = ca_z*ca_y + cb_z*cb_y + fa_z*fa_y + fb_z*fb_y + fc_z*fc_y + fd_z*fd_y;
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;
        ca_f = ca->fs->vel.x; cb_f = cb->fs->vel.x;
	        fa_f = fa->fs->vel.x; fb_f = fb->fs->vel.x; fc_f = fc->fs->vel.x; fd_f = fd->fs->vel.x;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dudx = x1[1];
                vtx->dudy = x1[2];
                vtx->dudz = x1[3];
	ca_f = ca->fs->vel.y; cb_f = cb->fs->vel.y;
	        fa_f = fa->fs->vel.y; fb_f = fb->fs->vel.y; fc_f = fc->fs->vel.y; fd_f = fd->fs->vel.y;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dvdx = x1[1];
                vtx->dvdy = x1[2];
                vtx->dvdz = x1[3];
	ca_f = ca->fs->vel.z; cb_f = cb->fs->vel.z;
	        fa_f = fa->fs->vel.z; fb_f = fb->fs->vel.z; fc_f = fc->fs->vel.z; fd_f = fd->fs->vel.z;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dwdx = x1[1];
                vtx->dwdy = x1[2];
                vtx->dwdz = x1[3];
	ca_f = ca->fs->tke; cb_f = cb->fs->tke;
	        fa_f = fa->fs->tke; fb_f = fb->fs->tke; fc_f = fc->fs->tke; fd_f = fd->fs->tke;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dtkedx = x1[1];
                vtx->dtkedy = x1[2];
                vtx->dtkedz = x1[3];
	ca_f = ca->fs->omega; cb_f = cb->fs->omega;
	        fa_f = fa->fs->omega; fb_f = fb->fs->omega; fc_f = fc->fs->omega; fd_f = fd->fs->omega;
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->domegadx = x1[1];
                vtx->domegady = x1[2];
                vtx->domegadz = x1[3];

        for ( size_t itm=0; itm<ntm; ++itm ) {
            ca_f = ca->fs->gas->T[itm]; cb_f = cb->fs->gas->T[itm];
	        fa_f = fa->fs->gas->T[itm]; fb_f = fb->fs->gas->T[itm]; fc_f = fc->fs->gas->T[itm]; fd_f = fd->fs->gas->T[itm];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dTdx[itm] = x1[1];
                vtx->dTdy[itm] = x1[2];
                vtx->dTdz[itm] = x1[3];
        }
        for( size_t isp = 0; isp < nsp; ++isp ) {
            ca_f = ca->fs->gas->massf[isp]; cb_f = cb->fs->gas->massf[isp];
	        fa_f = fa->fs->gas->massf[isp]; fb_f = fb->fs->gas->massf[isp]; fc_f = fc->fs->gas->massf[isp]; fd_f = fd->fs->gas->massf[isp];
                Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;
                A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->dfdx[isp] = x1[1];
                vtx->dfdy[isp] = x1[2];
                vtx->dfdz[isp] = x1[3];
        }
    }
    return SUCCESS;
} // end viscous_derivatives_edge_3D()

