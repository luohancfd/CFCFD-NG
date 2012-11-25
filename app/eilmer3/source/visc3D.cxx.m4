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
define(`for_each_scalar_element', `$1(vel.x,vel.x,u)
	$1(vel.y,vel.y,v)
	$1(vel.z,vel.z,w)
	$1(tke,tke,tke)
	$1(omega,omega,omega)
')dnl

define(`for_each_scalar_derivative', `$1(dudx)
       $1(dudy)
       $1(dudz)
       $1(dvdx)
       $1(dvdy)
       $1(dvdz)
       $1(dwdx)
       $1(dwdy)
       $1(dwdz)
       $1(dtkedx)
       $1(dtkedy)
       $1(dtkedz)
       $1(domegadx)
       $1(domegady)
       $1(domegadz)
')dnl


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
    int i, j, k;
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

    int nsp = gmodel->get_number_of_species();
    if( dfdx.size() == 0 ) {
	dfdx.resize(nsp); 
	dfdy.resize(nsp);
	dfdz.resize(nsp);
	jx.resize(nsp);
	jy.resize(nsp);
	jz.resize(nsp);
    }
    
    int ntm = gmodel->get_number_of_modes();
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
		    define(`select_upwind_scalar_value', 
                           `if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { $1 = Vtx1->$1; } else { $1 = Vtx4->$1; }
                           } else {
		               if ( vt2dp >= 0.0 ) { $1 = Vtx2->$1; } else { $1 = Vtx3->$1; }
                           }')dnl
		    define(`select_upwind_array_element', 
                           `if ( vt1dp >= 0.0 ) {
		               if ( vt2dp >= 0.0 ) { $1[$2] = Vtx1->$1[$2]; } else { $1[$2] = Vtx4->$1[$2]; }
                           } else {
		               if ( vt2dp >= 0.0 ) { $1[$2] = Vtx2->$1[$2]; } else { $1[$2] = Vtx3->$1[$2]; }
                           }')dnl
		    for_each_scalar_derivative(`select_upwind_scalar_value')
		    for ( int itm=0; itm<ntm; ++itm ) {
                        select_upwind_array_element(dTdx,itm)
                        select_upwind_array_element(dTdy,itm)
                        select_upwind_array_element(dTdz,itm)
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( int isp = 0; isp < nsp; ++isp ) {
                            select_upwind_array_element(dfdx,isp)
                            select_upwind_array_element(dfdy,isp)
                            select_upwind_array_element(dfdz,isp)
		        }
                    }
                } else {
                    // Symmetric average.
		    define(`average_scalar_over_face',
                           `$1 = 0.25*(Vtx1->$1+Vtx2->$1+Vtx3->$1+Vtx4->$1);')dnl
		    define(`average_array_element_over_face', 
                           `$1[$2] = 0.25*(Vtx1->$1[$2]+Vtx2->$1[$2]+Vtx3->$1[$2]+Vtx4->$1[$2]);')dnl
		    for_each_scalar_derivative(`average_scalar_over_face')
		    for ( int itm=0; itm<ntm; ++itm ) {
                        average_array_element_over_face(dTdx,itm)
                        average_array_element_over_face(dTdy,itm)
                        average_array_element_over_face(dTdz,itm)
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( int isp = 0; isp < nsp; ++isp ) {
                            average_array_element_over_face(dfdx,isp)
                            average_array_element_over_face(dfdy,isp)
                            average_array_element_over_face(dfdz,isp)
		        }
                    }
                }		    
                k_eff[0] = viscous_factor * (fs.gas->k[0] + fs.k_t);
		for ( int itm=1; itm<ntm; ++itm ) {
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
		    for( int isp = 0; isp < nsp; ++isp ) {
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
	        for ( int itm=1; itm<ntm; ++itm ) {
		    qx[itm] = k_eff[itm] * dTdx[itm];
		    qy[itm] = k_eff[itm] * dTdy[itm];
		    qz[itm] = k_eff[itm] * dTdz[itm];
		    qx[0] += qx[itm];
		    qy[0] += qy[itm];
		    qz[0] += qz[itm];
	        }
	        if( get_diffusion_flag() == 1 ) {
		    for( int isp = 0; isp < nsp; ++isp ) {
		    	double h = gmodel->enthalpy(*(fs.gas), isp);
		        qx[0] -= jx[isp] * h;
		        qy[0] -= jy[isp] * h;
		        qz[0] -= jz[isp] * h;
		        for ( int itm=1; itm<ntm; ++itm ) {
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
	  	    for( int isp = 0; isp < nsp; ++isp ) {
		        F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
		    }
	        }
	        // Modal energy flux (skipping first mode as this is handled by total energy)
	        for ( int itm=1; itm<ntm; ++itm ) {
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
		    for_each_scalar_derivative(`select_upwind_scalar_value')
		    for ( int itm=0; itm<ntm; ++itm ) {
                        select_upwind_array_element(dTdx,itm)
                        select_upwind_array_element(dTdy,itm)
                        select_upwind_array_element(dTdz,itm)
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( int isp = 0; isp < nsp; ++isp ) {
                            select_upwind_array_element(dfdx,isp)
                            select_upwind_array_element(dfdy,isp)
                            select_upwind_array_element(dfdz,isp)
		        }
                    }
                } else {
                    // Symmetric average.
 		    for_each_scalar_derivative(`average_scalar_over_face')
		    for ( int itm=0; itm<ntm; ++itm ) {
                        average_array_element_over_face(dTdx,itm)
                        average_array_element_over_face(dTdy,itm)
                        average_array_element_over_face(dTdz,itm)
                    }
 	            if( get_diffusion_flag() == 1 ) {
                        // derivatives needed for diffusion model, below
		        for( int isp = 0; isp < nsp; ++isp ) {
                            average_array_element_over_face(dfdx,isp)
                            average_array_element_over_face(dfdy,isp)
                            average_array_element_over_face(dfdz,isp)
		        }
                    }
                }
                k_eff[0] = viscous_factor * (fs.gas->k[0] + fs.k_t);
		for ( int itm=1; itm<ntm; ++itm ) {
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
		    for( int isp = 0; isp < nsp; ++isp ) {
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
	        for ( int itm=1; itm<ntm; ++itm ) {
		    qx[itm] = k_eff[itm] * dTdx[itm];
		    qy[itm] = k_eff[itm] * dTdy[itm];
		    qz[itm] = k_eff[itm] * dTdz[itm];
		    qx[0] += qx[itm];
		    qy[0] += qy[itm];
		    qz[0] += qz[itm];
	        }
	        if( get_diffusion_flag() == 1 ) {
		    for( int isp = 0; isp < nsp; ++isp ) {
		    	double h = gmodel->enthalpy(*(fs.gas), isp);
		        qx[0] -= jx[isp] * h;
		        qy[0] -= jy[isp] * h;
		        qz[0] -= jz[isp] * h;
		        for ( int itm=1; itm<ntm; ++itm ) {
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
	  	    for( int isp = 0; isp < nsp; ++isp ) {
		        F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
		    }
	        }
	        // Modal energy flux (skipping first mode as this is handled by total energy)
	        for ( int itm=1; itm<ntm; ++itm ) {
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
		    for_each_scalar_derivative(`select_upwind_scalar_value')
		    for ( int itm=0; itm<ntm; ++itm ) {
                        select_upwind_array_element(dTdx,itm)
                        select_upwind_array_element(dTdy,itm)
                        select_upwind_array_element(dTdz,itm)
                    }
	            if( get_diffusion_flag() == 1 ) {
                        // Needed for diffusion model, below.
		        for( int isp = 0; isp < nsp; ++isp ) {
                            select_upwind_array_element(dfdx,isp)
                            select_upwind_array_element(dfdy,isp)
                            select_upwind_array_element(dfdz,isp)
		        }
                    }
                } else {
                    // Symmetric average.
		    for_each_scalar_derivative(`average_scalar_over_face')
		    for ( int itm=0; itm<ntm; ++itm ) {
                        average_array_element_over_face(dTdx,itm)
                        average_array_element_over_face(dTdy,itm)
                        average_array_element_over_face(dTdz,itm)
                    }
 	            if( get_diffusion_flag() == 1 ) {
                        // derivatives needed for diffusion model, below
		        for( int isp = 0; isp < nsp; ++isp ) {
                            average_array_element_over_face(dfdx,isp)
                            average_array_element_over_face(dfdy,isp)
                            average_array_element_over_face(dfdz,isp)
		        }
                    }
                }
                k_eff[0] = viscous_factor * (fs.gas->k[0] + fs.k_t);
		for ( int itm=1; itm<ntm; ++itm ) {
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
		    for( int isp = 0; isp < nsp; ++isp ) {
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
	        for ( int itm=1; itm<ntm; ++itm ) {
		    qx[itm] = k_eff[itm] * dTdx[itm];
		    qy[itm] = k_eff[itm] * dTdy[itm];
		    qz[itm] = k_eff[itm] * dTdz[itm];
		    qx[0] += qx[itm];
		    qy[0] += qy[itm];
		    qz[0] += qz[itm];
	        }
	        if( get_diffusion_flag() == 1 ) {
		    for( int isp = 0; isp < nsp; ++isp ) {
		    	double h = gmodel->enthalpy(*(fs.gas), isp);
		        qx[0] -= jx[isp] * h;
		        qy[0] -= jy[isp] * h;
		        qz[0] -= jz[isp] * h;
		        for ( int itm=1; itm<ntm; ++itm ) {
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
	  	    for( int isp = 0; isp < nsp; ++isp ) {
		        F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
		    }
	        }
	        // Modal energy flux (skipping first mode as this is handled by total energy)
	        for ( int itm=1; itm<ntm; ++itm ) {
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
    int i, j, k;
    double q_e, q_w, q_n, q_s, q_top, q_bottom;
    double vol_inv;
    FV_Vertex *sec_ctr;
    // The ABC... notation is from Ian Johnston, and was used by Andrew Denman.
    // We'll keep using it to minimize the change from Elmer2 to Elmer3.
    FV_Cell *cA, *cB, *cC, *cD, *cE, *cF, *cG, *cH;
    FV_Interface *fA, *fB, *fC, *fD, *fE, *fF, *fG, *fH;
    FV_Interface *north, *east, *south, *west, *bottom, *top;

    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int ntm = gmodel->get_number_of_modes();

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
		define(`get_secondary_cell_interfaces', 
                       `// Secondary-cell interfaces
		       north = A->get_sifj(i,j+1,k);
		       east = A->get_sifi(i+1,j,k);
		       south = A->get_sifj(i,j,k);
		       west = A->get_sifi(i,j,k);
		       bottom = A->get_sifk(i,j,k);
		       top = A->get_sifk(i,j,k+1);')
		get_secondary_cell_interfaces
                define(`div_integral_expr',
                       `vol_inv * (q_e*east->area*east->n.$1 - q_w*west->area*west->n.$1 +
		        q_n*north->area*north->n.$1 - q_s*south->area*south->n.$1 +
		        q_top*top->area*top->n.$1 - q_bottom*bottom->area*bottom->n.$1)')dnl
		define(`div_theorem_core_scalar_element', 
                       `// Average property value on each face.
		       q_e = 0.25 * (cB->fs->$1 + cC->fs->$1 + cF->fs->$1 + cG->fs->$1);
		       q_w = 0.25 * (cA->fs->$1 + cD->fs->$1 + cH->fs->$1 + cE->fs->$1);
		       q_n = 0.25 * (cE->fs->$1 + cF->fs->$1 + cG->fs->$1 + cH->fs->$1);
		       q_s = 0.25 * (cA->fs->$1 + cB->fs->$1 + cC->fs->$1 + cD->fs->$1);
		       q_top = 0.25 * (cA->fs->$1 + cB->fs->$1 + cE->fs->$1 + cF->fs->$1);
		       q_bottom = 0.25 * (cD->fs->$1 + cC->fs->$1 + cG->fs->$1 + cH->fs->$1);
		       // Apply the divergence theorem.
		       sec_ctr->d$3dx = div_integral_expr(x);
		       sec_ctr->d$3dy = div_integral_expr(y);
		       sec_ctr->d$3dz = div_integral_expr(z);')
		for_each_scalar_element(`div_theorem_core_scalar_element')
		define(`div_theorem_core_array_element', 
                       `// Average property value on each face.
		       q_e = 0.25 * (cB->fs->$1[$4] + cC->fs->$1[$4] + cF->fs->$1[$4] + cG->fs->$1[$4]);
		       q_w = 0.25 * (cA->fs->$1[$4] + cD->fs->$1[$4] + cH->fs->$1[$4] + cE->fs->$1[$4]);
		       q_n = 0.25 * (cE->fs->$1[$4] + cF->fs->$1[$4] + cG->fs->$1[$4] + cH->fs->$1[$4]);
		       q_s = 0.25 * (cA->fs->$1[$4] + cB->fs->$1[$4] + cC->fs->$1[$4] + cD->fs->$1[$4]);
		       q_top = 0.25 * (cA->fs->$1[$4] + cB->fs->$1[$4] + cE->fs->$1[$4] + cF->fs->$1[$4]);
		       q_bottom = 0.25 * (cD->fs->$1[$4] + cC->fs->$1[$4] + cG->fs->$1[$4] + cH->fs->$1[$4]);
		       // Apply the divergence theorem.
		       sec_ctr->d$3dx[$4] = div_integral_expr(x);
		       sec_ctr->d$3dy[$4] = div_integral_expr(y);
		       sec_ctr->d$3dz[$4] = div_integral_expr(z);')
	        // Apply the divergence theorem to the primary temperature derivatives only.
                for ( int itm=0; itm<ntm; ++itm ) {
		    div_theorem_core_array_element(gas->T,gas->T,T,itm)
                }
  	        for( int isp = 0; isp < nsp; ++isp ) {
		    div_theorem_core_array_element(gas->massf,gas->massf,f,isp)
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
	    get_secondary_cell_interfaces
	    define(`div_theorem_east_scalar_element', 
                   `// Average property value on each face.
	           q_e = 0.25 * (fB->fs->$2 + fC->fs->$2 + fF->fs->$2 + fG->fs->$2);
	           q_w = 0.25 * (cA->fs->$1 + cD->fs->$1 + cH->fs->$1 + cE->fs->$1);
	           q_n = 0.25 * (cE->fs->$1 + fF->fs->$2 + fG->fs->$2 + cH->fs->$1);
	           q_s = 0.25 * (cA->fs->$1 + fB->fs->$2 + fC->fs->$2 + cD->fs->$1);
	           q_top = 0.25 * (cA->fs->$1 + fB->fs->$2 + cE->fs->$1 + fF->fs->$2);
	           q_bottom = 0.25 * (cD->fs->$1 + fC->fs->$2 + fG->fs->$2 + cH->fs->$1);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx = div_integral_expr(x);
	           sec_ctr->d$3dy = div_integral_expr(y);
	           sec_ctr->d$3dz = div_integral_expr(z);')
	    for_each_scalar_element(`div_theorem_east_scalar_element')
	    define(`div_theorem_east_array_element',
                   `// Average property value on each face.
	           q_e = 0.25 * (fB->fs->$2[$4] + fC->fs->$2[$4] + fF->fs->$2[$4] + fG->fs->$2[$4]);
	           q_w = 0.25 * (cA->fs->$1[$4] + cD->fs->$1[$4] + cH->fs->$1[$4] + cE->fs->$1[$4]);
	           q_n = 0.25 * (cE->fs->$1[$4] + fF->fs->$2[$4] + fG->fs->$2[$4] + cH->fs->$1[$4]);
	           q_s = 0.25 * (cA->fs->$1[$4] + fB->fs->$2[$4] + fC->fs->$2[$4] + cD->fs->$1[$4]);
	           q_top = 0.25 * (cA->fs->$1[$4] + fB->fs->$2[$4] + cE->fs->$1[$4] + fF->fs->$2[$4]);
	           q_bottom = 0.25 * (cD->fs->$1[$4] + fC->fs->$2[$4] + fG->fs->$2[$4] + cH->fs->$1[$4]);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx[$4] = div_integral_expr(x);
	           sec_ctr->d$3dy[$4] = div_integral_expr(y);
	           sec_ctr->d$3dz[$4] = div_integral_expr(z);')
            for ( int itm=0; itm<ntm; ++itm ) {
                div_theorem_east_array_element(gas->T,gas->T,T,itm)
            }
  	    for( int isp = 0; isp < nsp; ++isp ) {
		div_theorem_east_array_element(gas->massf,gas->massf,f,isp)
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
	    get_secondary_cell_interfaces
	    define(`div_theorem_west_scalar_element', 
                   `// Average property value on each face.
	           q_e = 0.25 * (cB->fs->$1 + cC->fs->$1 + cF->fs->$1 + cG->fs->$1);
	           q_w = 0.25 * (fA->fs->$2 + fD->fs->$2 + fH->fs->$2 + fE->fs->$2);
	           q_n = 0.25 * (fE->fs->$2 + cF->fs->$1 + cG->fs->$1 + fH->fs->$2);
	           q_s = 0.25 * (fA->fs->$2 + cB->fs->$1 + cC->fs->$1 + fD->fs->$2);
	           q_top = 0.25 * (fA->fs->$2 + cB->fs->$1 + fE->fs->$2 + cF->fs->$1);
	           q_bottom = 0.25 * (fD->fs->$2 + cC->fs->$1 + cG->fs->$1 + fH->fs->$2);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx = div_integral_expr(x);
	           sec_ctr->d$3dy = div_integral_expr(y);
	           sec_ctr->d$3dz = div_integral_expr(z);')
	    for_each_scalar_element(`div_theorem_west_scalar_element')
	    define(`div_theorem_west_array_element', `// Average property value on each face.
	           q_e = 0.25 * (cB->fs->$1[$4] + cC->fs->$1[$4] + cF->fs->$1[$4] + cG->fs->$1[$4]);
	           q_w = 0.25 * (fA->fs->$2[$4] + fD->fs->$2[$4] + fH->fs->$2[$4] + fE->fs->$2[$4]);
	           q_n = 0.25 * (fE->fs->$2[$4] + cF->fs->$1[$4] + cG->fs->$1[$4] + fH->fs->$2[$4]);
	           q_s = 0.25 * (fA->fs->$2[$4] + cB->fs->$1[$4] + cC->fs->$1[$4] + fD->fs->$2[$4]);
	           q_top = 0.25 * (fA->fs->$2[$4] + cB->fs->$1[$4] + fE->fs->$2[$4] + cF->fs->$1[$4]);
	           q_bottom = 0.25 * (fD->fs->$2[$4] + cC->fs->$1[$4] + cG->fs->$1[$4] + fH->fs->$2[$4]);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx[$4] = div_integral_expr(x);
	           sec_ctr->d$3dy[$4] = div_integral_expr(y);
	           sec_ctr->d$3dz[$4] = div_integral_expr(z);')
            for ( int itm=0; itm<ntm; ++itm ) {
	        div_theorem_west_array_element(gas->T,gas->T,T,itm)
            }
  	    for( int isp = 0; isp < nsp; ++isp ) {
		div_theorem_west_array_element(gas->massf,gas->massf,f,isp)
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
	    get_secondary_cell_interfaces
	    define(`div_theorem_north_scalar_element', 
                   `// Average property value on each face.
	           q_e = 0.25 * (cB->fs->$1 + cC->fs->$1 + fF->fs->$2 + fG->fs->$2);
	           q_w = 0.25 * (cA->fs->$1 + cD->fs->$1 + fH->fs->$2 + fE->fs->$2);
	           q_n = 0.25 * (fE->fs->$2 + fF->fs->$2 + fG->fs->$2 + fH->fs->$2);
	           q_s = 0.25 * (cA->fs->$1 + cB->fs->$1 + cC->fs->$1 + cD->fs->$1);
	           q_top = 0.25 * (cA->fs->$1 + cB->fs->$1 + fE->fs->$2 + fF->fs->$2);
	           q_bottom = 0.25 * (cD->fs->$1 + cC->fs->$1 + fG->fs->$2 + fH->fs->$2);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx = div_integral_expr(x);
	           sec_ctr->d$3dy = div_integral_expr(y);
	           sec_ctr->d$3dz = div_integral_expr(z);')
	    for_each_scalar_element(`div_theorem_north_scalar_element')
	    define(`div_theorem_north_array_element',
                   `// Average property value on each face.
	           q_e = 0.25 * (cB->fs->$1[$4] + cC->fs->$1[$4] + fF->fs->$2[$4] + fG->fs->$2[$4]);
	           q_w = 0.25 * (cA->fs->$1[$4] + cD->fs->$1[$4] + fH->fs->$2[$4] + fE->fs->$2[$4]);
	           q_n = 0.25 * (fE->fs->$2[$4] + fF->fs->$2[$4] + fG->fs->$2[$4] + fH->fs->$2[$4]);
	           q_s = 0.25 * (cA->fs->$1[$4] + cB->fs->$1[$4] + cC->fs->$1[$4] + cD->fs->$1[$4]);
	           q_top = 0.25 * (cA->fs->$1[$4] + cB->fs->$1[$4] + fE->fs->$2[$4] + fF->fs->$2[$4]);
	           q_bottom = 0.25 * (cD->fs->$1[$4] + cC->fs->$1[$4] + fG->fs->$2[$4] + fH->fs->$2[$4]);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx[$4] = div_integral_expr(x);
	           sec_ctr->d$3dy[$4] = div_integral_expr(y);
	           sec_ctr->d$3dz[$4] = div_integral_expr(z);')
            for ( int itm=0; itm<ntm; ++itm ) {
	        div_theorem_north_array_element(gas->T,gas->T,T,itm)
            }
  	    for( int isp = 0; isp < nsp; ++isp ) {
		div_theorem_north_array_element(gas->massf,gas->massf,f,isp)
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
	    get_secondary_cell_interfaces
	    define(`div_theorem_south_scalar_element',
                   `// Average property value on each face.
	           q_e = 0.25 * (fB->fs->$2 + fC->fs->$2 + cF->fs->$1 + cG->fs->$1);
	           q_w = 0.25 * (fA->fs->$2 + fD->fs->$2 + cH->fs->$1 + cE->fs->$1);
	           q_n = 0.25 * (cE->fs->$1 + cF->fs->$1 + cG->fs->$1 + cH->fs->$1);
	           q_s = 0.25 * (fA->fs->$2 + fB->fs->$2 + fC->fs->$2 + fD->fs->$2);
	           q_top = 0.25 * (fA->fs->$2 + fB->fs->$2 + cE->fs->$1 + cF->fs->$1);
	           q_bottom = 0.25 * (fD->fs->$2 + fC->fs->$2 + cG->fs->$1 + cH->fs->$1);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx = div_integral_expr(x);
	           sec_ctr->d$3dy = div_integral_expr(y);
	           sec_ctr->d$3dz = div_integral_expr(z);')
	    for_each_scalar_element(`div_theorem_south_scalar_element')
	    define(`div_theorem_south_array_element', 
                   `// Average property value on each face.
	           q_e = 0.25 * (fB->fs->$2[$4] + fC->fs->$2[$4] + cF->fs->$1[$4] + cG->fs->$1[$4]);
	           q_w = 0.25 * (fA->fs->$2[$4] + fD->fs->$2[$4] + cH->fs->$1[$4] + cE->fs->$1[$4]);
	           q_n = 0.25 * (cE->fs->$1[$4] + cF->fs->$1[$4] + cG->fs->$1[$4] + cH->fs->$1[$4]);
	           q_s = 0.25 * (fA->fs->$2[$4] + fB->fs->$2[$4] + fC->fs->$2[$4] + fD->fs->$2[$4]);
	           q_top = 0.25 * (fA->fs->$2[$4] + fB->fs->$2[$4] + cE->fs->$1[$4] + cF->fs->$1[$4]);
	           q_bottom = 0.25 * (fD->fs->$2[$4] + fC->fs->$2[$4] + cG->fs->$1[$4] + cH->fs->$1[$4]);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx[$4] = div_integral_expr(x);
	           sec_ctr->d$3dy[$4] = div_integral_expr(y);
	           sec_ctr->d$3dz[$4] = div_integral_expr(z);')
            for ( int itm=0; itm<ntm; ++itm ) {
	        div_theorem_south_array_element(gas->T,gas->T,T,itm)
            }
  	    for( int isp = 0; isp < nsp; ++isp ) {
		div_theorem_south_array_element(gas->massf,gas->massf,f,isp)
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
	    get_secondary_cell_interfaces
	    define(`div_theorem_top_scalar_element',
                   `// Average property value on each face.
	           q_e = 0.25 * (fB->fs->$2 + cC->fs->$1 + fF->fs->$2 + cG->fs->$1);
	           q_w = 0.25 * (fA->fs->$2 + cD->fs->$1 + cH->fs->$1 + fE->fs->$2);
	           q_n = 0.25 * (fE->fs->$2 + fF->fs->$2 + cG->fs->$1 + cH->fs->$1);
	           q_s = 0.25 * (fA->fs->$2 + fB->fs->$2 + cC->fs->$1 + cD->fs->$1);
	           q_top = 0.25 * (fA->fs->$2 + fB->fs->$2 + fE->fs->$2 + fF->fs->$2);
	           q_bottom = 0.25 * (cD->fs->$1 + cC->fs->$1 + cG->fs->$1 + cH->fs->$1);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx = div_integral_expr(x);
	           sec_ctr->d$3dy = div_integral_expr(y);
	           sec_ctr->d$3dz = div_integral_expr(z);')
	    for_each_scalar_element(`div_theorem_top_scalar_element')
	    define(`div_theorem_top_array_element', 
                   `// Average property value on each face.
	           q_e = 0.25 * (fB->fs->$2[$4] + cC->fs->$1[$4] + fF->fs->$2[$4] + cG->fs->$1[$4]);
	           q_w = 0.25 * (fA->fs->$2[$4] + cD->fs->$1[$4] + cH->fs->$1[$4] + fE->fs->$2[$4]);
	           q_n = 0.25 * (fE->fs->$2[$4] + fF->fs->$2[$4] + cG->fs->$1[$4] + cH->fs->$1[$4]);
	           q_s = 0.25 * (fA->fs->$2[$4] + fB->fs->$2[$4] + cC->fs->$1[$4] + cD->fs->$1[$4]);
	           q_top = 0.25 * (fA->fs->$2[$4] + fB->fs->$2[$4] + fE->fs->$2[$4] + fF->fs->$2[$4]);
	           q_bottom = 0.25 * (cD->fs->$1[$4] + cC->fs->$1[$4] + cG->fs->$1[$4] + cH->fs->$1[$4]);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx[$4] = div_integral_expr(x);
	           sec_ctr->d$3dy[$4] = div_integral_expr(y);
	           sec_ctr->d$3dz[$4] = div_integral_expr(z);')
            for ( int itm=0; itm<ntm; ++itm ) {
	        div_theorem_top_array_element(gas->T,gas->T,T,itm)
            }
  	    for( int isp = 0; isp < nsp; ++isp ) {
		div_theorem_top_array_element(gas->massf,gas->massf,f,isp)
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
	    get_secondary_cell_interfaces
	    define(`div_theorem_bottom_scalar_element', 
                   `// Average property value on each face.
	           q_e = 0.25 * (cB->fs->$1 + fC->fs->$2 + cF->fs->$1 + fG->fs->$2);
	           q_w = 0.25 * (cA->fs->$1 + fD->fs->$2 + fH->fs->$2 + cE->fs->$1);
	           q_n = 0.25 * (cE->fs->$1 + cF->fs->$1 + fG->fs->$2 + fH->fs->$2);
	           q_s = 0.25 * (cA->fs->$1 + cB->fs->$1 + fC->fs->$2 + fD->fs->$2);
	           q_top = 0.25 * (cA->fs->$1 + cB->fs->$1 + cE->fs->$1 + cF->fs->$1);
	           q_bottom = 0.25 * (fD->fs->$2 + fC->fs->$2 + fG->fs->$2 + fH->fs->$2);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx = div_integral_expr(x);
	           sec_ctr->d$3dy = div_integral_expr(y);
	           sec_ctr->d$3dz = div_integral_expr(z);')
	    for_each_scalar_element(`div_theorem_bottom_scalar_element')
	    define(`div_theorem_bottom_array_element',
                   `// Average property value on each face.
	           q_e = 0.25 * (cB->fs->$1[$4] + fC->fs->$2[$4] + cF->fs->$1[$4] + fG->fs->$2[$4]);
	           q_w = 0.25 * (cA->fs->$1[$4] + fD->fs->$2[$4] + fH->fs->$2[$4] + cE->fs->$1[$4]);
	           q_n = 0.25 * (cE->fs->$1[$4] + cF->fs->$1[$4] + fG->fs->$2[$4] + fH->fs->$2[$4]);
	           q_s = 0.25 * (cA->fs->$1[$4] + cB->fs->$1[$4] + fC->fs->$2[$4] + fD->fs->$2[$4]);
	           q_top = 0.25 * (cA->fs->$1[$4] + cB->fs->$1[$4] + cE->fs->$1[$4] + cF->fs->$1[$4]);
	           q_bottom = 0.25 * (fD->fs->$2[$4] + fC->fs->$2[$4] + fG->fs->$2[$4] + fH->fs->$2[$4]);
	           // Apply the divergence theorem.
	           sec_ctr->d$3dx[$4] = div_integral_expr(x);
	           sec_ctr->d$3dy[$4] = div_integral_expr(y);
	           sec_ctr->d$3dz[$4] = div_integral_expr(z);')
            for ( int itm=0; itm<ntm; ++itm ) {
	        div_theorem_bottom_array_element(gas->T,gas->T,T,itm)
            }
  	    for( int isp = 0; isp < nsp; ++isp ) {
		div_theorem_bottom_array_element(gas->massf,gas->massf,f,isp)
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
    int i, j, k;
    FV_Vertex *vtx;
    FV_Interface *a, *b, *d;
    FV_Cell *c;
    double xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd;
    double fa, fb, fc, fd, denom;

    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int ntm = gmodel->get_number_of_modes();

    // South-West-Bottom corner [0]
    i = bdp->imin; j = bdp->jmin; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j,k);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k);
    define(`denominator_for_corners', 
           `xa = a->pos.x; ya = a->pos.y; za = a->pos.z;
            xb = b->pos.x; yb = b->pos.y; zb = b->pos.z;
            xc = c->pos.x; yc = c->pos.y; zc = c->pos.z;
            xd = d->pos.x; yd = c->pos.y; zd = d->pos.z;
            denom = xa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+xb*(yc*zd-yd*zc)
                +ya*(xc*zd+xb*(zc-zd)-xd*zc+(xd-xc)*zb)+yb*(xd*zc-xc*zd)
                +(xc*yd-xd*yc)*zb+(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)*za;')dnl
    denominator_for_corners
    define(`numerator_x',
           `(fa*(yb*(zd-zc)-yc*zd+yd*zc+(yc-yd)*zb)+fb*(yc*zd-yd*zc)
             +ya*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+yb*(fd*zc-fc*zd)
             +(fc*yd-fd*yc)*zb+(fb*(yd-yc)-fc*yd+fd*yc+(fc-fd)*yb)*za)')dnl
    define(`numerator_y',
           `-(fa*(xb*(zd-zc)-xc*zd+xd*zc+(xc-xd)*zb)+fb*(xc*zd-xd*zc)
              +xa*(fc*zd+fb*(zc-zd)-fd*zc+(fd-fc)*zb)+xb*(fd*zc-fc*zd)
              +(fc*xd-fd*xc)*zb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*za)')dnl
    define(`numerator_z',
           `(fa*(xb*(yd-yc)-xc*yd+xd*yc+(xc-xd)*yb)+fb*(xc*yd-xd*yc)
             +xa*(fc*yd+fb*(yc-yd)-fd*yc+(fd-fc)*yb)+xb*(fd*yc-fc*yd)
             +(fc*xd-fd*xc)*yb+(fb*(xd-xc)-fc*xd+fd*xc+(fc-fd)*xb)*ya)')dnl
    define(`corner_derivatives_scalar_element', 
           `fa = a->fs->$2; fb = b->fs->$2; fc = c->fs->$1; fd = d->fs->$2;
           vtx->d$3dx = numerator_x / denom;
           vtx->d$3dy = numerator_y / denom;
           vtx->d$3dz = numerator_z / denom;')dnl
    for_each_scalar_element(`corner_derivatives_scalar_element')
    define(`corner_derivatives_array_element', 
           `fa = a->fs->$2[$4]; fb = b->fs->$2[$4]; fc = c->fs->$1[$4]; fd = d->fs->$2[$4];
           vtx->d$3dx[$4] = numerator_x / denom;
           vtx->d$3dy[$4] = numerator_y / denom;
           vtx->d$3dz[$4] = numerator_z / denom;')dnl
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    // South-East-Bottom corner [1]
    i = bdp->imax; j = bdp->jmin; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j,k);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k);
    denominator_for_corners
    for_each_scalar_element(`corner_derivatives_scalar_element')
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    // North-East-Bottom corner [2]
    i = bdp->imax; j = bdp->jmax; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j+1,k);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k);
    denominator_for_corners
    for_each_scalar_element(`corner_derivatives_scalar_element')
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    // North-West-Bottom corner [3]
    i = bdp->imin; j = bdp->jmax; k = bdp->kmin;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j+1,k);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k);
    denominator_for_corners
    for_each_scalar_element(`corner_derivatives_scalar_element')
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    // South-West-Top corner [4]
    i = bdp->imin; j = bdp->jmin; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j,k+1);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k+1);
    denominator_for_corners
    for_each_scalar_element(`corner_derivatives_scalar_element')
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    // South-East-Top corner [5]
    i = bdp->imax; j = bdp->jmin; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j,k+1);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j,k);
    d = bdp->get_ifk(i,j,k+1);
    denominator_for_corners
    for_each_scalar_element(`corner_derivatives_scalar_element')
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    // North-East-Top corner [6]
    i = bdp->imax; j = bdp->jmax; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i+1,j+1,k+1);
    a = bdp->get_ifi(i+1,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k+1);
    denominator_for_corners
    for_each_scalar_element(`corner_derivatives_scalar_element')
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    // North-West-Top corner [7]
    i = bdp->imin; j = bdp->jmax; k = bdp->kmax;
    c = bdp->get_cell(i,j,k);
    vtx = bdp->get_vtx(i,j+1,k+1);
    a = bdp->get_ifi(i,j,k);
    b = bdp->get_ifj(i,j+1,k);
    d = bdp->get_ifk(i,j,k+1);
    denominator_for_corners
    for_each_scalar_element(`corner_derivatives_scalar_element')
    for ( int itm=0; itm<ntm; ++itm ) {
        corner_derivatives_array_element(gas->T,gas->T,T,itm)
    }
    for( int isp = 0; isp < nsp; ++isp ) {
        corner_derivatives_array_element(gas->massf,gas->massf,f,isp)
    }

    return SUCCESS;
} // end viscous_derivatives_corners_3D()


/** \brief Edge derivatives from a least-squares fit of near-by data.
 *
 * See workbook pages dated 07-Jan-2010.
 */
int viscous_derivatives_edge_3D(Block *bdp)
{
    int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
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
    int nsp = gmodel->get_number_of_species();
    int ntm = gmodel->get_number_of_modes();

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
	define(`assemble_least_squares_matrix', 
               `ca_x = ca->pos.x; ca_y = ca->pos.y; ca_z = ca->pos.z;
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
                Szz = ca_z*ca_z + cb_z*cb_z + fa_z*fa_z + fb_z*fb_z + fc_z*fc_z + fd_z*fd_z;')dnl
        assemble_least_squares_matrix
        define(`get_rhs_scalar_elements',
               `ca_f = ca->fs->$1; cb_f = cb->fs->$1;
	        fa_f = fa->fs->$2; fb_f = fb->fs->$2; fc_f = fc->fs->$2; fd_f = fd->fs->$2;')dnl
        define(`get_rhs_array_elements',
               `ca_f = ca->fs->$1[$4]; cb_f = cb->fs->$1[$4];
	        fa_f = fa->fs->$2[$4]; fb_f = fb->fs->$2[$4]; fc_f = fc->fs->$2[$4]; fd_f = fd->fs->$2[$4];')dnl
	define(`assemble_least_squares_rhs', 
               `Sf = ca_f + cb_f + fa_f + fb_f + fc_f + fd_f;
	        Sfx = ca_f*ca_x + cb_f*cb_x + fa_f*fa_x + fb_f*fb_x + fc_f*fc_x + fd_f*fd_x;
	        Sfy = ca_f*ca_y + cb_f*cb_y + fa_f*fa_y + fb_f*fb_y + fc_f*fc_y + fd_f*fd_y;
	        Sfz = ca_f*ca_z + cb_f*cb_z + fa_f*fa_z + fb_f*fb_z + fc_f*fc_z + fd_f*fd_z;')dnl
        define(`solve_lsq_scalar_element',
               `A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->d$3dx = x1[1];
                vtx->d$3dy = x1[2];
                vtx->d$3dz = x1[3];')dnl 
        define(`solve_lsq_array_element',
               `A1.set(0,0,S1);  A1.set(0,1,Sx);  A1.set(0,2,Sy);  A1.set(0,3,Sz);
                A1.set(1,0,Sx);  A1.set(1,1,Sxx);  A1.set(1,2,Syx);  A1.set(1,3,Szx);
                A1.set(2,0,Sy);  A1.set(2,1,Syx);  A1.set(2,2,Syy);  A1.set(2,3,Szy);
                A1.set(3,0,Sz);  A1.set(3,1,Szx);  A1.set(3,2,Szy);  A1.set(3,3,Szz);
                B1[0] = Sf; B1[1] = Sfx; B1[2] = Sfy; B1[3] = Sfz;
                gaussian_elimination( A1, x1, B1 );
                vtx->d$3dx[$4] = x1[1];
                vtx->d$3dy[$4] = x1[2];
                vtx->d$3dz[$4] = x1[3];')dnl
        define(`assemble_least_squares_rhs_and_solve_scalar_element',
               `get_rhs_scalar_elements($1,$2)
                assemble_least_squares_rhs
                solve_lsq_scalar_element($1,$2,$3)')dnl
        define(`assemble_least_squares_rhs_and_solve_array_element',
               `get_rhs_array_elements($1,$2,$3,$4)
                assemble_least_squares_rhs
                solve_lsq_array_element($1,$2,$3,$4)')dnl
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
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
        assemble_least_squares_matrix
        for_each_scalar_element(`assemble_least_squares_rhs_and_solve_scalar_element')
        for ( int itm=0; itm<ntm; ++itm ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->T,gas->T,T,itm)
        }
        for( int isp = 0; isp < nsp; ++isp ) {
            assemble_least_squares_rhs_and_solve_array_element(gas->massf,gas->massf,f,isp)
        }
    }
    return SUCCESS;
} // end viscous_derivatives_edge_3D()

