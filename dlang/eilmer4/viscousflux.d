/**
 * viscousflux.d
 * Viscous-Flux calculation, where the fluxes are driven by molecular-transport effects.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-05-02: port essentials from Eilmer3 and refactor (a lot).
 */

module viscousflux;

import std.math;
import std.stdio;
import std.conv;
import geom;
import gas;
import flowstate;
import conservedquantities;
import fvcore;
import fvinterface;
import fvvertex;
import globalconfig;

class ViscousFluxData {
public:
    double[][] grad_vel;
    Vector3 grad_T, grad_tke, grad_omega;
    Vector3[] grad_f;
    double[] jx, jy, jz;
    size_t nsp;

    this() {
	nsp = GlobalConfig.gmodel.n_species;
	grad_f.length = nsp; 
	if ( GlobalConfig.diffusion ) {
	    jx.length = nsp;
	    jy.length = nsp;
	    jz.length = nsp;
	}
    }

    @nogc
    void viscous_flux_calc(ref FVInterface IFace)
    // Unified 2D and 3D viscous-flux calculation.
    // Note that the gradient values need to be in place before calling this procedure.
    {
	double viscous_factor = GlobalConfig.viscous_factor;
	FlowState fs = IFace.fs;
        double k_eff = viscous_factor * (fs.gas.k[0] + fs.k_t);
	double mu_eff =  viscous_factor * (fs.gas.mu + fs.mu_t);
	double lmbda = -2.0/3.0 * mu_eff;
	if ( GlobalConfig.diffusion ) {
	    // Apply a diffusion model
	    double D_t = 0.0;
	    if ( GlobalConfig.turbulence_model != TurbulenceModel.none ) {
		double Sc_t = GlobalConfig.turbulence_schmidt_number;
		D_t = fs.mu_t / (fs.gas.rho * Sc_t);
	    }
	    // [TODO] Rowan, calculate_diffusion_fluxes(fs.gas, D_t, grad_f, jx, jy, jz);
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		jx[isp] = 0.0;
		jy[isp] = 0.0;
		jz[isp] = 0.0;
	    }
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		jx[isp] *= viscous_factor;
		jy[isp] *= viscous_factor;
		jz[isp] *= viscous_factor;
	    }
	}
	double tau_xx = 0.0;
	double tau_yy = 0.0;
	double tau_zz = 0.0;
	double tau_xy = 0.0;
	double tau_xz = 0.0;
	double tau_yz = 0.0;
	if (GlobalConfig.dimensions == 3) {
	    double dudx = grad_vel[0][0];
	    double dudy = grad_vel[0][1];
	    double dudz = grad_vel[0][2];
	    double dvdx = grad_vel[1][0];
	    double dvdy = grad_vel[1][1];
	    double dvdz = grad_vel[1][2];
	    double dwdx = grad_vel[2][0];
	    double dwdy = grad_vel[1][1];
	    double dwdz = grad_vel[2][2];
	    // 3-dimensional planar stresses.
	    tau_xx = 2.0*mu_eff*dudx + lmbda*(dudx + dvdy + dwdz);
	    tau_yy = 2.0*mu_eff*dvdy + lmbda*(dudx + dvdy + dwdz);
	    tau_zz = 2.0*mu_eff*dwdz + lmbda*(dudx + dvdy + dwdz);
	    tau_xy = mu_eff * (dudy + dvdx);
	    tau_xz = mu_eff * (dudz + dwdx);
	    tau_yz = mu_eff * (dvdz + dwdy);
	} else {
	    // 2D
	    double dudx = grad_vel[0][0];
	    double dudy = grad_vel[0][1];
	    double dvdx = grad_vel[1][0];
	    double dvdy = grad_vel[1][1];
	    if (GlobalConfig.axisymmetric) {
		// Viscous stresses at the mid-point of the interface.
		// Axisymmetric terms no longer include the radial multiplier
		// as that has been absorbed into the interface area calculation.
		double ybar = IFace.Ybar;
                if (ybar > 1.0e-10) { // something very small for a cell height
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
	}
	// Thermal conductivity (NOTE: q is total energy flux)
	double qx = k_eff * grad_T.x;
	double qy = k_eff * grad_T.y;
	double qz = k_eff * grad_T.z;
	if ( GlobalConfig.diffusion ) {
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		double h = 0.0; // [TODO] Rowan, transport of species enthalpies?
		// double h = GlobalConfig.gmodel.enthalpy(fs.gas, isp);
		qx -= jx[isp] * h;
		qy -= jy[isp] * h;
		qz -= jz[isp] * h;
		// [TODO] Rowan, modal enthalpies ?
	    }
	}
	double tau_kx = 0.0;
	double tau_ky = 0.0;
	double tau_kz = 0.0;
	double tau_wx = 0.0;
	double tau_wy = 0.0;
	double tau_wz = 0.0;
	if ( GlobalConfig.turbulence_model == TurbulenceModel.k_omega ) {
	    // Turbulence contribution to the shear stresses.
	    tau_xx -= 0.66667 * fs.gas.rho * fs.tke;
	    tau_yy -= 0.66667 * fs.gas.rho * fs.tke;
	    if (GlobalConfig.dimensions == 3) { tau_zz -= 0.66667 * fs.gas.rho * fs.tke; }
	    // Turbulence contribution to heat transfer.
	    double sigma_star = 0.6;
	    double mu_effective = fs.gas.mu + sigma_star * fs.mu_t;
	    qx += mu_effective * grad_tke.x;
	    qy += mu_effective * grad_tke.y;
	    if (GlobalConfig.dimensions == 3) { qz += mu_effective * grad_tke.z; }
	    // Turbulence transport of the turbulence properties themselves.
	    tau_kx = mu_effective * grad_tke.x; 
	    tau_ky = mu_effective * grad_tke.y;
	    if (GlobalConfig.dimensions == 3) { tau_kz = mu_effective * grad_tke.z; }
	    double sigma = 0.5;
	    mu_effective = fs.gas.mu + sigma * fs.mu_t;
	    tau_wx = mu_effective * grad_omega.x; 
	    tau_wy = mu_effective * grad_omega.y; 
	    if (GlobalConfig.dimensions == 3) { tau_wz = mu_effective * grad_omega.z; } 
	}
	// Combine into fluxes: store as the dot product (F.n).
	ConservedQuantities F = IFace.F;
	double nx = IFace.n.x;
	double ny = IFace.n.y;
	double nz = IFace.n.z;
	// Mass flux -- NO CONTRIBUTION, unless there's diffusion (below)
	F.momentum.refx -= tau_xx*nx + tau_xy*ny + tau_xz*nz;
	F.momentum.refy -= tau_xy*nx + tau_yy*ny + tau_yz*nz;
	F.momentum.refz -= tau_xz*nx + tau_yz*ny + tau_zz*nz;
	F.total_energy -=
	    (tau_xx*fs.vel.x + tau_xy*fs.vel.y + tau_xz*fs.vel.z + qx)*nx +
	    (tau_xy*fs.vel.x + tau_yy*fs.vel.y + tau_yz*fs.vel.z + qy)*ny +
	    (tau_xz*fs.vel.x + tau_yz*fs.vel.y + tau_zz*fs.vel.z + qz)*nz;
	if (GlobalConfig.turbulence_model == TurbulenceModel.k_omega) {
	    F.tke -= tau_kx * nx + tau_ky * ny + tau_kz * nz;
	    F.omega -= tau_wx * nx + tau_wy * ny + tau_wz * nz;
	}
	if (GlobalConfig.diffusion) {
	    // Species mass flux
	    // [TODO] Rowan, what happens with user-defined flux?
	    for( size_t isp = 0; isp < nsp; ++isp ) {
		F.massf[isp] += jx[isp]*nx + jy[isp]*ny + jz[isp]*nz;
	    }
	}
	// [TODO] Rowan, Modal energy flux?
    } // end viscous_flux_calc()

    @nogc
    void average_vertex_values_2D(const FVVertex vtx1, const FVVertex vtx2)
    {
	grad_vel[0][0] = 0.5*(vtx1.grad_vel[0][0] + vtx1.grad_vel[0][0]); // du/dx
	grad_vel[0][1] = 0.5*(vtx1.grad_vel[0][1] + vtx2.grad_vel[0][1]); // du/dy
	grad_vel[0][2] = 0.0; // du/dz
	grad_vel[1][0] = 0.5*(vtx1.grad_vel[1][0] + vtx2.grad_vel[1][0]); // dv/dx
	grad_vel[1][1] = 0.5*(vtx1.grad_vel[1][1] + vtx2.grad_vel[1][1]); // dv/dy
	grad_vel[1][2] = 0.0; // dv/dz
	grad_vel[2][0] = 0.0; // dw/dx
	grad_vel[2][1] = 0.0; // dw/dy
	grad_vel[2][2] = 0.0; // dw/dz
	grad_tke.refx = 0.5*(vtx1.grad_tke.x + vtx2.grad_tke.x);
	grad_tke.refy = 0.5*(vtx1.grad_tke.y + vtx2.grad_tke.y);
	grad_tke.refz = 0.0;
	grad_omega.refx = 0.5*(vtx1.grad_omega.x + vtx2.grad_omega.x);
	grad_omega.refy = 0.5*(vtx1.grad_omega.y + vtx2.grad_omega.y);
	grad_omega.refz = 0.0;
	grad_T.refx = 0.5*(vtx1.grad_T.x + vtx2.grad_T.x);
	grad_T.refy = 0.5*(vtx1.grad_T.y + vtx2.grad_T.y);
	grad_T.refz = 0.0;
	for (size_t isp = 0; isp < nsp; ++isp) {
	    grad_f[isp].refx = 0.5*(vtx1.grad_f[isp].x + vtx2.grad_f[isp].x);
	    grad_f[isp].refy = 0.5*(vtx1.grad_f[isp].y + vtx2.grad_f[isp].y);
	    grad_f[isp].refz = 0.0;
	}
    } // end average_vertex_values_2D()

    @nogc
    void average_vertex_values_3D(const FVVertex vtx1, const FVVertex vtx2,
				  const FVVertex vtx3, const FVVertex vtx4)
    {
	grad_vel[0][0] = 0.25*(vtx1.grad_vel[0][0] + vtx2.grad_vel[0][0] +
			       vtx3.grad_vel[0][0] + vtx4.grad_vel[0][0]); // du/dx
	grad_vel[0][1] = 0.25*(vtx1.grad_vel[0][1] + vtx2.grad_vel[0][1] +
			       vtx3.grad_vel[0][1] + vtx4.grad_vel[0][1]); // du/dy
	grad_vel[0][2] = 0.25*(vtx1.grad_vel[0][2] + vtx2.grad_vel[0][2] +
			       vtx3.grad_vel[0][2] + vtx4.grad_vel[0][2]); // du/dz
	grad_vel[1][0] = 0.25*(vtx1.grad_vel[1][0] + vtx2.grad_vel[1][0] + 
			       vtx3.grad_vel[1][0] + vtx4.grad_vel[1][0]); // dv/dx
	grad_vel[1][1] = 0.25*(vtx1.grad_vel[1][1] + vtx2.grad_vel[1][1] + 
			       vtx3.grad_vel[1][1] + vtx4.grad_vel[1][1]); // dv/dy
	grad_vel[1][2] = 0.25*(vtx1.grad_vel[1][2] + vtx2.grad_vel[1][2] + 
			       vtx3.grad_vel[1][2] + vtx4.grad_vel[1][2]); // dv/dz
	grad_vel[2][0] = 0.25*(vtx1.grad_vel[1][0] + vtx2.grad_vel[1][0] + 
			       vtx3.grad_vel[1][0] + vtx4.grad_vel[1][0]); // dw/dx
	grad_vel[2][1] = 0.25*(vtx1.grad_vel[2][1] + vtx2.grad_vel[2][1] + 
			       vtx3.grad_vel[2][1] + vtx4.grad_vel[2][1]); // dw/dy
	grad_vel[2][2] = 0.25*(vtx1.grad_vel[2][2] + vtx2.grad_vel[2][2] + 
			       vtx3.grad_vel[2][2] + vtx4.grad_vel[2][2]); // dw/dz
	grad_tke.refx = 0.25*(vtx1.grad_tke.x + vtx2.grad_tke.x + 
			      vtx3.grad_tke.x + vtx4.grad_tke.x);
	grad_tke.refy = 0.25*(vtx1.grad_tke.y + vtx2.grad_tke.y +
			      vtx3.grad_tke.y + vtx4.grad_tke.y);
	grad_tke.refz = 0.25*(vtx1.grad_tke.z + vtx2.grad_tke.z +
			      vtx3.grad_tke.z + vtx4.grad_tke.z);
	grad_omega.refx = 0.25*(vtx1.grad_omega.x + vtx2.grad_omega.x +
				vtx3.grad_omega.x + vtx4.grad_omega.x);
	grad_omega.refy = 0.25*(vtx1.grad_omega.y + vtx2.grad_omega.y +
				vtx3.grad_omega.y + vtx4.grad_omega.y);
	grad_omega.refz = 0.25*(vtx1.grad_omega.z + vtx2.grad_omega.z +
				vtx3.grad_omega.z + vtx4.grad_omega.z);
	grad_T.refx = 0.25*(vtx1.grad_T.x + vtx2.grad_T.x +
			    vtx3.grad_T.x + vtx4.grad_T.x);
	grad_T.refy = 0.25*(vtx1.grad_T.y + vtx2.grad_T.y +
			    vtx3.grad_T.y + vtx4.grad_T.y);
	grad_T.refz = 0.25*(vtx1.grad_T.z + vtx2.grad_T.z +
			    vtx3.grad_T.z + vtx4.grad_T.z);
	for (size_t isp = 0; isp < nsp; ++isp) {
	    grad_f[isp].refx = 0.25*(vtx1.grad_f[isp].x + vtx2.grad_f[isp].x +
				     vtx3.grad_f[isp].x + vtx4.grad_f[isp].x);
	    grad_f[isp].refy = 0.25*(vtx1.grad_f[isp].y + vtx2.grad_f[isp].y +
				     vtx3.grad_f[isp].y + vtx4.grad_f[isp].y);
	    grad_f[isp].refz = 0.25*(vtx1.grad_f[isp].z + vtx2.grad_f[isp].z +
				     vtx3.grad_f[isp].z + vtx4.grad_f[isp].z);
	}
    } // end average_vertex_values_3D()

} // end class ViscousFluxData
