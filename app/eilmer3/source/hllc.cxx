/** \file hllc.cxx
 * \ingroup eilmer3
 * \brief HLLC flux calculator.
 * \Author: Xin Kang
 * 11-Sep-2015
 * Reference: Riemann Solvers and Numerical Methods for Fluid Dynamics:
 * A practical Introduction, Toro et al.
 */

#include <stdio.h>
#include <math.h>
#include <numeric>
#include "../../../lib/util/source/useful.h"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"

/** \brief Compute the fluxes across an interface. */
int hllc(FlowState &Lft, FlowState &Rght, FV_Interface &IFace)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();

    double rL, rR;
    double pL, pR;
    double uL, uR;
    double vL, vR;
    double wL, wR;
    double aL, aR;
    double ieL, ieR;
    double keL, keR;
    double eL, eR;
    double HL, HR;

    // Step I: pressure estimation
    rL = Lft.gas->rho;
    pL = Lft.gas->p;
    uL = Lft.vel.x;
    vL = Lft.vel.y;
    wL = Lft.vel.z;
    aL = Lft.gas->a;
    ieL = accumulate(Lft.gas->e.begin(), Lft.gas->e.end(), 0.0);
    keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    eL = ieL + keL;
    HL = ieL + pL/rL + keL;

    rR = Rght.gas->rho;
    pR = Rght.gas->p;
    uR = Rght.vel.x;
    vR = Rght.vel.y;
    wR = Rght.vel.z;
    aR = Rght.gas->a;
    ieR = accumulate(Rght.gas->e.begin(), Rght.gas->e.end(), 0.0);
    keR = 0.5 * (uR * uR + vR * vR + wR * wR);
    eR = ieR + keR;
    HR = ieR + pR/rR + keR;

    double rho_mean;
    double a_mean;
    double p_pvrs, p_star;
    rho_mean = 0.5*(rL+rR);
    a_mean = 0.5*(aL+aR);
    p_pvrs = 0.5*(pL+pR)-0.5*(uR-uL)*rho_mean*a_mean;
    p_star = max(0.0, p_pvrs);


    // Step II: Wave speed estimation
    // Derive the gas "constants" from the local conditions.
    double cvL, RgasL;
    cvL = gmodel->Cv(*Lft.gas, IFace.status);
    RgasL = gmodel->R(*Lft.gas, IFace.status);

    double cvR, RgasR;
    cvR = gmodel->Cv(*Rght.gas, IFace.status);
    RgasR = gmodel->R(*Rght.gas, IFace.status);

    double rLsqrt, rRsqrt, alpha;
    rLsqrt = sqrt(rL);
    rRsqrt = sqrt(rR);
    alpha = rLsqrt / (rLsqrt + rRsqrt);

    double cv, Rgas;
    cv = alpha * cvL + (1.0 - alpha) * cvR;
    Rgas = alpha * RgasL + (1.0 - alpha) * RgasR;

    double cp, gam;
    cp = cv + Rgas;
    gam = cp / cv;   

    double SL, SR;
    double S_star;
    double qL, qR;
    // Left wave speed SL
    if (p_star <= pL) {
        qL = 1.0;
    } else {
        qL = sqrt(1.0+(gam+1.0)/2.0/gam*(p_star/pL-1.0));
    }
    SL = uL - aL*qL;

    // Right wave speed SR
    if (p_star <= pR) {
        qR = 1.0;
    } else {
        qR = sqrt(1.0+(gam+1.0)/2.0/gam*(p_star/pR-1.0));
    }
    SR = uR + aR*qR;

    // Intermediate wave speed S_star
    S_star = rL*(SL-uL) - rR*(SR-uR);
    S_star = (pR-pL+rL*uL*(SL-uL)-rR*uR*(SR-uR))/S_star;

   
    // Step III: HLLC flux

    ConservedQuantities &F = *(IFace.F);

    double U_factor_L, U_factor_R;
    U_factor_L = rL*(SL-uL)/(SL-S_star);
    U_factor_R = rR*(SR-uR)/(SR-S_star);

    // Mass flux
    double fmsL, fmsR;
    fmsL = rL*uL;
    fmsR = rR*uR;

    double fmsL_star, fmsR_star;
    double rL_star, rR_star;
    rL_star = U_factor_L;
    rR_star = U_factor_R;
    fmsL_star = fmsL + SL*(rL_star-rL);
    fmsR_star = fmsR + SR*(rR_star-rR);

    // Momentum flux
    // x direction
    double fmomxL, fmomxR;
    fmomxL = fmsL*uL + pL;
    fmomxR = fmsR*uR + pR;

    double fmomxL_star, fmomxR_star;
    double ruL_star, ruR_star;
    ruL_star = U_factor_L*S_star;
    ruR_star = U_factor_R*S_star;
    fmomxL_star = fmomxL + SL*(ruL_star-rL*uL);
    fmomxR_star = fmomxR + SR*(ruR_star-rR*uR);

    // y direction
    double fmomyL, fmomyR;
    fmomyL = fmsL*vL;
    fmomyR = fmsR*vR;

    double fmomyL_star, fmomyR_star;
    double rvL_star, rvR_star;
    rvL_star = U_factor_L*vL;
    rvR_star = U_factor_R*vR;
    fmomyL_star = fmomyL + SL*(rvL_star-rL*vL);
    fmomyR_star = fmomyR + SR*(rvR_star-rR*vR);

    // z direction
    double fmomzL, fmomzR;
    fmomzL = fmsL*wL;
    fmomzR = fmsR*wR;

    double fmomzL_star, fmomzR_star;
    double rwL_star, rwR_star;
    rwL_star = U_factor_L*wL;
    rwR_star = U_factor_R*wR;
    fmomzL_star = fmomzL + SL*(rwL_star-rL*wL);
    fmomzR_star = fmomzR + SR*(rwR_star-rR*wR);

    // Energy flux
    double fenergyL, fenergyR;
    fenergyL = fmsL*HL;
    fenergyR = fmsR*HR;

    double fenergyL_star, fenergyR_star;
    double reL_star, reR_star;
    reL_star = U_factor_L*(eL+(S_star-uL)*(S_star+pL/rL/(SL-uL)));
    reR_star = U_factor_R*(eR+(S_star-uR)*(S_star+pR/rR/(SR-uR)));
    fenergyL_star = fenergyL + SL*(reL_star-rL*eL);
    fenergyR_star = fenergyR + SR*(reR_star-rR*eR);

    // Species mass fluxes 
    vector<double> fmassfL;
    fmassfL.resize(nsp, 0.0);
    for ( int isp = 0; isp < nsp; ++isp ) {
	 fmassfL[isp] = fmsL*Lft.gas->massf[isp];
    }    
    vector<double> fmassfR;
    fmassfR.resize(nsp, 0.0);
    for ( int isp = 0; isp < nsp; ++isp ) {
	 fmassfR[isp] = fmsR*Rght.gas->massf[isp];
    }

    vector<double> fmassfL_star, fmassfR_star;
    fmassfL_star.resize(nsp, 0.0);
    fmassfR_star.resize(nsp, 0.0);
    vector<double> rmassfL_star, rmassfR_star;
    rmassfL_star.resize(nsp, 0.0);
    rmassfR_star.resize(nsp, 0.0);
    for ( int isp = 0; isp < nsp; ++isp ) {
    rmassfL_star[isp] = U_factor_L*Lft.gas->massf[isp];
    }
    for ( int isp = 0; isp < nsp; ++isp ) {
    rmassfR_star[isp] = U_factor_R*Rght.gas->massf[isp];
    }

    for ( int isp = 0; isp < nsp; ++isp ) {
    fmassfL_star[isp] = fmassfL[isp] + SL*(rmassfL_star[isp]-rL*Lft.gas->massf[isp]);
    }
    for ( int isp = 0; isp < nsp; ++isp ) {
    fmassfR_star[isp] = fmassfR[isp] + SR*(rmassfR_star[isp]-rR*Rght.gas->massf[isp]);
    }

    if(SL >= 0.0) {
       F.mass = fmsL;
       F.momentum.x = fmomxL;
       F.momentum.y = fmomyL;
       F.momentum.z = fmomzL;
       F.total_energy = fenergyL;
       for ( int isp = 0; isp < nsp; ++isp ) {
	    F.massf[isp] = fmassfL[isp];
       }
    }
    else if ( (SL < 0.0) && (S_star >= 0.0) ) {
       F.mass = fmsL_star;
       F.momentum.x = fmomxL_star;
       F.momentum.y = fmomyL_star;
       F.momentum.z = fmomzL_star;
       F.total_energy = fenergyL_star;
       for ( int isp = 0; isp < nsp; ++isp ) {
	    F.massf[isp] = fmassfL_star[isp];
       }
    }
    else if ( (S_star < 0.0) && (SR >= 0.0) ) {
       F.mass = fmsR_star;
       F.momentum.x = fmomxR_star;
       F.momentum.y = fmomyR_star;
       F.momentum.z = fmomzR_star;
       F.total_energy = fenergyR_star;
       for ( int isp = 0; isp < nsp; ++isp ) {
	    F.massf[isp] = fmassfR_star[isp];
       }
    }
    else {
       F.mass = fmsR;
       F.momentum.x = fmomxR;
       F.momentum.y = fmomyR;
       F.momentum.z = fmomzR;
       F.total_energy = fenergyR;
       for ( int isp = 0; isp < nsp; ++isp ) {
	    F.massf[isp] = fmassfR[isp];
       }
    }

    /* 
     * Individual energies 
     */
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if (F.mass >= 0.0) {
	/* Wind is blowing from the left */
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = F.mass * Lft.gas->e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += F.mass * Lft.gas->p_e / Lft.gas->rho;
    } else {
	/* Wind is blowing from the right */
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = F.mass * Rght.gas->e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += F.mass * Rght.gas->p_e / Rght.gas->rho;
    }

    return SUCCESS;
}   /* end of hllc() */

/*--------------------------------------------------------------*/


