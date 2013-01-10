/** \file ausmdv.cxx
 * \ingroup eilmer3
 * \brief Wada and Liou's flux calculator.
 * 
 * Implemented from details in their AIAA paper 
 * with hints from Ian Johnston.
 * \verbatim
 * Y. Wada and M. -S. Liou (1994)
 * A flux splitting scheme with high-resolution and 
 * robustness for discontinuities.
 * AIAA-94-0083.
 * \endverbatim
 *
 * \author P. A. Jacobs
 *     Department of Mechanical Engineering
 *     The University of Queensland
 *
 *
 * \version
 * \verbatim
 * 26-Jun-97: Initial coding.
 * 08-Jul-97: Fix entropy fix (ironic, eh?)
 *            And finally debugged by Ian J.
 * 12-Oct-97: data structure for fluxes, array data
 * 15-Oct-97: vectorise for Cray
 * \endverbatim
 * 
 * \todo Really should get rid of the vector loop some time.
 */

#include <stdio.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"

/*------------------------------------------------------------*/

#define  K_SWITCH  10.0
#define  C_EFIX    0.125

/*------------------------------------------------------------*/

/** \brief Compute the fluxes across an interface. */
int ausmdv(FlowState &Lft, FlowState &Rght, FV_Interface &IFace)
{
    double rL, rR;
    double pL, pR;
    double uL, uR;
    double aL, aR;
    double vL, vR;
    double wL, wR;
    double HL, HR;
    double pLrL, pRrR;
    double ML, MR;
    double eL, eR;
    double keL, keR;
    double alphaL, alphaR, am;
    double pLplus, pRminus;
    double uLplus, uRminus;
    double duL, duR;

    double p_half, ru_half, ru2_half;
    double dp, s, ru2_AUSMV, ru2_AUSMD;

    int caseA, caseB;
    double d_ua;

    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();
    
    /*
     * Unpack the flow-state vectors for either side of the interface.
     * Store in work vectors, those quantities that will be neede later.
     */
    rL = Lft.gas->rho;
    pL = Lft.gas->p;
    pLrL = pL / rL;
    uL = Lft.vel.x - IFace.vel.x;
    vL = Lft.vel.y - IFace.vel.y;
    wL = Lft.vel.z - IFace.vel.z;
    eL = Lft.gas->e[0];
    aL = Lft.gas->a;
    keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    HL = eL + pLrL + keL;

    rR = Rght.gas->rho;
    pR = Rght.gas->p;
    pRrR = pR / rR;
    uR = Rght.vel.x - IFace.vel.x;
    vR = Rght.vel.y - IFace.vel.y;
    wR = Rght.vel.z - IFace.vel.z;
    eR = Rght.gas->e[0];
    aR = Rght.gas->a;
    keR = 0.5 * (uR * uR + vR * vR + wR * wR);
    HR = eR + pRrR + keR;

    /*
     * This is the main part of the flux calculator.
     */
    /*
     * Weighting parameters (eqn 32) for velocity splitting.
     */
    alphaL = 2.0 * pLrL / (pLrL + pRrR);
    alphaR = 2.0 * pRrR / (pLrL + pRrR);
    /*
     * Common sound speed (eqn 33) and Mach numbers.
     */
    am = MAXIMUM(aL, aR);
    ML = uL / am;
    MR = uR / am;
    /*
     * Left state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 30)
     */
    duL = 0.5 * (uL + FABS(uL));
    if (FABS(ML) <= 1.0) {
	pLplus = pL * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) * 0.25;
	uLplus =
	    alphaL * ((uL + am) * (uL + am) / (4.0 * am) - duL) +
	    duL;
    } else {
	pLplus = pL * duL / uL;
	uLplus = duL;
    }
    /*
     * Right state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 31)
     */
    duR = 0.5 * (uR - FABS(uR));
    if (FABS(MR) <= 1.0) {
	pRminus = pR * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) * 0.25;
	uRminus =
	    alphaR * (-(uR - am) * (uR - am) / (4.0 * am) - duR) +
	    duR;
    } else {
	pRminus = pR * duR / uR;
	uRminus = duR;
    }
    /*
     * Mass Flux (eqn 29)
     */
    ru_half = uLplus * rL + uRminus * rR;
    /*
     * Pressure flux (eqn 34)
     */
    p_half = pLplus + pRminus;
    /*
     * Momentum flux: normal direction
     *
     * Compute blending parameter s (eqn 37),
     * the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
     * and blend (eqn 36).
     */
    dp = pL - pR;
    dp = K_SWITCH * FABS(dp) / MINIMUM(pL, pR);
    s = 0.5 * MINIMUM(1.0, dp);

    ru2_AUSMV = uLplus * rL * uL + uRminus * rR * uR;
    ru2_AUSMD = 0.5 * (ru_half * (uL + uR) -
		       FABS(ru_half) * (uR - uL));

    ru2_half = (0.5 + s) * ru2_AUSMV + (0.5 - s) * ru2_AUSMD;

    /*
     * Assemble components of the flux vector.
     */
    ConservedQuantities &F = *(IFace.F);
    F.mass = ru_half;
    F.momentum.x = ru2_half + p_half;
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	F.momentum.y = ru_half * vL;
	F.momentum.z = ru_half * wL;
	F.total_energy = ru_half * HL;
	F.tke = ru_half * Lft.tke;
	F.omega = ru_half * Lft.omega;
    } else {
	/* Wind is blowing from the right */
	F.momentum.y = ru_half * vR;
	F.momentum.z = ru_half * wR;
	F.total_energy = ru_half * HR;
	F.tke = ru_half * Rght.tke;
	F.omega = ru_half * Rght.omega;
    }
    /* 
     * Species mass fluxes 
     */
    for ( int isp = 0; isp < nsp; ++isp ) {
	if (ru_half >= 0.0) {
	    /* Wind is blowing from the left */
	    F.massf[isp] = ru_half * Lft.gas->massf[isp];
	} else {
	    /* Wind is blowing from the right */
	    F.massf[isp] = ru_half * Rght.gas->massf[isp];
	}
    } /* isp loop */
    /* 
     * Individual energies 
     */
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = ru_half * Lft.gas->e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += ru_half * Lft.gas->p_e / Lft.gas->rho;
    } else {
	/* Wind is blowing from the right */
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = ru_half * Rght.gas->e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += ru_half * Rght.gas->p_e / Rght.gas->rho;
    }

    /*
     * Apply entropy fix (section 3.5 in Wada and Liou's paper)
     */
    caseA = ((uL - aL) < 0.0) && ((uR - aR) > 0.0);
    caseB = ((uL + aL) < 0.0) && ((uR + aR) > 0.0);

    d_ua = 0.0;
    if (caseA && !caseB)
	d_ua = C_EFIX * ((uR - aR) - (uL - aL));
    if (caseB && !caseA)
	d_ua = C_EFIX * ((uR + aR) - (uL + aL));

    if (d_ua != 0.0) {
	F.mass -= d_ua * (rR - rL);
	F.momentum.x -= d_ua * (rR * uR - rL * uL);
	F.momentum.y -= d_ua * (rR * vR - rL * vL);
	F.momentum.z -= d_ua * (rR * wR - rL * wL);
	F.total_energy -= d_ua * (rR * HR - rL * HL);
	F.tke -= d_ua * (rR * Rght.tke - rL * Lft.tke);
	F.omega -= d_ua * (rR * Rght.omega - rL * Lft.omega);
    }   /* end of entropy fix (d_ua != 0) */

    for ( int isp = 0; isp < nsp; ++isp ) {
	if (d_ua != 0.0) F.massf[isp] -= d_ua * (rR * Rght.gas->massf[isp] - 
						       rL * Lft.gas->massf[isp]);
    }

    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    for ( int imode = 1; imode < nmodes; ++imode ) {
	if (d_ua != 0.0)
	    F.energies[imode] -= d_ua * (rR * Rght.gas->e[imode] - rL * Lft.gas->e[imode]);
    }

    return SUCCESS;
}   /* end of ausmdv() */

/*--------------------------------------------------------------*/
