/** \file ausm_plus_up.cxx
 * \ingroup eilmer3
 * \brief Liou's 2006 AUSM+up flux calculator
 *
 * A new version of the AUSM-family schemes, based 
 * on the low Mach number asymptotic analysis.
 * Ironically, this flux calculator causes simulations
 * initialised with 0.0 m/s velocities to crash.
 *
 * \verbatim
 * M. -S. Liou (2006)
 * A sequel to AUSM, Part II: AUSM+-up for all speeds
 * Journal of Computational Physics, Vol 214, pp 137-170
 * \endverbatim
 *
 * \author W. Y. K. Chan & P. A. Jacobs
 *     Department of Mechanical Engineering
 *     The University of Queensland
 *
 *
 * \version
 * \verbatim
 * 16-Dec-10: Initial coding.
 * \endverbatim
 */

#include <stdio.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"

/*------------------------------------------------------------*/

#define  KP     0.25
#define  KU     0.75
#define  SIGMA  1.0

/*------------------------------------------------------------*/

/** \brief Compute the fluxes across an interface. */
int ausm_plus_up(FlowState &Lft, FlowState &Rght, FV_Interface &IFace)
{
    // if ( get_shock_fitting_flag() ) {
    // 	cerr << "Error, we have not implemented AUSM_PLUS_UP with shock fitting. Please use AUSMDV." << endl;
    // 	exit(NOT_IMPLEMENTED_ERROR);
    // }
    double rL, rR;
    double pL, pR;
    double uL, uR;
    double aL, aR;
    double vL, vR;
    double wL, wR;
    double HL, HR;
    double eL, eR;
    double keL, keR;
    double a_half;
    double ML, MR;
    double MbarSq; 
    double MinfL, MinfR;
    double MoSq, MinfSq;
    double fa, alpha, beta;
    double M1plus, M1minus;  
    double M2plus, M2minus;
    double M4plus, M4minus;  
    double P5plus, P5minus;
    double Mp, Pu;

    double m_half, p_half, ru_half, ru2_half;

    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();

    /*
     * Unpack the flow-state vectors for either side of the interface.
     * Store in work vectors, those quantities that will be needed later.
     */
    rL = Lft.gas->rho;
    pL = Lft.gas->p;
    uL = Lft.vel.x;
    vL = Lft.vel.y;
    wL = Lft.vel.z;
    aL = Lft.gas->a;
    eL = Lft.gas->e[0];
    keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    HL = eL + pL/rL + keL;

    rR = Rght.gas->rho;
    pR = Rght.gas->p;
    uR = Rght.vel.x;
    vR = Rght.vel.y;
    wR = Rght.vel.z;
    aR = Rght.gas->a;
    eR = Rght.gas->e[0];
    keR = 0.5 * (uR * uR + vR * vR + wR * wR);
    HR = eR + pR/rR + keR;

    /*
     * This is the main part of the flux calculator.
     */
    /*
     * Interface sound speed (eqns 28 & 30). 
     * An approximation is used instead of these equations as
     * suggested by Liou in his paper (see line below eqn 69).
     */
    a_half = 0.5 * (aR + aL);
    /*
     * Left and right state Mach numbers (eqn 69).
     */
    ML = uL / a_half;
    MR = uR / a_half;
    /*
     * Mean local Mach number (eqn 70).
     */
    MbarSq = (uL*uL + uR*uR) / (2.0 * a_half *a_half);
    /*
     * Reference Mach number (eqn 71).
     */
    /* Left and right state freestream Mach number */
    MinfL = sqrt(uL*uL + vL*vL + wL*wL) / aL;
    MinfR = sqrt(uR*uR + vR*vR + wR*wR) / aR;
    /* Mean of square of freestream Mach number */
    MinfSq = 0.5 * (MinfL*MinfL + MinfR*MinfR);
    /* Reference Mach number */
    MoSq = MINIMUM(1.0, MAXIMUM(MbarSq, MinfSq));
    /*
     * Some additional parameters.
     */
    fa = sqrt(MoSq) * (2.0 - sqrt(MoSq));   // eqn 72
    alpha = 0.1875 * (-4.0 + 5 * fa * fa);  // eqn 76
    beta = 0.125;                           // eqn 76
    /*
     * Left and right Mach number splitting (eqns 18 & 19).
     * Preliminary calculation of 1st and 2nd degree
     * polynomial functions.
     */    
    M1plus = 0.5 * (ML + FABS(ML)); 
    M1minus = 0.5 * (MR - FABS(MR));
    M2plus = 0.25 * (ML + 1.0) * (ML + 1.0);
    M2minus = -0.25 * (MR - 1.0) * (MR - 1.0);
    /*
     * Left state: 
     * Mach number splitting (eqn 20).
     * and pressure splitting (eqn 24). 
     */
    if (FABS(ML) < 1.0) {
        M4plus = M2plus * (1.0 - 16.0 * beta * M2minus);
	P5plus = M2plus * ((2.0 - ML) - 16.0 * alpha * ML * M2minus);
    } else {
        M4plus = M1plus;
	P5plus = M1plus / ML;
    }
    /*
     * Right state: 
     * Mach number splitting (eqn 20).
     * and pressure splitting (eqn 24).
     */
    if (FABS(MR) < 1.0) {
        M4minus = M2minus * (1.0 + 16.0 * beta * M2plus);
	P5minus = M2minus * ((-2.0 - MR) + 16.0 * alpha * MR * M2plus);
    } else {
        M4minus = M1minus;
	P5minus = M1minus / MR;
    }
    /*
     * Pressure diffusion modification for 
     * mass flux (eqn 73) and pressure flux (eqn 75).
     */
    Mp = -KP / fa * MAXIMUM((1.0 - SIGMA * MbarSq), 0.0) * 
         2.0 * (pR - pL) / (rL + rR) / (a_half * a_half);
    Pu = -KU * P5plus * P5minus * (rL + rR) * fa * a_half * (uR - uL);
    /*
     * Mass Flux (eqns 73 & 74).
     */
    m_half = M4plus + M4minus + Mp;
    ru_half = a_half * m_half;
    if (m_half > 0) {
       ru_half = ru_half * rL;
    } else {
       ru_half = ru_half * rR;
    }
    /*
     * Pressure flux (eqn 75).
     */
    p_half = P5plus * pL + P5minus * pR + Pu;
    /*
     * Momentum flux: normal direction
     */
    if (ru_half > 0) {
       ru2_half = ru_half * uL;
    } else {
       ru2_half = ru_half * uR;
    }
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

    return SUCCESS;
}   /* end of ausm_plus_up() */

/*--------------------------------------------------------------*/
