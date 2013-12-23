/** \file ausm_plus_up.cxx
 * \ingroup eilmer3
 * \brief Liou's 2006 AUSM+up flux calculator
 *
 * A new version of the AUSM-family schemes, based 
 * on the low Mach number asymptotic analysis.
 * Ironically, this flux calculator causes simulations
 * initialised with 0.0 m/s velocities to crash.
 *
 * RJG -- 26-Apr-2013
 * Added a (+ EPSILON) to help with any divide by zero problems.
 * That being said, I'm not sure this helps with the 
 * crashes at zero velocity because it would seem that the flow
 * of control would pass through a different branch for these
 * cases.
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
#include <numeric>
#include "../../../lib/util/source/useful.h"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"

/*------------------------------------------------------------*/

constexpr double KP = 0.25;
constexpr double KU = 0.75;
constexpr double SIGMA = 1.0;
// Choose a value for M_INF that is good for low Mach numbers.
// To be strictly correct, we should set this at run time
// if an M_INF value is easily defined.
double M_INF = 0.01;
double set_M_inf(double M)
{
    M_INF = M;
    return M_INF;
}
double get_M_inf()
{
    return M_INF;
}

/*------------------------------------------------------------*/

// Some helper functions
double M1plus(double M)
{
    return 0.5*(M + fabs(M));
}

double M1minus(double M)
{
    return 0.5*(M - fabs(M));
}

double M2plus(double M)
{
    return 0.25*(M + 1.0)*(M + 1.0);
}

double M2minus(double M)
{
    return -0.25*(M - 1.0)*(M - 1.0);
}

double M4plus(double M, double beta)
{
    if ( fabs(M) >= 1.0 ) {
	return M1plus(M);
    }
    else {
	double M2p = M2plus(M);
	double M2m = M2minus(M);
	return M2p*(1.0 - 16.0*beta*M2m);
    }
}

double M4minus(double M, double beta)
{
    if ( fabs(M) >= 1.0 ) {
	return M1minus(M);
    }
    else {
	double M2p = M2plus(M);
	double M2m = M2minus(M);
	return M2m*(1.0 + 16.0*beta*M2p);
    }
}

double P5plus(double M, double alpha)
{
    if ( fabs(M) >= 1.0 ) {
	return (1.0/M)*M1plus(M);
    }
    else {
	double M2p = M2plus(M);
	double M2m = M2minus(M);
	return M2p*((2.0 - M) - 16.0*alpha*M*M2m);
    }
}

double P5minus(double M, double alpha)
{
    if ( fabs(M) >= 1.0 ) {
	return (1.0/M)*M1minus(M);
    }
    else {
	double M2p = M2plus(M);
	double M2m = M2minus(M);
	return M2m*((-2.0 - M) + 16.0*alpha*M*M2p);
    }
}

/** \brief Compute the fluxes across an interface. */
int ausm_plus_up(FlowState &Lft, FlowState &Rght, FV_Interface &IFace)
{
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
    double MbarSq, M0Sq; 
    double fa, alpha, beta;
    double M4plus_ML, M4minus_MR;  
    double P5plus_ML, P5minus_MR;
    double Mp, Pu;

    double M_half, p_half, ru_half, ru2_half;

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
    eL = accumulate(Lft.gas->e.begin(), Lft.gas->e.end(), 0.0);
    keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    HL = eL + pL/rL + keL;

    rR = Rght.gas->rho;
    pR = Rght.gas->p;
    uR = Rght.vel.x;
    vR = Rght.vel.y;
    wR = Rght.vel.z;
    aR = Rght.gas->a;
    eR = accumulate(Rght.gas->e.begin(), Rght.gas->e.end(), 0.0);
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
    M0Sq = min(1.0, max(MbarSq, M_INF));
    /*
     * Some additional parameters.
     */
    fa = sqrt(M0Sq) * (2.0 - sqrt(M0Sq));   // eqn 72
    alpha = 0.1875 * (-4.0 + 5 * fa * fa);  // eqn 76
    beta = 0.125;                           // eqn 76

    /*
     * Left state: 
     * M4plus(ML)
     * P5plus(ML)
     */
    M4plus_ML = M4plus(ML, beta);
    P5plus_ML = P5plus(ML, alpha);
    
    /*
     * Right state: 
     * M4minus(MR)
     * P5minus(MR)
     */
    M4minus_MR = M4minus(MR, beta);
    P5minus_MR = P5minus(MR, alpha);

    /*
     * Pressure diffusion modification for 
     * mass flux (eqn 73) and pressure flux (eqn 75).
     */
    double r_half = 0.5*(rL + rR);
    Mp = -KP / fa * max((1.0 - SIGMA * MbarSq), 0.0) * (pR - pL) / (r_half*a_half*a_half);
    Pu = -KU * P5plus_ML * P5minus_MR * (rL + rR) * fa * a_half * (uR - uL);
    /*
     * Mass Flux (eqns 73 & 74).
     */
    M_half = M4plus_ML + M4minus_MR + Mp;
    ru_half = a_half * M_half;
    if ( M_half > 0.0 ) {
       ru_half *= rL;
    } else {
       ru_half *= rR;
    }
    /*
     * Pressure flux (eqn 75).
     */
    p_half = P5plus_ML*pL + P5minus_MR*pR + Pu;

    /*
     * Momentum flux: normal direction
     */
    if ( ru_half >= 0.0 ) {
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
