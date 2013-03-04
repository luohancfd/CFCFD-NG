/** \file efm.cxx
 * \ingroup eilmer3
 *
 * \brief Mike Macrossan's Equilibrium Flux Method via Paul Petrie.
 *
 * \author Paul Petrie
 * \author PJ
 *
 * \version 02-May-96: Code obtained from Paul Petrie and adapted for mb_cns.
 * \version 03-Nov-96: Update species mass flux code.
 *            Have to fudge the gas constants.
 * \version 26-Jun-97: density weighting of gas constants
 * \version 12-Oct-97: flux_data_1D used to return the flux data
 *            array data
 *
 */

#include <stdio.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"

#if _CRAY
#  pragma _CRI inline exxef
#endif

/*------------------------------------------------------------*/

const double PHI = 1.0;

/** \brief Compute the fluxes across an interface using
 * the Equilibrium Flux	Method of Macrossan & Pullin
 *
 * \param Lft    : IN     : array of Left flow states
 *     (with velocities in local frame of reference)
 * \param Rght   : IN     : array of Right flow state
 * \param IF     : IN/OUT : array of interface flux data structures
 *
 * \verbatim
 * interface data -- contains...
 *     flux of mass across the interface (kg/s/m**2)
 *     flux of normal momentum
 *     flux of tangential momentum
 *     flux of energy
 *     array of species fluxes 
 *     vibrational energies
 *     free-electron energy
 * \endverbatim 
 */
int efmflx(FlowState &Lft, FlowState &Rght, FV_Interface &IFace)
{
    // if ( get_shock_fitting_flag() ) {
    // 	cerr << "Error, we have not implemented EFM with shock fitting. Please use AUSMDV." << endl;
    // 	exit(NOT_IMPLEMENTED_ERROR);
    // }
    /*
     * Local variable names reflect the names used in the original
     * FORTRAN code by MNM.
     */
    double vnL, vpL, vnR, vpR, vqL, vqR;
    double rtL, cmpL, rtR, cmpR;
    double hvsqL, hvsqR;
    double wL, wR, dL, dR;
    double rhoL, rhoR, presL, presR, tR, tL;
    double hL, hR;
    double snL, snR, exL, exR, efL, efR;
    double fmsL, fmsR;
    double cv, cp, con, gam, Rgas;
    double cvL, cvR, RgasL, RgasR;
    double rLsqrt, rRsqrt, alpha;
    Gas_model *gmodel = get_gas_model_ptr();
    int statusf;

    /* Collect the global gas constants for later use. */
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();

    /*
     * Calculate Constants
     */
    /* dtwspi = 1.0 / (2.0 * sqrt ( 3.14159265359 )); */
    const double dtwspi = 0.282094792;

    /*
     * Unpack Left flow state
     */
    rhoL = Lft.gas->rho;
    presL = Lft.gas->p;
    hL = Lft.gas->e[0] + presL/rhoL;
    tL = Lft.gas->T[0];
    vnL = Lft.vel.x;
    vpL = Lft.vel.y;
    vqL = Lft.vel.z;

    /*
     * Unpack Right flow state
     */
    rhoR = Rght.gas->rho;
    presR = Rght.gas->p;
    hR = Rght.gas->e[0] + presR/rhoR;
    tR = Rght.gas->T[0];
    vnR = Rght.vel.x;
    vpR = Rght.vel.y;
    vqR = Rght.vel.z;

    /* Derive the gas "constants" from the local conditions. */
    cvL = gmodel->Cv(*(Lft.gas), statusf);
    RgasL = presL / (rhoL * tL);
    cvR = gmodel->Cv(*(Rght.gas), statusf);
    RgasR = presR / (rhoR * tR);

    rLsqrt = sqrt(rhoL);
    rRsqrt = sqrt(rhoR);
    alpha = rLsqrt / (rLsqrt + rRsqrt);

    cv = alpha * cvL + (1.0 - alpha) * cvR;
    Rgas = alpha * RgasL + (1.0 - alpha) * RgasR;

    cp = cv + Rgas;
    gam = cp / cv;

    /*
     * Start EFM calculation proper.
     */
    con = 0.5 * (gam + 1.0) / (gam - 1.0);

    rtL = Rgas * tL;
    cmpL = sqrt(2.0 * rtL);
    hvsqL = 0.5 * (vnL * vnL + vpL * vpL + vqL * vqL);

    snL = vnL / (PHI * cmpL);
    exxef(snL, exL, efL);

    wL = 0.5 * (1.0 + efL);
    dL = exL * dtwspi;

    rtR = presR / rhoR;
    cmpR = sqrt(2.0 * rtR);
    hvsqR = 0.5 * (vnR * vnR + vpR * vpR + vqR * vqR);

    snR = vnR / (PHI * cmpR);
    exxef(snR, exR, efR);

    wR = 0.5 * (1.0 - efR);
    dR = -exR * dtwspi;

    /*
     * combine the fluxes
     */
    fmsL = (wL * rhoL * vnL) + (dL * cmpL * rhoL);
    fmsR = (wR * rhoR * vnR) + (dR * cmpR * rhoR);

    ConservedQuantities &F = *(IFace.F);
    F.mass = fmsL + fmsR;

    F.momentum.x = fmsL * vnL + fmsR * vnR + wL * presL + wR * presR;
    F.momentum.y = fmsL * vpL + fmsR * vpR;
    F.momentum.z = fmsL * vqL + fmsR * vqR;

    F.total_energy = (wL * rhoL * vnL) * (hvsqL + hL) +
	(wR * rhoR * vnR) * (hvsqR + hR) +
	(dL * cmpL * rhoL) * (hvsqL + con * rtL) +
	(dR * cmpR * rhoR) * (hvsqR + con * rtR);

    if ((F.mass) > 0.0) {
	F.tke = (F.mass) * Lft.tke;
	F.omega = (F.mass) * Lft.omega;
    } else {
	F.tke = (F.mass) * Rght.tke;
	F.omega = (F.mass) * Rght.omega;
    }

    /* 
     * Species mass flux.
     * Presently, this is implemented by assuming that
     * the wind is blowing one way or the other and then
     * picking the appropriate side for the species fractions.
     * Such an approach may not be fully compatible with the
     * EFM approach where there can be fluxes from both sides.
     */
    for ( int isp = 0; isp < nsp; ++isp ) {
	if ((F.mass) > 0.0)
	    F.massf[isp] = (F.mass) * Lft.gas->massf[isp];
	else
	    F.massf[isp] = (F.mass) * Rght.gas->massf[isp];
    }   /* isp loop */

    // Individual energies.
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if ((F.mass) > 0.0) {
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = (F.mass) * Lft.gas->e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += (F.mass) * Lft.gas->p_e / Lft.gas->rho;
    } else {
	for ( int imode = 1; imode < nmodes; ++imode ) {
	    F.energies[imode] = (F.mass) * Rght.gas->e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += (F.mass) * Rght.gas->p_e / Rght.gas->rho;
    }
    return 0;
}   /* end efmflx() */

/*--------------------------------------------------------------*/

/** \brief Compute exp(-x**2) and erf(x) with a polynomial approximation.
 *
 * \param sn   : IN  : x
 * \param &exx : OUT : exp(x**2)
 * \param &ef  : OUT : erf(x)  error function
 */
int exxef( double sn, double &exx, double &ef )
{
    double snsq, ef1, y;

    const double P = 0.327591100;
    const double A1 = 0.254829592;
    const double A2 = -0.284496736;
    const double A3 = 1.421413741;
    const double A4 = -1.453152027;
    const double A5 = 1.061405429;
    const double LIMIT = 5.0;
    const double EXLIM = 0.138879e-10;
    const double EFLIM = 1.0;

#   define DSIGN(val,sgn) ( (sgn >= 0.0)? fabs(val): -fabs(val) )

    if (fabs(sn) > LIMIT) {
        exx = EXLIM;
        ef1 = EFLIM;
    } else {
        snsq = sn * sn;
        exx = exp(-snsq);
        y = 1.0 / (1.0 + P * fabs(sn));
        ef1 = 1.0 - y * (A1 + y * (A2 + y * (A3 + y * (A4 + A5 * y)))) * exx;
    }
    ef = DSIGN(ef1, sn);
    return 0;
}   /* end exxef() */
