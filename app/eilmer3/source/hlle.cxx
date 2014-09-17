/** \file hlle.cxx
 * \ingroup eilmer3
 * \brief HLLE fluxes for MHD.
 *
 * From V. Wheatley Matlab implementation

 * \author D. M. Bond
 *     Department of Mechanical Engineering
 *     The University of Queensland
 *
 *
 * \version
 * \verbatim
 * 03-Sep-2012: Initial coding.
 * \endverbatim
 */

#include <stdio.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "cell.hh"
#include "flux_calc.hh"
#include "kernel.hh"

/*------------------------------------------------------------*/

/*------------------------------------------------------------*/

#define SAFESQRT(x) ( (fabs(x)>1.0e-14) ? (sqrt(x)) : (0.0) )


/** \brief Compute the fluxes across an interface. */
int hlle(FlowState &Lft, FlowState &Rght, FV_Interface &IFace)
{
    global_data &G = *get_global_data_ptr();
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();

    /*
     * Unpack the flow-state vectors for either side of the interface.
     * Store in work vectors, those quantities that will be neede later.
     */
    double rL, pL, uL, vL, wL, BxL, ByL, BzL, psiL;

    rL = Lft.gas->rho;
    pL = Lft.gas->p;
    uL = Lft.vel.x;
    vL = Lft.vel.y;
    wL = Lft.vel.z;
    BxL = Lft.B.x;
    ByL = Lft.B.y;
    BzL = Lft.B.z;
    psiL = Lft.psi;

    double rR, pR, uR, vR, wR, BxR, ByR, BzR, psiR;

    rR = Rght.gas->rho;
    pR = Rght.gas->p;
    uR = Rght.vel.x;
    vR = Rght.vel.y;
    wR = Rght.vel.z;
    BxR = Rght.B.x;
    ByR = Rght.B.y;
    BzR = Rght.B.z;
    psiR = Rght.psi;

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

    // Compute Roe Average State (currently simple average)

    double rho, p, u, Bx, By, Bz;

    rho = 0.5*(rL+rR);
    p   = 0.5*(pL+pR);
    u   = 0.5*(uL+uR);
    //v   = 0.5*(vL+vR);
    //w   = 0.5*(wL+wR);
    Bx  = BxL + 0.5*(BxR-BxL) - 0.5*(psiR-psiL)/G.c_h; //0.5*(BxL+BxR);
    By  = 0.5*(ByL+ByR);
    Bz  = 0.5*(BzL+BzR);

    // Compute Eigenvalues of Roe Matrix

    double a2;
    double Bx2, Bt2, BB, ca2;
    double alf, als, cf2, cf;

    //u2=u*u;
    //v2=v*v;
    //w2=w*w;
    //uu=u2+v2+w2;
    a2 = gam*p/rho;
    Bx2 = Bx*Bx;
    Bt2 = By*By + Bz*Bz;
    BB = Bx2 + Bt2;
    ca2 = Bx2/rho;
    alf = a2+BB/rho;
    als = SAFESQRT(alf*alf-4.0*a2*ca2);
    cf2 = 0.5*(alf+als);
    cf = sqrt(cf2);

    double wp, wm;
    wp = u+cf;
    wm = u-cf;

    // Compute the Jump in Conserved Variables between L and R
    double BxL2, BtL2, BBL, ptL;
    double uL2, uuL, aL2, caL2;
    double alfL, alsL, cfL2, cfL;
    double wmL;

    BxL2 = BxL*BxL;
    BtL2 = ByL*ByL + BzL*BzL;
    BBL = BxL2 + BtL2;
    ptL = pL + 0.5*BBL;
    uL2 = uL*uL;
    uuL = uL2 + vL*vL + wL*wL;
    aL2 = gam*pL/rL;
    caL2 = BxL2/rL;
    alfL = aL2+BBL/rL;
    alsL = alfL*alfL-4.0*aL2*caL2;
    alsL = SAFESQRT(alfL*alfL-4.0*aL2*caL2);
    cfL2 = 0.5*(alfL+alsL);
    cfL = sqrt(cfL2);
    //wpL = uL+cfL;
    wmL = uL-cfL;

    double BxR2, BtR2, BBR, ptR;
    double uR2, uuR, aR2, caR2;
    double alfR, alsR, cfR2, cfR;
    double wpR;

    BxR2 = BxR*BxR;
    BtR2 = ByR*ByR + BzR*BzR;
    BBR = BxR2 + BtR2;
    ptR = pR + 0.5*BBR;
    uR2 = uR*uR;
    uuR = uR2 + vR*vR + wR*wR;
    aR2 = gam*pR/rR;
    caR2 = BxR2/rR;
    alfR = aR2+BBR/rR;
    alsR = SAFESQRT(alfR*alfR-4.0*aR2*caR2);
    cfR2 = 0.5*(alfR+alsR);
    cfR = sqrt(cfR2);
    wpR = uR+cfR;
    //wmR = uR-cfR;

    double dU[8];

    dU[0] = rR - rL;
    dU[1] = rR*uR - rL*uL;
    dU[2] = rR*vR - rL*vL;
    dU[3] = rR*wR - rL*wL;
    dU[4] = BxR - BxL;
    dU[5] = ByR - ByL;
    dU[6] = BzR - BzL;
    dU[7] = (pR - pL)/(gam-1.0) + 0.5*(rR*uuR+BBR) - 0.5*(rL*uuL+BBL);

    double bl, br, blm, brp;
    bl = min(wmL, wm);
    br = max(wpR, wp);
    blm = min(bl, 0.0);
    brp = max(br, 0.0);

    double fmassL, fmassR;
    fmassL = rL*uL;
    fmassR = rR*uR;

    double fmomxL, fmomxR;
    fmomxL = rL*uL2 - BxL2 + ptL;
    fmomxR = rR*uR2 - BxR2 + ptR;

    double fmomyL, fmomyR;
    fmomyL = rL*uL*vL - BxL*ByL;
    fmomyR = rR*uR*vR - BxR*ByR;

    double fmomzL, fmomzR;
    fmomzL = rL*uL*wL - BxL*BzL;
    fmomzR = rR*uR*wR - BxR*BzR;

    double fBxL, fBxR;
    fBxL = 0.0;
    fBxR = 0.0;

    double fByL, fByR;
    fByL = uL*ByL - vL*BxL;
    fByR = uR*ByR - vR*BxR;

    double fBzL, fBzR;
    fBzL = uL*BzL - wL*BxL;
    fBzR = uR*BzR - wR*BxR;

    double fenergyL, fenergyR;
    fenergyL = (pL/(gam-1.0)+0.5*(rL*uuL+BBL)+ptL)*uL - (uL*BxL+vL*ByL+wL*BzL)*BxL;
    fenergyR = (pR/(gam-1.0)+0.5*(rR*uuR+BBR)+ptR)*uR - (uR*BxR+vR*ByR+wR*BzR)*BxR;

    double iden, fac1;
    iden = 1.0/(brp - blm);
    fac1 = brp*blm;

    ConservedQuantities &F = *(IFace.F);

    F.mass = (brp*fmassL   - blm*fmassR   + fac1*dU[0])*iden;

    F.momentum.x = (brp*fmomxL   - blm*fmomxR   + fac1*dU[1])*iden;
    F.momentum.y = (brp*fmomyL   - blm*fmomyR   + fac1*dU[2])*iden;
    F.momentum.z = (brp*fmomzL   - blm*fmomzR   + fac1*dU[3])*iden;

    F.B.x = (brp*fBxL     - blm*fBxR     + fac1*dU[4])*iden;
    F.B.y = (brp*fByL     - blm*fByR     + fac1*dU[5])*iden;
    F.B.z = (brp*fBzL     - blm*fBzR     + fac1*dU[6])*iden;

    F.total_energy = (brp*fenergyL - blm*fenergyR + fac1*dU[7])*iden;

    /*
     * Species mass fluxes
     */
    for ( int isp = 0; isp < nsp; ++isp ) {
	if (F.mass >= 0.0) {
	    /* Wind is blowing from the left */
	    F.massf[isp] = F.mass * Lft.gas->massf[isp];
	} else {
	    /* Wind is blowing from the right */
	    F.massf[isp] = F.mass * Rght.gas->massf[isp];
	}
    } /* isp loop */
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
}   /* end of hlle() */

/*--------------------------------------------------------------*/
