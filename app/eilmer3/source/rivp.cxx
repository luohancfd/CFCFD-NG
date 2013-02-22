/** \file rivp.cxx
 * \ingroup basic_gas_dyn
 * \brief One-dimensional Riemann Initial Value Problem.
 *
 * \todo Really should get rid of the vector loop some time.
 * \todo Need to rework interpolation and remove the assumption 
 *       that e = Cv.T
 */

/*------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "flux_calc.hh"
#include "kernel.hh"
#include "../../../lib/nm/source/qd_power.h"

/*------------------------------------------------------------*/

/** \brief Compute flux using a solution to the 1D Riemann problem.
 *
 *
 * Given the initial LEFT and RIGHT states (QL and QR)
 * either side of an interface, compute the solution to
 * the one-dimensional Riemann problem when the interface
 * is removed.  This solution is returned as values of
 * the flow quantities at the interface (QIF) and the
 * velocities of the left and right waves.
 *
 * \param &QL   : IN : reference to the LEFT state
 *                     see "flux_calc.hh" for a definition of the components
 * \param &QR   : IN : reference to the RIGHT state
 * \param &QIF  : OUT : Flow state at the interface
 * \param &WSL  : OUT : Left wave velocity
 * \param &WSR  : OUT : Right wave velocity
 *
 * \version
 * \verbatim
 * 1.0,  20 Dec 90
 * 1.01, 20-Feb-91: include "compilers.h"
 * 1.10, 17-Mar-91: installed quick and dirty power routine
 * 1.2,  25-Apr-91: vacuum case installed, WSL and WSR added
 * 1.3,  26-Apr-91: vector version
 * 1.4,  27-May-91: strong shock solution added.
 * 1.4.1 28-May-91: Newton steps added as a second stage
 * 1.4.2 05-Jun-91: changed the way intermediate variables are
 *                  backed out
 * 2.0   18-Nov-92: Split the interpolation section off into
 *                  the new function rivp_stage_3().
 *                  Made the stage numbering consistent with
 *                  the AIAA Journal paper.
 * 3.0   04-Mar-93: Split the main loop into smaller parts
 *                  for the Fujitsu compiler
 *                  Promote scalar temporaries to vectors.
 * 3.01  10-Mar-93: Have removed the 4 Newton-step code but
 *                  have left the two step code.
 *                  Have introduced a new temporary variable
 *                  into the calculation of the wave speeds.
 * 4.0   17-Mar-93: Added Toro's linearized Riemann solver.
 * 5.0   13-Aug-93: Added general-strength shock expressions 
 *                  to stage 2.
 * 5.01  16-Aug-93: Finally got the signs right in the Newton steps.
 * 5.1   29-Oct-93: Fixed a couple of problems for cases where
 *                  the Newton steps run into trouble.  Also
 *                  put the number of Newton steps back to 4.
 * 6.0   29-Jan-96: Thermo properties now come via service functions.
 * 6.01  29-Apr-96: Promoted DEBUG output to level 4 so that not so much
 *                  is written.
 * 6.02  03-Nov-96: Removed Toro's linearised solver and
 *                  the strong-shock code which hasn't been
 *                  used for ages.
 * --------------------------------------------------------------------------
 * 08-Sep-2002 : Started to remove the e = C_v * T dependency from the
 *               code here but it is deeply embedded in the later stages
 *               of the calculations.  Maybe it just isn't worth bothering
 *               to keep it up-to-date.  Consider these functions deprecated!
 * --------------------------------------------------------------------------
 *\endverbatim
 *
 * \author P.A. Jacobs
 * ICASE
 * Mail Stop 123C
 * NASA Langley Rearch Centre
 * Hampton VA 23665.
 * 
 * \author PJ
 * Department of Mechanical Engineering
 * The University of Queensland
 * Brisbane 4072
 *
 * \verbatim
 * References ...
 * ----------
 *
 * Most of this solver is now documented in
 * P. A. Jacobs "An approximate Riemann solver for hypervelocity
 * flows."  AIAA Journal vol. 30(10) 1992.
 *
 * S.R. Chakravarthy "Development of upwind schemes for the
 * Euler equations"
 * NASA Contractor Report 4043, 1987.
 *
 * P. Colella "Glimm's method for gas dynamics"
 * SIAM J. Sci. Stat. Comput. vol.3, 76-110, 1982.
 *
 * T. A. Edwards (1988)
 * "The effect of exhaust plume/afterbody interaction on
 * installed scramjet performance."
 * NASA Technical Memorandum 101033
 *
 * H.M. Glaz & A.B. Wardlaw "A high-order Godunov scheme for
 * steady supersonic gas dynamics" J. Comput. Phys. vol.58,
 * pp157-187, 1985.
 *
 * J.J. Gottlieb & C.P.T. Groth "Assessment of Riemann solvers
 * for unsteady one-dimensional inviscid flows of perfect gases"
 * J. Comput. Phys. vol.78(2), 437-458, 1988.
 *
 * B. van Leer "On the relation between the upwind-differencing
 * schemes of Godunov, Engquist-Osher & Roe"
 * SIAM J. Sci. Stat. Comput. vol.5(1), 1-20, 1984.
 *
 * B. van Leer, W-T. Lee & K.G. Powell "Sonic point capturing"
 * AIAA-89-1945-CP (AIAA 9th CFD conference), 1989.
 *
 * H.W. Liepmann & A.Roshko "Elements of Gas Dynamics"
 * Wiley, 1957.
 *
 * S. Osher & F. Solomon "Upwind difference schemes for
 * hyperbolic systems of conservation laws"
 * Mathematics of Computation, vol.38, 339-, 1982.
 *
 * P.L. Roe "Some contributions to the modelling of
 * discontinuous flows" Lectures in Applied Mathematics,
 * Vol. 22 (part 2), 163-193, 1985
 * 
 * E. F. Toro "A linearized Riemann solver for the time-
 * dependent Euler equations of gas-dynamics"
 * Proc. R. Soc. Lond. A Vol. 434, 683-693 (1991)
 *
 * M. Vinokur & Y. Liu (1988)
 * "Equilibrium gas flow computations II: An analysis of
 * numerical formulations of conservation laws."
 * AIAA Paper 88-0127
 * \endverbatim
 */
int rivp(FlowState &QL, FlowState &QR, FlowState &QIF, double &WSL, double &WSR)
{
    // if ( get_shock_fitting_flag() ) {
    // 	cerr << "Error, we have not implemented RIVP with shock fitting. Please use AUSMDV." << endl;
    // 	exit(NOT_IMPLEMENTED_ERROR);
    // }
    Gas_model *gmodel = get_gas_model_ptr();
    int statusf;
    FlowState QLstar(gmodel);
    FlowState QRstar(gmodel);
    /* Intermediate flow states */
    double geff; /* effective GAMMA          */
    double gL, gR, sqrL, sqrR, alpha;
    double gm1, gp1, z;
    double uLbar, uRbar, pstar, ustar;
    double rhoLstar, rhoRstar, eLstar, eRstar;
    double aLstar, aRstar, TLstar, TRstar;
    double wspeedL, wspeedR;
    int vacuum, toro_solve;

    double base, expon, pwr;    /* for qd_power MACRO */
    double temporary;
    int iter;
    double term1, term2, F, dFdpstar, delp;

    /* double p_min, p_max, rho_bar, a_bar;   FIX-ME CHECK-ME */
    double PMIN, TMIN, AMIN, RHOMIN;

    /*
     * Set the following macro to 1 if the quick-and-dirty
     * evaluation of the power function is to be used.
     */
#   define  QUICK_AND_DIRTY  1

    RHOMIN = 0.0;
    PMIN = 0.0;
    TMIN = 1.0;
    AMIN = 0.0;

    sqrL = sqrt(QL.gas->rho);
    sqrR = sqrt(QR.gas->rho);
    vacuum = 0;
    toro_solve = 0;

    /*************************************************************
     * Compute a Roe-averaged effective-gamma for each interface *
     * as used by Edwards (1988) and Vinokur & Liu (1988).       *
     *************************************************************/
    gL = gmodel->gamma(*(QL.gas), statusf);
    gR = gmodel->gamma(*(QR.gas), statusf);
    alpha = sqrL / (sqrL + sqrR);
    geff = alpha * gL + (1.0 - alpha) * gR;
    gm1 = geff - 1.0;
    gp1 = geff + 1.0;

    /*****************************************************************
     * Solve one-dimensional problem in a frame of reference         *
     * for which QL(i,uq) is the normal velocity while QL.v and   *
     * QL.w are the transverse velocities.                        *
     *                                                               *
     * First, solve for the intermediate states.                     *
     * STAGE 1:                                                      *
     * We will use an Osher-type approximate Riemann solver.         *
     * One of the original descriptions can be found in              *
     * Osher & Solomon (1982) while the more recent report           *
     * by Chakravarthy (1987) contains a brief description.          *
     * Both waves are considered to be isentropic whether            *
     * they are compressions or expansions.                          *
     * We use the relations:                                         *
     * [1] p / rho**gamma = constant                                 *
     * [2] p = rho * e * (gamma - 1.0)                               *
     * [3] J+ = constant = u + 2 a / (gamma - 1) for the (u-a) wave  *
     * (Note that we integrate along a (u+a) characteristic          *
     * through this (u-a) wave.)                                     *
     * [4] J- = constant = u - 2 a / (gamma - 1) for the (u+a) wave  *
     *                                                               *
     * STAGE 2:                                                      *
     * If the intermediate pressure is much larger than that for     *
     * both initial states, the solution is replaced by the strong   *
     * shock solution (consisting of 2 strong shocks).               *
     * If we have a combination of one strong shock and one          *
     * isentropic wave then we take a couple of Newton steps to get  *
     * a better estimate for pstar, ustar.                           *
     * Explicitly coding the steps allows the main loop to           *
     * vectorize (at least on the Cray Standard C compiler).         *
     *                                                               *
     * STAGE 2a:                                                      *
     * Back out the other intermediate variables using either the    *
     * shock relations or the isentropic wave relations.             *
     *                                                               *
     *****************************************************************/

    /*
     * -----------------------------------------------------
     * STAGE 1: Explicit solution using two isentropic waves.
     * -----------------------------------------------------
     */
#   if DEBUG_FLUX >= 1
    printf("STAGE 1\n");
#   endif

    if (!toro_solve) {
	/*
	 * Intermediate variable. 
	 */
	base = QL.gas->p / QR.gas->p;
	expon = gm1 / (2.0 * geff);
#       if (QUICK_AND_DIRTY == 1)
	/* We expect exponent to be <= 0.2 */
	if (base < 0.12 || base > 8.4) {
	    /* Use the full power function */
	    pwr = pow(base, expon);
	} else {
	    /* Use the quick and dirty approach */
	    QD_POWER(base, expon, pwr);
	}
#       else
	pwr = pow(base, expon);
#       endif
	z = (QR.gas->a / QL.gas->a) * pwr;
    }   /* end if (!toro_solve... */

    /* 
     * Riemann invariants. 
     */
    if (!toro_solve) {
	uLbar = QL.vel.x + 2.0 / gm1 * QL.gas->a;
	uRbar = QR.vel.x - 2.0 / gm1 * QR.gas->a;
    }

    if (!toro_solve && (uLbar - uRbar) <= 0.0) {
	/*
	 * We have a situation in which a (near) vacuum is formed
	 * between the left and right states.
	 */
	vacuum = 1;

	ustar = 0.0;
	pstar = PMIN;

	eLstar = 1.0;  /* arbitrary, positive value */ 
	eRstar = 1.0;

	rhoLstar = RHOMIN;
	rhoRstar = RHOMIN;

	TLstar = TMIN;
	TRstar = TMIN;

	aLstar = AMIN;
	aRstar = AMIN;

    }   /* End of Vacuum solution */

    if (!toro_solve && !vacuum) {
	/*
	 * Positive-pressure solution.
	 */

	/* First, compute ustar. */
	ustar = (uLbar * z + uRbar) / (1.0 + z);

	/* Then, pstar. */
	base = (0.5 * gm1 * (uLbar - uRbar) /
		(QL.gas->a * (1.0 + z)));
	expon = 2.0 * geff / gm1;
#       if (QUICK_AND_DIRTY == 1)
	/* Expect 5.0 <= expon <= 12 */
	if (base < 0.75 || base > 1.5) {
	    /* Use the full power function */
	    pwr = pow(base, expon);
	    pstar = QL.gas->p * pwr;
	} else {
	    /* Reduce the exponent and use the quick and dirty 
	     * approach 
	     */
	    expon -= 8.0;
	    QD_POWER(base, expon, pwr);
	    temporary = base * base;
	    temporary *= temporary;
	    temporary *= temporary;
	    pwr *= temporary;
	    pstar = QL.gas->p * pwr;
	}
#       else
	/* Use the full power function */
	pwr = pow(base, expon);
	pstar = QL.gas->p * pwr;
#       endif

	/*
	 * ------------------------------------------
	 * End of isentropic-wave (explicit) solution.
	 * ------------------------------------------
	 */
#       if DEBUG_FLUX >= 1
	printf("Isentropic wave solution: pstar = %e, ustar = %e\n", pstar, ustar);
#       endif
    }   /* if !vacuum ... */

    /*
     * --------------------------------------------------------
     * STAGE 2: Apply Newton Iterations for the shock relations.
     * --------------------------------------------------------
     */
#   define  SHOCK_RATIO  1.5

    if (!toro_solve && !vacuum) {
	if (pstar > SHOCK_RATIO * QL.gas->p &&
	    pstar > SHOCK_RATIO * QR.gas->p) {
	    /*
	     * Both waves are shocks.
	     */
	    for ( iter = 1; iter <= 4; ++iter ) {
		term1 = sqrt(gp1 / (2.0 * geff) * pstar / QL.gas->p
			     + gm1 / (2.0 * geff));
		term2 = sqrt(gp1 / (2.0 * geff) * pstar / QR.gas->p
			     + gm1 / (2.0 * geff));
		F = QL.vel.x
		    - QL.gas->a / geff * (pstar / QL.gas->p - 1.0) / term1
		    - QR.vel.x
		    - QR.gas->a / geff * (pstar / QR.gas->p - 1.0) / term2;
		dFdpstar =
		    -gp1 * QL.gas->a / (4.0 * geff * geff * QL.gas->p)
		    * (pstar / QL.gas->p + (3.0 * geff - 1.0) / gp1)
		    / (term1 * term1 * term1)
		    - gp1 * QR.gas->a / (4.0 * geff * geff * QR.gas->p)
		    * (pstar / QR.gas->p + (3.0 * geff - 1.0) / gp1)
		    / (term2 * term2 * term2);
		delp = F / dFdpstar;
		if (delp > pstar)
		    pstar *= 0.1;
		else
		    pstar -= delp;
	    } /* end for iter */

	    /* -------------- velocity -------------- */
	    ustar = QL.vel.x
		- QL.gas->a / geff * (pstar / QL.gas->p - 1.0) / term1;

#           if DEBUG_FLUX >= 1
	    printf("Two shocks: pstar = %e, ustar = %e\n", pstar, ustar);
#           endif
	} else if (pstar > SHOCK_RATIO * QR.gas->p) {
	    /*
	     * Treat the right-moving wave as a shock, the
	     * left-moving wave as an isentropic wave, and take
	     * two Newton steps to improve the guess for pstar, ustar.
	     */
	    for ( iter = 1; iter <= 4; ++iter ) {
		term1 = pow(pstar / QL.gas->p, gm1 / (2.0 * geff));
		term2 = sqrt(gp1 / (2.0 * geff) * pstar / QR.gas->p
			     + gm1 / (2.0 * geff));
		F = QL.vel.x - 2.0 * QL.gas->a / gm1 * (term1 - 1.0)
		    - QR.vel.x
		    - QR.gas->a / geff * (pstar / QR.gas->p - 1.0) / term2;
		dFdpstar = -QL.gas->a / (geff * pstar) * term1
		    - gp1 * QR.gas->a / (4.0 * geff * geff * QR.gas->p)
		    * (pstar / QR.gas->p + (3.0 * geff - 1.0) / gp1)
		    / (term2 * term2 * term2);
		delp = F / dFdpstar;
		if (delp > pstar)
		    pstar *= 0.1;
		else
		    pstar -= delp;
	    } /* end for iter */

	    /* ---------------- velocity ------------------- */
	    ustar = QR.vel.x +
		QR.gas->a / geff * (pstar / QR.gas->p - 1.0) / term2;

#           if DEBUG_FLUX >= 1
	    printf("Right shock: pstar = %e, ustar = %e\n", pstar, ustar);
#           endif
	} else if (pstar > SHOCK_RATIO * QL.gas->p) {
	    /*
	     * Treat the left-moving wave as a shock, the
	     * right-moving wave as an isentropic wave, and take
	     * two Newton steps to improve the guess for pstar, ustar.
	     */
	    for ( iter = 1; iter <= 4; ++iter ) {
		term1 = sqrt(gp1 / (2.0 * geff) * pstar / QL.gas->p
			     + gm1 / (2.0 * geff));
		term2 = pow(pstar / QR.gas->p, gm1 / (2.0 * geff));
		F = QL.vel.x
		    - QL.gas->a / geff * (pstar / QL.gas->p - 1.0) / term1
		    - QR.vel.x - 2.0 * QR.gas->a / gm1 * (term2 - 1.0);
		dFdpstar =
		    -gp1 * QL.gas->a / (4.0 * geff * geff * QL.gas->p)
		    * (pstar / QL.gas->p + (3.0 * geff - 1.0) / gp1)
		    / (term1 * term1 * term1)
		    - QR.gas->a / (geff * pstar) * term2;
		delp = F / dFdpstar;
		if (delp > pstar)
		    pstar *= 0.1;
		else
		    pstar -= delp;
	    } /* end for iter */

	    /* ---------------- velocity ------------------- */
	    ustar = QL.vel.x -
		QL.gas->a / geff * (pstar / QL.gas->p - 1.0) / term1;

#           if DEBUG_FLUX >= 1
	    printf("Left shock: pstar = %e, ustar = %e\n", pstar, ustar);
#           endif
	}
    }   /* if !vacuum ... */
    /* End of Newton Iterations for normal shocks */

    /*
     * ---------------------------------------------
     * STAGE 2a: Back out the other quantities using
     *           thermodynamic relations.
     * ---------------------------------------------
     */

    if (!toro_solve && !vacuum) {
	if (pstar > QL.gas->p) {
	    /*
	     * Back out values using the shock relations.
	     * Density -- from the Rankine-Hugoniot relations 
	     */
	    rhoLstar = QL.gas->rho *
		(gp1 * pstar + gm1 * QL.gas->p) /
		(gp1 * QL.gas->p + gm1 * pstar);
	    /* 
	     * Specific energy -- from the Equation of state 
	     */
	    eLstar = pstar / (gm1 * rhoLstar);
	    /* 
	     * Local speed of sound -- Perfect gas version.
	     */
	    aLstar = sqrt(geff * gm1 * eLstar);
	} else {
	    /* 
	     * Use the isentropic-wave relations. 
	     * Local speed of sound -- Riemann invariants. 
	     */
	    aLstar = (uLbar - ustar) * 0.5 * gm1;
	    /* 
	     * Specific energy -- sound speed 
	     */
	    eLstar = aLstar * aLstar / (geff * gm1);
	    /* 
	     * Density -- equation of state
	     */
	    rhoLstar = pstar / (gm1 * eLstar);
	}

	if (pstar > QR.gas->p) {
	    /*
	     * Back out values using the shock relations.
	     * Density -- from the Rankine-Hugoniot relations 
	     */
	    rhoRstar = QR.gas->rho *
		(gp1 * pstar + gm1 * QR.gas->p) /
		(gp1 * QR.gas->p + gm1 * pstar);
	    /* 
	     * Specific energy -- from the Equation of state 
	     */
	    eRstar = pstar / (gm1 * rhoRstar);
	    /*
	     * Local speed of sound -- Perfect gas version. 
	     */
	    aRstar = sqrt(geff * gm1 * eRstar);
	} else {
	    /* 
	     * Use the isentropic-wave relations. 
	     * Local speed of sound -- Riemann invariants. 
	     */
	    aRstar = (ustar - uRbar) * 0.5 * gm1;
	    /* 
	     * Specific energy -- sound speed
	     */
	    eRstar = aRstar * aRstar / (geff * gm1);
	    /* 
	     * Density -- equation of state
	     */
	    rhoRstar = pstar / (gm1 * eRstar);
	}

	/* 
	 * Temperatures -- equation of state also.
	 */
	TLstar = QL.gas->T[0] * (pstar * QL.gas->rho) /
	    (QL.gas->p * rhoLstar);
	TRstar = QR.gas->T[0] * (pstar * QR.gas->rho) /
	    (QR.gas->p * rhoRstar);

    }
    /* if !vacuum ... */

    /*
     * --------------------------------------------------------------
     * At this point, we should have a good estimate for pstar, ustar,
     * and the thermodynamic properties.
     * --------------------------------------------------------------
     */

    /* 
     * ***********
     * Wave speeds. 
     * ***********
     */
    if ((pstar > QL.gas->p) && !vacuum && !toro_solve) {
	/* Left wave is a shock. */
	temporary = 0.5 * gp1 * QL.gas->p / QL.gas->rho *
	    (pstar / QL.gas->p + gm1 / gp1);
	wspeedL = QL.vel.x - sqrt(temporary);
    } else {
	/* Left wave is an expansion fan. */
	wspeedL = QL.vel.x - QL.gas->a;
    }

    if ((pstar > QR.gas->p) && !vacuum && !toro_solve) {
	/* Right wave is a shock. */
	temporary = 0.5 * gp1 * QR.gas->p / QR.gas->rho *
	    (pstar / QR.gas->p + gm1 / gp1);
	wspeedR = QR.vel.x + sqrt(temporary);
    } else {
	/* Right wave is an expansion fan. */
	wspeedR = QR.vel.x + QR.gas->a;
    }

    /*
     * Record solution in vector.
     * This bit of double handling is required only for the 
     * small TopSpeed compiler since it could not cope with
     * the inperpolation section being in the same function.
     */
    QLstar.gas->rho = rhoLstar;
    QRstar.gas->rho = rhoRstar;

    QLstar.vel.x = ustar;
    QRstar.vel.x = ustar;

    QLstar.gas->e[0] = eLstar;
    QRstar.gas->e[0] = eRstar;

    QLstar.gas->p = pstar;
    QRstar.gas->p = pstar;

    QLstar.gas->a = aLstar;
    QRstar.gas->a = aRstar;

    QLstar.gas->T[0] = TLstar;
    QRstar.gas->T[0] = TRstar;

    WSL = wspeedL;
    WSR = wspeedR;


    /***************************************************************
     * Now, interpolate to obtain the interface values.            *
     ***************************************************************/

    rivp_stage_3(QL, QR, QLstar, QRstar, WSL, WSR, geff, QIF);


    return 0;
}   /* end of rivp() */

/*------------------------------------------------------------*/

/** \brief Interpolate interface state from intermediate states.
 *
 *
 * Given the initial LEFT and RIGHT states (QL and QR)
 * either side of an interface and the intermediate
 * states, interpolate (or select) the interface state.
 *
 * \param &QL    : IN  : reference to the LEFT state
 *     see "flux_calc.hh" for a definition of the components
 * \param &QR    : IN  : RIGHT states
 * \param &QLstar, &QRstar  : IN : the intermediate state
 * \param geff  : IN  : effective gamma
 * \param WSL    : IN  : Left wave velocity
 * \param WSR    : IN  : Right wave velocity
 * \param &QIF   : OUT : Flow state at the interface
 *
 * \version 1.0,  18-Nov-92: Split off from rivp()
 * \version 1.01, 03-Nov-96: Update multiple-species code.
 *
 */
int rivp_stage_3(FlowState &QL, FlowState &QR,
		 FlowState &QLstar, FlowState &QRstar,
		 double WSL, double WSR,
		 double geff, FlowState &QIF )
{
    double VA, VB, frac, Tratio;
    int option;
    Gas_model *gmodel = get_gas_model_ptr();
    int nsp = gmodel->get_number_of_species();
    int nmodes = gmodel->get_number_of_modes();

    /***************************************************************
     * Now, interpolate to obtain the interface values.            *
     *                                                             *
     * This is mainly a problem of deciding which way the waves    *
     * went once the imaginary diaphragm was removed.  However,    *
     * if an expansion wave straddles the interface (i.e. the      *
     * so called sonic expansion), we linearly                     *
     * interpolate the flow velocity at the interface and then     *
     * compute the rest of the quantities within the fan from      *
     * the characteristic equations and isentropic relations.      *
     * See for example sec. 3.11 in Liepmann & Roshko (1957).      *
     * The performance of the method seems to be nearly as good    *
     * when linear interpolation is used for all of the flow       *
     * quantities.                                                 *
     *                                                             *
     * The importance of spreading expansion fans is suggested in  *
     * van Leer (1984), Colella (1982) and Glaz & Wardlaw (1985).  *
     * A graphic demonstration is provided by Roe (1985) for flow  *
     * over a forward facing step.  If the expansion waves are     *
     * not spread, an expansion shock forms across the flow.       *
     * Standard interpolation within the expansion fans prevents   *
     * the formation of strong expansion shocks but a small        *
     * amount of extra dissipation is needed to eliminate the      *
     * small "glitch" that may occur in a sonic rarefaction.       *
     *                                                             *
     ***************************************************************/

    /* The default interpolation within the expansion fan.
     * Either linear for all variables or characteristic eqns. */
#   define  LINEAR  1

    /* Decide which way the waves are going. */
    option = 0;

    if (QLstar.vel.x > 0.0) {
	/* The left wave and the contact discontinuity determine
	 * the cell interface quantities. */

	if (QLstar.gas->p - QL.gas->p >= 0.0) {
	    /* The left wave is a shock. */
	    if (WSL >= 0.0) {
		/* All waves have gone into the right cell.
		 * The values are taken from the left cell only. */
		option = 1;
	    } else {
		/* The values are taken from behind the left shock. */
		option = 2;
	    }
	} else {
	    /* The left wave is an expansion. */
	    VA = QL.vel.x - QL.gas->a;
	    VB = QLstar.vel.x - QLstar.gas->a;
	    if (VA >= 0.0) {
		/* All waves have gone into the right cell.
		 * The values are taken from the left cell only. */
		option = 3;
	    } else if (VB > 0.0) {
		/* The left rarefaction straddles the cell interface.
		 * Interpolate velocity inside the left rarefaction. */
		option = 4;
	    } else {
		/* The values come from behind the left rarefaction. */
		option = 5;
	    }
	}
    } else {
	/* The right wave and the contact discontinuity determine
	 * the cell interface quantities. */

	if (QRstar.gas->p - QR.gas->p >= 0.0) {
	    /* The right wave is a shock. */
	    if (WSR <= 0.0) {
		/* All waves have gone into the left cell.
		 * The values are taken from the right cell only. */
		option = 11;
	    } else {
		/* The values are taken from behind the right shock. */
		option = 12;
	    }
	} else {
	    /* The right wave is an expansion. */
	    VA = QR.vel.x + QR.gas->a;
	    VB = QRstar.vel.x + QRstar.gas->a;
	    if (VA <= 0.0) {
		/* All waves have gone into the left cell.
		 * The values are taken from the right cell only. */
		option = 13;
	    } else if (VB < 0.0) {
		/* The right rarefaction straddles the cell interface.
		 * Interpolate velocity inside the right rarefaction. */
		option = 14;
	    } else {
		/* The values come from behind the right rarefaction. */
		option = 15;
	    }
	}
    }

    /************************************************
     *  Now, copy or interpolate the relevant data. *
     ************************************************/

    if (option == 1 || option == 3) {
	/* All waves have gone into the right cell.
	 * The values are taken from the left cell only. */
	QIF.gas->rho = QL.gas->rho;
	QIF.vel.x = QL.vel.x;
	QIF.gas->p = QL.gas->p;
	QIF.gas->e[0] = QL.gas->e[0];
	QIF.gas->a = QL.gas->a;
	QIF.gas->T[0] = QL.gas->T[0];
    }

    if (option == 2 || option == 5) {
	/* The values are taken from behind the left wave. */
	QIF.gas->rho = QLstar.gas->rho;
	QIF.vel.x = QLstar.vel.x;
	QIF.gas->p = QLstar.gas->p;
	QIF.gas->e[0] = QLstar.gas->e[0];
	QIF.gas->a = QLstar.gas->a;
	QIF.gas->T[0] = QLstar.gas->T[0];
    }

    if (option == 4) {
	/* The left wave is an expansion. */
	VA = QL.vel.x - QL.gas->a;
	VB = QLstar.vel.x - QLstar.gas->a;
	/* The left rarefaction straddles the cell interface.
	 * Interpolate velocity inside the left rarefaction. */
	frac = (-VA) / (VB - VA);
	QIF.vel.x = QL.vel.x - frac * (QL.vel.x - QLstar.vel.x);
#       if (LINEAR != 1)
	/* Get local speed of sound by following a (u+a)
	 * characteristic line through the fan. */
	gm1 = geff - 1.0;
	QIF.gas->a = QL.gas->a + 0.5 * gm1 * (QL.vel.x - QIF.vel.x);
	/* Now use the isentropic relations to get the rest. */
	Tratio = (2.0 + gm1 * (QL.vel.x / QL.gas->a) * (QL.vel.x / QL.gas->a)) /
	    (2.0 + gm1 * (QIF.vel.x / QIF.gas->a) * (QIF.vel.x / QIF.gas->a));
	QIF.gas->rho = QL.gas->rho * pow(Tratio, 1.0 / gm1);
	QIF.gas->p = QL.gas->p * pow(Tratio, geff / gm1);
	QIF.gas->e[0] = QL.gas->e[0] * Tratio;
	QIF.gas->T[0] = QL.gas->T[0] * Tratio;
#       else
	/* Take the easy way out and linearly interpolate. */
	UNUSED_VARIABLE(Tratio);
	QIF.gas->a = QL.gas->a - frac * (QL.gas->a - QLstar.gas->a);
	QIF.gas->rho = QL.gas->rho - frac * (QL.gas->rho - QLstar.gas->rho);
	QIF.gas->p = QL.gas->p - frac * (QL.gas->p - QLstar.gas->p);
	QIF.gas->e[0] = QL.gas->e[0] - frac * (QL.gas->e[0] - QLstar.gas->e[0]);
	QIF.gas->T[0] = QL.gas->T[0] - frac * (QL.gas->T[0] - QLstar.gas->T[0]);
#       endif
    }

    /* For the options below ...
     * The right wave and the contact discontinuity determine
     * the cell interface quantities. */

    if (option == 11 || option == 13) {
	/* All waves have gone into the left cell.
	 * The values are taken from the right cell only. */
	QIF.gas->rho = QR.gas->rho;
	QIF.vel.x = QR.vel.x;
	QIF.gas->p = QR.gas->p;
	QIF.gas->e[0] = QR.gas->e[0];
	QIF.gas->a = QR.gas->a;
	QIF.gas->T[0] = QR.gas->T[0];
    }

    if (option == 12 || option == 15) {
	/* The values are taken from behind the right wave. */
	QIF.gas->rho = QRstar.gas->rho;
	QIF.vel.x = QRstar.vel.x;
	QIF.gas->p = QRstar.gas->p;
	QIF.gas->e[0] = QRstar.gas->e[0];
	QIF.gas->a = QRstar.gas->a;
	QIF.gas->T[0] = QRstar.gas->T[0];
    }

    if (option == 14) {
	/* The right wave is an expansion. */
	VA = QR.vel.x + QR.gas->a;
	VB = QRstar.vel.x + QRstar.gas->a;
	/* The right rarefaction straddles the cell interface.
	 * Interpolate velocity inside the right rarefaction. */
	frac = (-VB) / (VA - VB);
	QIF.vel.x = QRstar.vel.x + frac * (QR.vel.x - QRstar.vel.x);
#       if (LINEAR != 1)
	/* Get local speed of sound by following a (u-a)
	 * characteristic line through the fan. */
	gm1 = geff - 1.0;
	QIF.gas->a = QR.gas->a + 0.5 * gm1 * (QIF.vel.x - QR.vel.x);
	/* Now use the isentropic relations to get the rest. */
	Tratio = (2.0 + gm1 * (QR.vel.x / QR.gas->a) * (QR.vel.x / QR.gas->a)) /
	    (2.0 + gm1 * (QIF.vel.x / QIF.gas->a) * (QIF.vel.x / QIF.gas->a));
	QIF.gas->rho = QR.gas->rho * pow(Tratio, 1.0 / gm1);
	QIF.gas->p = QR.gas->p * pow(Tratio, geff / gm1);
	QIF.gas->e[0] = QR.gas->e[0] * Tratio;
	QIF.gas->T[0] = QR.gas->T[0] * Tratio;
#       else
	/* Take the easy way out and linearly interpolate. */
	QIF.gas->a = QRstar.gas->a + frac * (QR.gas->a - QRstar.gas->a);
	QIF.gas->rho = QRstar.gas->rho + frac * (QR.gas->rho - QRstar.gas->rho);
	QIF.gas->p = QRstar.gas->p + frac * (QR.gas->p - QRstar.gas->p);
	QIF.gas->e[0] = QRstar.gas->e[0] + frac * (QR.gas->e[0] - QRstar.gas->e[0]);
	QIF.gas->T[0] = QRstar.gas->T[0] + frac * (QR.gas->T[0] - QRstar.gas->T[0]);
#       endif
    }

    /* ******************
     * Passive Quantities.
     * ******************
     *
     * We assume that the transverse velocity is unaffected by
     * the normal interactions.  We only need to select the
     * correct value.  This is assumed so for species mass
     * fraction also.
     */
    if (QIF.vel.x < 0.0) {
	QIF.vel.y = QR.vel.y;
	QIF.vel.z = QR.vel.z;
    } else {
	QIF.vel.y = QL.vel.y;
	QIF.vel.z = QL.vel.z;
    }
    for ( int jspec = 0; jspec < nsp; ++jspec ) {
	if (QIF.vel.x < 0.0)
	    QIF.gas->massf[jspec] = QR.gas->massf[jspec];
	else
	    QIF.gas->massf[jspec] = QL.gas->massf[jspec];
    }
    for ( int imode = 0; imode < nmodes; ++imode ) {
	if (QIF.vel.x < 0.0)
	    QIF.gas->e[imode] = QR.gas->e[imode];
	else
	    QIF.gas->e[imode] = QL.gas->e[imode];
    }
    return SUCCESS;
}   /* end of rivp_stage_3() */
