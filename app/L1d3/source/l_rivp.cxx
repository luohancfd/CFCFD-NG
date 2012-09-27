/** \file l_rivp.cxx
 * \ingroup l1d3
 * \brief One-dimensional Riemann Initial Value Problem for l1d.cxx.
 *
 * \author PA Jacobs
 */

/*------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

#include "../../../lib/util/source/useful.h"
#include "../../../lib/nm/source/qd_power.h"
#include "l1d.hh"
#include "l_kernel.hh"
#include "l_cell.hh"

/*------------------------------------------------------------*/

int L_rivp(struct L_flow_state QL[], struct L_flow_state QR[],
           double ustar[], double pstar[],
           int first, int last, int end1, int end2)
{
    /*
     * Purpose ...
     * -------
     * Given the initial LEFT and RIGHT states (QL and QR)
     * either side of an interface, compute the solution to
     * the one-dimensional Riemann problem when the interface
     * is removed.  This solution is returned as values of
     * the interface pressures and velocities.
     *
     *
     * Input ...
     * -----
     * *QL     : pointer to the array of LEFT states
     * 	     see "physics.h" for a definition of the components
     * *QR     : RIGHT states
     * first   : first element to be computed
     * last    : last element to be computed
     * end1    : type of boundary condition at first interface
     *           0 : standard solution with L and R states
     *           1 : ustar and R state is specified, find pstar only
     * end2    : type of boundary condition at last interface
     *           0 : standard solution with L and R states
     *           1 : ustar and L state is specified, find pstar only
     * ustar   : the first and/or the last elements may be specified
     *           if the appropraite boundary condition flag is set
     *
     * Output ...
     * ------
     * pstar   : interface pressures
     * ustar   : interface velocities
     *
     */

    /*---------------------------------------------------------------*/

    int i, ip;

    static double geff[NDIM];   /* effective GAMMA          */
    double g, gL, gR;
    static double sqrL[NDIM], sqrR[NDIM];
    double alpha;
    double gm1, gp1, z, uLbar, uRbar;
    double base, expon, pwr;    /* for qd_power MACRO */
    double temporary;
    double term1, term2, F, dFdpstar;
    double PMIN, dp;
    Gas_model *gmodel = get_gas_model_ptr();

#   define  QUICK_AND_DIRTY  0

    /*************************************************************
     * Compute a Roe-averaged effective-gamma for each interface *
     * as used by Edwards (1988) and Vinokur & Liu (1988).       *
     *************************************************************/

    /* 
     * General (gaseous) equation of state.
     * Compute a Roe-averaged effective-gamma. 
     */
    for (i = first; i <= last; ++i) {
        if (end1 == 1 && i == first) {
            sqrR[i] = sqrt(QR[i].gas->rho);
            geff[i] = gmodel->gamma(*(QR[i].gas));
        } else if (end2 == 1 && i == last) {
            sqrL[i] = sqrt(QL[i].gas->rho);
            geff[i] = gmodel->gamma(*(QL[i].gas));
        } else {
            sqrL[i] = sqrt(QL[i].gas->rho);
            sqrR[i] = sqrt(QR[i].gas->rho);
            gL = gmodel->gamma(*(QL[i].gas));
            gR = gmodel->gamma(*(QR[i].gas));
            alpha = sqrL[i] / (sqrL[i] + sqrR[i]);
            geff[i] = alpha * gL + (1.0 - alpha) * gR;
        }
    }

    /*****************************************************************
     * Solve one-dimensional problem in a frame of reference         *
     * for which QL[i].u is the normal velocity.                     *
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
     * STAGE 3: OMITTED FOR THE LAGRANGIAN CASE                      *
     * Back out the other intermediate variables using either the    *
     * shock relations or the isentropic wave relations.             *
     *                                                               *
     *****************************************************************/

    for (i = first; i <= last; ++i) {

        g = geff[i];
        gm1 = g - 1.0;
        gp1 = g + 1.0;

        PMIN = 1.0e-6;

        if (end1 == 1 && i == first) {
            /*
             * ***************************
             * * Left boundary condition *
             * ***************************
             *
             * ustar is specified.
             */
#           if (DEBUG >= 3)
            printf("Left boundary rivp:\n");
#           endif
            uRbar = QR[i].u - 2.0 / gm1 * QR[i].gas->a;
            term1 = gm1 / (2.0 * sqrt(g));
            term2 = sqrt(QR[i].gas->rho / pow(QR[i].gas->p, 1.0 / g));
            base = (ustar[i] - uRbar) * term1 * term2;
            expon = 2.0 * g / gm1;
            pstar[i] = pow(base, expon);
#           if (DEBUG >= 3)
            printf("QR.u=%e, QR.a=%e, g=%e\n", QR[i].u, QR[i].gas->a, g);
            printf("uRbar=%e, term1=%e, term2=%e\n", uRbar, term1, term2);
            printf("base=%e, expon=%e\n", base, expon);
            printf("ustar=%e, pstar=%e\n", ustar[i], pstar[i]);
#           endif
        } else if (end2 == 1 && i == last) {
            /* 
             * ****************************
             * * Right boundary condition *
             * ****************************
             *
             * ustar is specified.
             */
#           if (DEBUG >= 3)
            printf("Right boundary rivp:\n");
#           endif
            uLbar = QL[i].u + 2.0 / gm1 * QL[i].gas->a;
            term1 = gm1 / (2.0 * sqrt(g));
            term2 = sqrt(QL[i].gas->rho / pow(QL[i].gas->p, 1.0 / g));
            base = (uLbar - ustar[i]) * term1 * term2;
            expon = 2.0 * g / gm1;
            pstar[i] = pow(base, expon);
#           if (DEBUG >= 3)
            printf("QL.u=%e, QL.a=%e, g=%e\n", QL[i].gas->u, QL[i].gas->a, g);
            printf("uLbar=%e, term1=%e, term2=%e\n", uLbar, term1, term2);
            printf("base=%e, expon=%e\n", base, expon);
            printf("ustar=%e, pstar=%e\n", ustar[i], pstar[i]);
#           endif
        } else {
            /* 
             * **********************************************
             * * Standard application of the Riemann solver * 
             * **********************************************
             *
             * We have valid Left and Right states.
             */

            /*
             * -----------------------------------------------------
             * STAGE 1: Explicit solution using two isentropic waves.
             * -----------------------------------------------------
             */
#           if (DEBUG >= 3)
            printf("STAGE 1\n");
#           endif

            /* 
             * Riemann invariants. 
             */
            uLbar = QL[i].u + 2.0 / gm1 * QL[i].gas->a;
            uRbar = QR[i].u - 2.0 / gm1 * QR[i].gas->a;

            if ((uLbar - uRbar) < 0.0) {
                /* 
                 * We have a situation in which a (near) vacuum is formed
                 * between the left and right states.
                 */

                ustar[i] = 0.0;
                pstar[i] = PMIN;

            }
            /* End of Vacuum solution */
            else {
                /*
                 * Positive-pressure solution.
                 */

                /* Intermediate variable. */
                base = QL[i].gas->p / QR[i].gas->p;
                expon = gm1 / (2.0 * g);
#               if (QUICK_AND_DIRTY == 1)
                /* We expect exponent to be <= 0.2 */
                if (base < 0.12 || base > 8.4) {
                    /* Use the full power function */
                    pwr = pow(base, expon);
                } else {
                    /* Use the quick and dirty approach */
                    QD_POWER(base, expon, pwr);
                }
#               else
                pwr = pow(base, expon);
#               endif
                z = (QR[i].gas->a / QL[i].gas->a) * pwr;

                /* First, compute ustar. */
                ustar[i] = (uLbar * z + uRbar) / (1.0 + z);

                /* Then, pstar. */
                base = (0.5 * gm1 * (uLbar - uRbar) / (QL[i].gas->a * (1.0 + z)));
                expon = 2.0 * g / gm1;
#               if (QUICK_AND_DIRTY == 1)
                /* Expect 5.0 <= expon <= 12 */
                if (base < 0.75 || base > 1.5) {
                    /* Use the full power function */
                    pwr = pow(base, expon);
                } else {
                    /* Reduce the exponent and use the fast macro */
                    expon -= 8.0;
                    QD_POWER(base, expon, pwr);
                    temporary = base * base;
                    temporary *= temporary;
                    temporary *= temporary;
                    pwr *= temporary;
                }
#               else
                /* Use the full power function */
                pwr = pow(base, expon);
		UNUSED_VARIABLE(temporary);
#               endif
                pstar[i] = QL[i].gas->p * pwr;

                /*
                 * ------------------------------------------
                 * End of isentropic-wave (explicit) solution.
                 * ------------------------------------------
                 */
#               if (DEBUG >= 3)
                printf("Isentropic wave solution: pstar = %e, ustar = %e\n",
                       pstar[i], ustar[i]);
#               endif


                /*
                 * ---------------------------------------------------
                 * STAGE 2: Apply the strong-shock relations if needed.
                 * ---------------------------------------------------
                 */
#               define  BIG_RATIO  10.0
#               define MAX_NEWTON_STEPS 12
#               define NEWTON_TOL 0.001
#               define TOO_MUCH_ERROR 0.1

                if (pstar[i] > BIG_RATIO * QL[i].gas->p
                    && pstar[i] > BIG_RATIO * QR[i].gas->p) {
                    /*
                     * Apply the STRONG SHOCK relations to get an explicit solution
                     * if both of the pressure jumps are large enough.
                     */
                    ustar[i] = (sqrL[i] * QL[i].u + sqrR[i] * QR[i].u) /
                        (sqrL[i] + sqrR[i]);
                    term1 =
                        sqrR[i] * (QL[i].u - QR[i].u) / (sqrL[i] + sqrR[i]);
                    pstar[i] = 0.5 * gp1 * QL[i].gas->rho * term1 * term1;
#                   if (DEBUG >= 3)
                    printf
                        ("Two strong shocks: pstar = %e, ustar = %e, gamma = %f\n",
                         pstar[i], ustar[i], g);
                    printf("   QL: u=%e, rho=%e, p=%e, T=%e\n", QL[i].u,
                           QL[i].gas->rho, QL[i].gas->p, QL[i].gas->T);
                    printf("   QR: u=%e, rho=%e, p=%e, T=%e\n", QR[i].u,
                           QR[i].gas->rho, QR[i].gas->p, QR[i].gas->T);
#                   endif
                } else if (pstar[i] > BIG_RATIO * QR[i].gas->p) {
                    /*
                     * Treat the right-moving wave as a strong shock, the
                     * left-moving wave as an isentropic wave, and take
                     * several Newton steps to improve the guess for pstar, ustar.
                     */
                    for (ip = 0; ip < MAX_NEWTON_STEPS; ++ip) {
                        term1 = QL[i].gas->a * 
			    pow(pstar[i] / QL[i].gas->p, gm1 / (2.0 * g));
                        term2 = sqrt(2.0 * pstar[i] / (QR[i].gas->rho * gp1));
                        F = uLbar - 2.0 * term1 / gm1 - QR[i].u - term2;
                        dFdpstar = -term1 / (g * pstar[i]) -
                            1.0 / (QR[i].gas->rho * gp1 * term2);
                        dp = -F / dFdpstar;
                        if (pstar[i] + dp > PMIN) {
                            pstar[i] = pstar[i] + dp;
                        } else {
                            pstar[i] *= 0.1;
                        }   /* end if */
			if (fabs(dp) / pstar[i] < NEWTON_TOL) break;
                    }   /* end for */

                    ustar[i] =
                        QR[i].u + sqrt(2.0 * pstar[i] / (QR[i].gas->rho * gp1));

                    if (fabs(dp) / pstar[i] > TOO_MUCH_ERROR) {
                        printf( "Right shock wave: Newton iteration is crap: \n" );
                        printf( "pstar=%e, dp=%e, ustar=%e\n",
				pstar[i], dp, ustar[i]);
                    }   /* end if */
#                   if (DEBUG >= 3)
                    printf("Right strong shock: pstar = %e, ustar = %e\n",
                           pstar[i], ustar[i]);
#                   endif
                } else if (pstar[i] > BIG_RATIO * QL[i].gas->p) {
                    /*
                     * Treat the left-moving wave as a strong shock, the
                     * right-moving wave as an isentropic wave, and take
                     * several Newton steps to improve the guess for pstar, ustar.
                     */
                    for (ip = 0; ip < MAX_NEWTON_STEPS; ++ip) {
                        term1 =
                            QR[i].gas->a * pow(pstar[i] / QR[i].gas->p,
                                          gm1 / (2.0 * g));
                        term2 = sqrt(2.0 * pstar[i] / (QL[i].gas->rho * gp1));
                        F = -uRbar - 2.0 * term1 / gm1 + QL[i].u - term2;
                        dFdpstar = -term1 / (g * pstar[i]) -
                            1.0 / (QL[i].gas->rho * gp1 * term2);
                        dp = -F / dFdpstar;
                        if (pstar[i] + dp > PMIN) {
                            pstar[i] = pstar[i] + dp;
                        } else {
                            pstar[i] *= 0.1;
                        }   /* end if */
			if (fabs(dp) / pstar[i] < NEWTON_TOL) break;
                    }   /* end for */

                    ustar[i] =
                        QL[i].u - sqrt(2.0 * pstar[i] / (QL[i].gas->rho * gp1));

                    if (fabs(dp) / pstar[i] > TOO_MUCH_ERROR) {
                        printf("Left shock wave: Newton iteration is crap: \n");
                        printf("i=%d pstar=%e, dp=%e, ustar=%e\n",
                             i, pstar[i], dp, ustar[i]);
                    }   /* end if */
#                   if (DEBUG >= 3)
                    printf("Left strong shock: pstar = %e, ustar = %e\n",
                           pstar[i], ustar[i]);
#                   endif
                }

            }   /* End of positive-pressure solution */

        }   /* end of standard application of the Riemann solver */

	/* 
	 * --------------------------------------------------------------
	 * At this point, we should have a good estimate for pstar, ustar.
	 * --------------------------------------------------------------
	 */
	if ( pstar[i] > 1.0e20 || pstar[i] < 0.0 || ustar[i] > 1.0e20 ) {
	    printf("L_rivp(): unreasonable value for pressure\n");
	    printf("    i=%d pstar=%g ustar=%g first=%d last=%d end1=%d end2=%d\n", i, pstar[i], ustar[i], first, last, end1, end2);
	    printf("    QL u=%g\n", QL[i].u);
	    QL[i].gas->print_values();
	    printf("    QR u=%g\n", QR[i].u);
	    QR[i].gas->print_values();
	    printf("    geff=%g, uLbar=%g uRbar=%g\n", geff[i], uLbar, uRbar);
	}
    } // end of i loop

    return SUCCESS;
} // end L_rivp()

