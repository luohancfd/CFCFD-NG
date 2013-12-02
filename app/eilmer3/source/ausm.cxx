/** \file ausm.cxx
 * \ingroup eilmer3
 *
 * \brief One-dimensional flux splitting according to Liou and Steffen.
 *
 * \todo Really should get rid of the vector loop some time.
 */

/*------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include "flux_calc.hh"
#include "kernel.hh"

/*------------------------------------------------------------*/

/** \brief Compute interface flux of mass, momentum and energy.
 *
 * Given the initial LEFT and RIGHT states (QL and QR)
 * either side of an interface, compute the interface state.
 * This solution is returned as values of
 * the flow quantities at the interface (QIF) and the
 * velocities of the left and right waves.
 *
 * \param QL      : IN : reference to the LEFT state
 *       see "flux_calc.hh" for a definition of the components
 * \param QR      : IN : RIGHT state
 * \param QIF     : OUT : Flow state at the interface
 * \param WSL     : OUT : Left wave velocity
 * \param WSR     : OUT : Right wave velocity
 *
 * \version 1.0,  20 Aug 91
 * \version 1.1,  18-Oct-96, second-order polynomial for pressure splitting
 *                  as per Liou and Steffen's J. Comput. Phys paper
 *
 * \author P.A. Jacobs
 * ICASE
 * Mail Stop 123C
 * NASA Langley Rearch Centre
 * Hampton VA 23665.
 *
 * \verbatim
 * References ...
 * M. -S. Liou and C. J. Steffen Jr.
 * A new flux splitting scheme.
 * NASA Technical Memorandum 104404
 * \endverbatim
 */
int ausm(FlowState &QL, FlowState &QR, FlowState &QIF, double &WSL, double &WSR)
{
    global_data &G = *get_global_data_ptr();
    if ( G.shock_fitting ) {
    	cerr << "Error, we have not implemented AUSM with shock fitting. Please use AUSMDV." << endl;
    	exit(NOT_IMPLEMENTED_ERROR);
    }
    Gas_model *gmodel = get_gas_model_ptr();

    /* Choose the polynomial form for pressure splitting by selecting one... */
#   define FIRST_ORDER  1
#   define SECOND_ORDER 2
#   define P_SPLIT FIRST_ORDER

    double ML, MR, MLplus, MRminus, Mhalf, advect;
    double PLplus, PRminus, Phalf;

    /*
     * Characteristic speeds.
     */
    ML = QL.vel.x / QL.gas->a;
    MR = QR.vel.x / QR.gas->a;

    /*
     * Apply the van Leer splitting.
     */
    if (ML > 1.0)
	MLplus = ML;
    else if (ML > -1.0)
	MLplus = 0.25 * (ML + 1.0) * (ML + 1.0);
    else
	MLplus = 0.0;

    if (MR > 1.0)
	MRminus = 0.0;
    else if (MR > -1.0)
	MRminus = -0.25 * (MR - 1.0) * (MR - 1.0);
    else
	MRminus = MR;

    /*
     * Advective velocity.
     */
    Mhalf = MLplus + MRminus;
    if (Mhalf >= 0.0)
	advect = QL.gas->a * Mhalf;
    else
	advect = QR.gas->a * Mhalf;

    /*
     * Rolf's mod is to limit the Mach numbers to 1
     * for the pressure splitting.
     */
    if (MR > 1.0) MR = 1.0;
    if (MR < -1.0) MR = -1.0;
    if (ML > 1.0) ML = 1.0;
    if (ML < -1.0) ML = -1.0;

#   if P_SPLIT == FIRST_ORDER
    /*
     * Split pressure with the linear expression.
     * Liou & Steffen, J. Comput Phys p26, eqn (8b)
     * Note: Bram seems to prefer the cubic expression.
     */
    PLplus = 0.5 * QL.gas->p * (1.0 + ML);
    PRminus = 0.5 * QR.gas->p * (1.0 - MR);
#   endif

#   if P_SPLIT == SECOND_ORDER
    /*
     * Split pressure with the second-order polynomial.
     * Liou & Steffen, J. Comput Phys p26, eqn (8a)
     */
    PLplus = 0.25 * QL.gas->p * (1.0 + ML) * (1.0 + ML) * (2.0 - ML);
    PRminus = 0.5 * QR.gas->p * (1.0 - MR) * (1.0 - MR) * (2.0 + MR);
#   endif

    Phalf = PLplus + PRminus;

    /*
     * Set the interface properties based on the advective
     * velocity. 
     * NOTE that the normal velocity saved here will result
     * in a slightly different flux calculation to that 
     * described in the original report.
     */
    QIF.vel.x = advect;
    QIF.gas->p = Phalf;

    if (advect >= 0.0) {
	/* Use Left state */
	QIF.gas->rho = QL.gas->rho;
	QIF.gas->a = QL.gas->a;
	for ( size_t itm=0; itm<QIF.gas->e.size(); ++itm ) {
	    QIF.gas->e[itm] = QL.gas->e[itm];
	    QIF.gas->T[itm] = QL.gas->T[itm];
	}
    } else {
	/* Use Right state. */
	QIF.gas->rho = QR.gas->rho;
	QIF.gas->a = QR.gas->a;
	for ( size_t itm=0; itm<QIF.gas->e.size(); ++itm ) {
	    QIF.gas->e[itm] = QR.gas->e[itm];
	    QIF.gas->T[itm] = QR.gas->T[itm];
	}
    }

    /*
     * Dummy wave-speeds.
     */
    WSL = QL.vel.x - QL.gas->a;
    WSR = QR.vel.x + QR.gas->a;

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
    int nsp = gmodel->get_number_of_species();
    for ( int jspec = 0; jspec < nsp; ++jspec ) {
	if (QIF.vel.x < 0.0)
	    QIF.gas->massf[jspec] = QR.gas->massf[jspec];
	else
	    QIF.gas->massf[jspec] = QL.gas->massf[jspec];
    }
    int nmodes = gmodel->get_number_of_modes();
    for ( int imode = 0; imode < nmodes; ++imode ) {
	if (QIF.vel.x < 0.0)
	    QIF.gas->e[imode] = QR.gas->e[imode];
	else
	    QIF.gas->e[imode] = QL.gas->e[imode];
    }
    
    return SUCCESS;
} /* end ausm() */
