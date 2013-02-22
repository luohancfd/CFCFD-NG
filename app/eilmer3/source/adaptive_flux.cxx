/** \file adaptive_flux.cxx
 * \ingroup eilmer3
 * 
 * \brief Adaptive flux calculator as suggested by James Quirk.
 *
 * \author PJ
 *
 * \version 13-Feb-01 : Finally pulled the finger out and implemented
 *                      a scheme combining EFM and AUSMDV.
 * \version 2010 : tweaks to suppress in the presence of significant shear.
 */

#include <stdio.h>
#include <math.h>
#include "../../../lib/util/source/useful.h"
#include "kernel.hh"
#include "cell.hh"
#include "flux_calc.hh"


/// \brief Compute the fluxes across an interface.
///
/// This adaptive flux calculator uses uses the Equilibrium Flux Method.
/// near shocks and AUSMDV away from shocks, however, we really don't want
/// EFM to be used across interfaces with strong shear.
/// EFM should still be able to do it's work as we really needed it for the
/// situations where the shock is closely aligned with the grid.
/// In that situation, we don't expect a stong shear at the interface.
///
/// The actual work is passed off to the original flux calculation functions.
///
/// \param Lft    : IN     : reference to Left flow state
///      (with velocities in local frame of reference)
/// \param Rght   : IN     : reference to Right flow state
/// \param IFace  : IN-OUT : reference to interface flux data structure
int adaptive_flux( FlowState &Lft, FlowState &Rght, FV_Interface &IFace )
{
    // if ( get_shock_fitting_flag() ) {
    // 	cerr << "Error, we have not implemented the ADAPTIVE flux calculator with shock fitting. Please use AUSMDV." << endl;
    // 	exit(NOT_IMPLEMENTED_ERROR);
    // }
    double sound_speed = 0.5 * (Lft.gas->a + Rght.gas->a);
    double shear_y = fabs(Lft.vel.y - Rght.vel.y) / sound_speed;
    double shear_z = fabs(Lft.vel.z - Rght.vel.z) / sound_speed;
    bool shear_is_small = MAXIMUM(shear_y, shear_z) <= get_shear_tolerance();
    if ( (Lft.S == 1 || Rght.S == 1) && shear_is_small ) {
	efmflx(Lft, Rght, IFace);
    } else {
	ausmdv(Lft, Rght, IFace);
    }
    return SUCCESS;
} // end adaptive_flux()
