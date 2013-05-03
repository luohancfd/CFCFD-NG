/** \file bc_defs.hh
 * \ingroup eilmer3
 * \brief Header boundary-condition index values.
 */


#ifndef BC_DEFS_HH
/* Go ahead and make some definitions... */

/** \brief Types of boundary conditions for blocks...
 *
 * \verbatim
 * ADJACENT   adjacent to another tile
 * SUP_IN     supersonic inflow
 * SUP_OUT    supersonic outflow
 * SLIP       slip/tangency (adiabatic)
 * ADIABATIC  no-slip, adiabatic
 * FIXED_T    no-slip, fixed T wall
 * SUBSONIC_IN  specify total pressure and temperature
 * SUBSONIC_OUT no-reflection BC for weak waves
 * TRANSIENT_UNI transient, uniform inflow (presumably)
 * TRANSIENT_PROF transient, profiled inflow
 * STATIC_PROF a profile is read from "profile.dat".
 * EULER_MANUFACTURED used in a verification case based on
 *    the Method of Manufactured Solutions
 * USER_DEFINED calls out to a Lua function to get the flow data
 * ADJACENT_PLUS_UDF exchanges data with another block then
 *    applies a user-defined Lua function
 *
 * SPECIAL    special purpose code has been included in
 *            the routines apply_inviscid_bc() and
 *            apply_viscous_bc().
 * \endverbatim
 */
const int ADJACENT = 0;
const int COMMON = 0;
const int SUP_IN = 1;
const int SUP_OUT = 2;
const int EXTRAPOLATE_OUT = 2;
const int SLIP = 3;
const int SLIP_WALL = 3;
const int ADIABATIC = 4;
const int FIXED_T = 5;
const int SUBSONIC_IN = 6;
const int SUBSONIC_OUT = 7;
const int TRANSIENT_UNI = 8;
const int TRANSIENT_PROF = 9;
const int STATIC_PROF = 10;
const int FIXED_P_OUT = 11;
const int RRM = 12;
const int TRANSIENT_T_WALL = 13;
// WAS const int EULER_MANUFACTURED = 14;
const int SEB = 15;
const int USER_DEFINED = 16;
const int ADJACENT_PLUS_UDF = 17;
const int ABLATING = 18;
const int SLIDING_T = 19;
const int FSTC = 20;
const int SHOCK_FITTING_IN = 21;

const int SPECIAL = -1;

/** \brief Types of wall catalycity available as boundary
 *         conditions for diffusive terms.
 *
 * \verbatim
 * NON_CATALYTIC         - as the name suggests (zero concentration gradient)
 * PARTIALLY_CATALYTIC   - some finite amount of recombination allowed
 * CATALYTIC             - equilibrium conditions at local p and and wall T enforced
 * FULLY_CATALYTIC       - same as CATALYTIC
 * SUPER_CATALYTIC       - free-stream conditions enforced
 * 
 * \endverbatim
 **/

const int NON_CATALYTIC = 22;
const int EQUIL_CATALYTIC = 23;
const int SUPER_CATALYTIC = 24;
const int PARTIALLY_CATALYTIC = 25;

/* Before leaving, leave a mark... */
#define BC_DEFS_HH
#endif
