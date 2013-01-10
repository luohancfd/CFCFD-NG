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
#define ADJACENT        0
#define COMMON          0
#define SUP_IN          1
#define SUP_OUT         2
#define EXTRAPOLATE_OUT 2
#define SLIP            3
#define SLIP_WALL       3
#define ADIABATIC       4
#define FIXED_T         5
#define SUBSONIC_IN     6
#define SUBSONIC_OUT    7
#define TRANSIENT_UNI   8
#define TRANSIENT_PROF  9
#define STATIC_PROF    10
#define FIXED_P_OUT    11
#define RRM            12
#define TRANSIENT_T_WALL 13
// WAS #define EULER_MANUFACTURED 14
#define SEB            15
#define USER_DEFINED   16
#define ADJACENT_PLUS_UDF 17
#define ABLATING       18
#define SLIDING_T      19
#define FSTC           20
#define SHOCK_FITTING_IN 21

#define SPECIAL        -1

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

#define NON_CATALYTIC       21
#define EQUIL_CATALYTIC     22
#define SUPER_CATALYTIC     23

/* Before leaving, leave a mark... */
#define BC_DEFS_HH
#endif
