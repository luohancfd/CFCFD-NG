/** \file useful.h
 * \brief Useful macro definitions for C code.
 * \ingroup util
 *
 * This header file contains a few definitions for convenience
 *
 * \author PJ
 *
 * \version 05-Aug-04 refactored from mb_cns.h, compiler.h, etc.
 */

#ifndef USEFUL_HEADER_ALREADY_INCLUDED
/* Go ahead and make some definitions... */


/** \brief Macro to eliminate warnings about unused variables.
 * \verbatim
 * from
 * Stephen J. Friedl
 * COMPILE TIME: Wizardly Words for Warnings
 * Linux Magazine August 2003
 * \endverbatim
 */
#define UNUSED_VARIABLE(x) ((void) (x))

/** Other Useful macros.
 */
#define FABS(a)        ( ((a) < 0.0) ? (-(a)) : (a) )
#define MAXIMUM(a,b)   ( ((a) > (b)) ? (a) : (b) )
#define MINIMUM(a,b)   ( ((a) < (b)) ? (a) : (b) )
#define SMALLEST(a,b)  ( (FABS(a) < FABS(b)) ? (a) : (b) )
#define LARGEST(a,b)   ( (FABS(a) > FABS(b)) ? (a) : (b) )

/* The following macro is to try to eliminate very small
 * floating-point values that seem to upset the Paraview program
 * when reading legacy VTK files.  It effectively sets the value
 * to zero when it is already very small.  Typically, it will be 
 * applied to velocities.
 */
#define UFLOWZ(x) ((fabs(x) < 1.0e-20) ? 0.0 : (x))

#ifdef BCC32
#   define inline
#endif

#define FAILURE                 1     /* Generic undefined failure */
#define SUCCESS                 0     /* Normal code exit. */
#define FILE_ERROR              102   /* Problems related to accessing files.*/
#define MEMORY_ERROR            103   /* Memory allocation problems. */
#define BAD_INPUT_ERROR         104   /* Errors of input from config files. */
#define NUMERICAL_ERROR         105   /* Irrecoverable numerical problems. */
#define NOT_INITIALISED_ERROR   106   /* Gas models not correctly initialised */
#define NOT_INITIALIZED_ERROR   106   /*     different spelling */
#define NULL_POINTER_ERROR      107   /* Someone gave us a NULL pointer. */
#define DUFF_EOS_ERROR          108   /* Something went wrong in call to EOS fn */
#define ZERO_DETERMINANT_ERROR  109
#define ITERATION_ERROR         110
#define BAD_TEMPERATURE_ERROR   111
#define MASS_FRACTION_ERROR     112   /* Typically mass fractions don't sum to 1.0 */
#define NOT_IMPLEMENTED_ERROR   113   /* Unfinished code, often in boundary condition functions */
#define VALUE_ERROR             114
#define BAD_CELLS_ERROR         115
#define DT_SEARCH_FAILED        116   /* The CFL check failed to find a good dt. */
#define BAD_TRANSPORT_ERROR     117
#define MISMATCHED_DIMENSIONS   118
#define INI_FILE_ERROR          119   /* Failed to load the dictionary from the INI file */
#define DUFF_UPDATE_ERROR       120
#define EOF_INDICATOR           121
#define PISTON_FAILURE          122
#define INDEX_OUT_OF_RANGE      123
#define LUA_ERROR               124   /* Associated with errors based on interaction with
				       * with the Lua interpreter. */
#define RECONSTRUCTION_ERROR    125   /* For problems with interpolation/reconstruction of
					 flow quantities. */
#define BAD_REACTION_RATE_ERROR 125
#define UDF_ERROR               126

#define USEFUL_HEADER_ALREADY_INCLUDED
#endif
/*-----------------------------------------------------------------*/


