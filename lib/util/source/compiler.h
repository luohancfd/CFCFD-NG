/** \file compiler.h
 * \brief Compiler Specific definitions for C code.
 * \ingroup util
 *
 * This header file contains a few definitions for convenience
 * when porting the CFCFD codes across different machines and compilers.
 * Most were to get around pre-ANSI C definitions, 16-bit and 32-bit systems,
 * and to cope with early vectorising compilers that had trouble with 
 * pragmas from others.  
 *
 * \author PA Jacobs
 *
 * \version 07-mar-95 : Added the IBM Cset++ compiler as an option
 * \version 08-apr-95 : Added the Borland C++ compiler for OS/2
 * \version 23-Nov-95 : Added the SGI Power Challenge and the IBM SP2
 *                      MAYBE_STATIC for small stack computers such as PC's
 * \version 16-Jan-96 : EMX_GCC environment added (replaces old GNU_C)
 * \version 05-Dec-97 : On UNIX workstations, have only two types of compilers
 *                      KR_C   : the old K&R style C
 *                      ANSI_C : the newer standard
 * \version 29-Mar-98 : Move to Linux for future development; assume ANSI
 *                      unless specifically stated
 */

#ifndef COMPILER_H
#define COMPILER_H 1

/*
 * Select (ONLY) one of the following switches.
 * They will enable/disable function prototypes and the
 * inclusion of the standard library header file.
 */

/* MS-DOS systems... */
#define  TURBO_C     0
#define  MSOFT_C     0

/* OS/2 2.x, Warp ... */
#define  EMX_GCC     0
#define  TOPSPEED_C  0
#define  IBMCSET     0
#define  BC_OS2      0

/* Workstations... */
#define  KR_C        0
#define  ANSI_C      1

/* Real computers with vector processors... */
#define  CRAY_C      0
#define  VP_2600_C   0

/* Parallel computers */
#define  SGI_CHALLENGE  0
#define  IBM_SP2        0

/*-------------------------------------------------------------*/

/*
 * If we are on a vector-processing machine
 * and we would like to ignore some of the vector dependencies,
 * set the following flags to 1.
 * It makes the code about 10 times faster.
 */
#if CRAY_C == 1
#  define  CRAY_VECTOR 1
#else
#  define  CRAY_VECTOR 0
#endif

#if VP_2600_C == 1
#  define  VP_VECTOR  1
#else
#  define  VP_VECTOR  0
#endif

/*---------------------------------------------------------------------*/

/*
 * Full function prototypes or the old K&R style...
 */
#if  KR_C == 1
   /* Don't use function prototypes for old compilers. */
#  define  PROTO    0
#else
   /* For everything else, assume ANSI prototypes. */
#  define  PROTO    1
#endif

/*---------------------------------------------------------------------*/

/*
 * Standard library header (ANSII) ...
 */
#if KR_C == 1 
#   define  STDLIBH  0
#else
#   define  STDLIBH  1
#endif


/*
 * Console I/O header ...
 */
#if TURBO_C || MSOFT_C || TOPSPEED_C
#   define  CONIOH  1
#else
#   define  CONIOH  0
#endif


/*
 * Huge allocation functions, header ...
 */
#if  TURBO_C || MSOFT_C || TOPSPEED_C
#   define  ALLOCH  1
#else
#   define  ALLOCH  0
#endif

/*---------------------------------------------------------------------*/

/*
 * Set up a constant indicating the default word-length of the
 * compiler/machine.
 */
#if  TURBO_C || MSOFT_C || TOPSPEED_C
#  define  WORD_LENGTH  16
#else
#  define  WORD_LENGTH  32
#endif

/*-----------------------------------------------------------------*/

/*
 * Memory allocation for arrays and data structures
 * larger then 64Kbytes.
 */
#if  TURBO_C
   /* 16-bit compiler. */
#  define  HALLOC     farmalloc
#  define  HFREE      farfree
#  define  FAR_PTR    far
#  define  HUGE_PTR   huge
#elif  TOPSPEED_C || MSOFT_C
   /* 16-bit compiler. */
#  define  HALLOC     halloc
#  define  HFREE      hfree
#  define  FAR_PTR    far
#  define  HUGE_PTR   huge
#else
   /* 32-bit compilers don't need anything special. */
#  define  HALLOC     malloc
#  define  HFREE      free
#  define  FAR_PTR
#  define  HUGE_PTR
#endif

/*--------------------------------------------------------------*/

/*
 * The following MACROs will be 32-bit integers (at least)
 */
#if  TURBO_C || MSOFT_C || TOPSPEED_C
   /* 16-bit environment. */
#  define  INTEGER  long int
#  define  CARDINAL unsigned long
#else
   /* 32-bit (or greater) environment. */
#  define  INTEGER  int
#  define  CARDINAL unsigned
#endif

/*----------------------------------------------------------------*/

/*
 * Some (small) computers have a small stack that cannot hold large
 * (local) work arrays.  The following macro allows those environments
 * to declare the work arrays as static while allowing the larger
 * computers to have the work arrays on the stack.  This will most
 * probably affect only parallel codes as separate threads will need
 * to have separate workspaces.
 */

#if TURBO_C || MSOFT_C || TOPSPEED_C || EMX_GCC || IBMCSET || BC_OS2
#  define  MAYBE_STATIC   static
#else
#  define  MAYBE_STATIC
#endif

/*----------------------------------------------------------------*/
#endif


