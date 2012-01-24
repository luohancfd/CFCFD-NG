/** \file bspatch.h
 * \brief  Header file for the interpolation codes
 *
 * - used to find x,y,z coordinates of parametric
 *   points within an inlet
 *
 * \author Chris Craddock
 * \date   8/9/93
 *
 * \version Adapted for the "sm_3d" Space-marching flow solver
 *          P. J. 29-Sep-93
 * \version 11-Mar-95 : Debug added to bs_calc() and bs_calc2().
 *          Detection of the include files added.
 */

/*------------------------------------------------------------*/

#ifndef BSPATCH_H
#define BSPATCH_H 1

#ifndef COMPILER_H
#  include "../../util/source/compiler.h"
#endif

#ifndef GEOM_H
#  include "../../plot/geom.h"
#endif

/* Set to following macro to 1 to activate debug statements in bspatch.c */
#define  DEBUG_BSPATCH  0


int bs_init_arrays ( char base_file_name[]);

int bs_calc ( int neta, int nzeta,
              double xi, double eta[], double zeta[],
              struct point_3D ed1[],
              struct point_3D ed2[],
              struct point_3D ed3[],
              struct point_3D ed4[]);

int bs_calc2 (struct point_3D **q,
             double u, double v,
             int nx, int ny,
             struct point_3D *pnt );

double bs_nval ( double x, double t );

#endif
