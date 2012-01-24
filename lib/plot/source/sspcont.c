/** \file sspcont.c
 * \ingroup plot
 * \brief Contouring routines for Simple Scientific Plotting
 *
 * \author PA Jacobs
 *
 * \version 1.0   : ??-Oct-91
 * \version 1.1   : 04-Jun-93 : Added ssp_ContourTriangle3D()
 * \version 1.2   : 08-Jun-93 : Added ssp_ContourQuad3D()
 *                  Improved the triangle contouring by adding special
 *                  cases if the contour is along one side.
 *                  Don't draw the triangle if all vertex values lie
 *                  on the level value.
 * \version 1.3   : 16-Aug-93 : Tighten the tolerance from 1.0e-15 to 1.0e-20
 *                  to see if we can get fewer gaps. -- didn't help.
 *                  The problem seems to have been in not taking the
 *                  absolute values of delA, delB and delC when
 *                  checking and scaling the offsets.
 * \version 1.4   : 01-Jan-94 : Colour-filled contours.
 * \version 1.5   : 12-Jun-94 : Grey-shaded contours as well as colour.
 *
 */

/*
 * ------------------
 * Global Definitions
 * ------------------
 */

#include "../../util/source/compiler.h"
#include "../../geometry/source/geom.h"
#include "ssp.h"

#include <stdio.h>
#include <math.h>

#if (STDLIBH)
#  include <stdlib.h>
#  include <stdarg.h>
#endif

#define  CON_TOL  1.0e-20

/*-------------------------------------------------------------*/

#if (PROTO)

int ssp_ContourTriangle (struct ssp_data_point *A,
		     struct ssp_data_point *B,
		     struct ssp_data_point *C,
		     double level)

#else

int ssp_ContourTriangle (A, B, C, level)
struct ssp_data_point *A, *B, *C;
double level;

#endif

/*
 * Purpose...
 * -------
 * Given the three vertices of a triangle, assume a linear
 * distribution and plot at the specified level.
 *
 * Input...
 * -----
 * A, B, C  : pointers to the corner data structures
 *            (Assume that they are in a counter-clockwise order.)
 * level    : value of the function at which the contour
 *            is to be drawn
 *
 *                     B -------- A
 *                      \       /
 *                        \   /
 *                          C
 *
 */

{  /* begin ssp_ContourTriangle () */
double delA, delB, delC, delmax;
double multAB, multBC, multCA;
double lambdaAB, lambdaBC, lambdaCA;
double xAB, yAB, xBC, yBC, xCA, yCA;

/*
 * Work in terms of the differences to the specified level.
 */
delA = A->value - level;
delB = B->value - level;
delC = C->value - level;

/* Suppresses a maybe uninitialized warning -- should fix */
xAB = 0.0; yAB = 0.0; xBC = 0.0; yBC = 0.0; xCA = 0.0; yCA = 0.0;

/*
 * Pick the maximum magnitude and scale if it is not zero.
 */
delmax = fabs(delA);
if (fabs(delB) > delmax) delmax = fabs(delB);
if (fabs(delC) > delmax) delmax = fabs(delC);
if (delmax > 0.0) {
   delA /= delmax;
   delB /= delmax;
   delC /= delmax;
}

if (delmax < CON_TOL) {
   /*
    * All points lie on the level; plot nothing.
    */
   return (0);
}

/*
 * Check for a complete miss.
 */
if (delA > 0.0 && delB > 0.0 && delC > 0.0) return (0);
if (delA < 0.0 && delB < 0.0 && delC < 0.0) return (0);

/*
 * Special cases where the contour is along one edge.
 */
if (fabs(delA) < CON_TOL && fabs(delB) < CON_TOL) {
   /* The contour lies along AB, plot it and leave. */
   ssp_Move (A->x, A->y);
   ssp_Plot (B->x, B->y);
   return (0);
}

if (fabs(delB) < CON_TOL && fabs(delC) < CON_TOL) {
   /* The contour lies along BC, plot it and leave. */
   ssp_Move (B->x, B->y);
   ssp_Plot (C->x, C->y);
   return (0);
}

if (fabs(delC) < CON_TOL && fabs(delA) < CON_TOL) {
   /* The contour lies along CA, plot it and leave. */
   ssp_Move (C->x, C->y);
   ssp_Plot (A->x, A->y);
   return (0);
}


/*
 * Locate intercepts along each side.
 */

multAB = delA * delB;
lambdaAB = -1.0;
if (multAB < 0.0) {
   lambdaAB = -delA / (delB - delA);
   xAB = A->x + lambdaAB * (B->x - A->x);
   yAB = A->y + lambdaAB * (B->y - A->y);
} else if (delA == 0.0) {
   lambdaAB = 0.0;
   xAB = A->x;
   yAB = A->y;
}

multBC = delB * delC;
lambdaBC = -1.0;
if (multBC < 0.0) {
   lambdaBC = -delB / (delC - delB);
   xBC = B->x + lambdaBC * (C->x - B->x);
   yBC = B->y + lambdaBC * (C->y - B->y);
} else if (delB == 0.0) {
   lambdaBC = 0.0;
   xBC = B->x;
   yBC = B->y;
}

multCA = delC * delA;
lambdaCA = -1.0;
if (multCA < 0.0) {
   lambdaCA = -delC / (delA - delC);
   xCA = C->x + lambdaCA * (A->x - C->x);
   yCA = C->y + lambdaCA * (A->y - C->y);
} else if (delC == 0.0) {
   lambdaCA = 0.0;
   xCA = C->x;
   yCA = C->y;
}

/*
 * Now plot the contours (straight lines) connecting the
 * intercepts.
 */

if ( (lambdaAB >= 0.0 && lambdaAB <= 1.0) &&
     (lambdaBC >= 0.0 && lambdaBC <= 1.0) ) {
   ssp_Move (xAB, yAB);
   ssp_Plot (xBC, yBC);
}

if ( (lambdaBC >= 0.0 && lambdaBC <= 1.0) &&
     (lambdaCA >= 0.0 && lambdaCA <= 1.0) ) {
   ssp_Move (xBC, yBC);
   ssp_Plot (xCA, yCA);
}

if ( (lambdaCA >= 0.0 && lambdaCA <= 1.0) &&
     (lambdaAB >= 0.0 && lambdaAB <= 1.0) ) {
   ssp_Move (xCA, yCA);
   ssp_Plot (xAB, yAB);
}


return (0);
}  /* end of ssp_ContourTriangle () */

/*-------------------------------------------------------------*/

#if (PROTO)

int ssp_ContourTriangle3D (struct point_3D *A,
		       struct point_3D *B,
		       struct point_3D *C,
		       double A_value,
		       double B_value,
		       double C_value,
		       double level)

#else

int ssp_ContourTriangle3D (A, B, C, A_value, B_value, C_value, level)
struct point_3D *A, *B, *C;
double A_value, B_value, C_value;
double level;

#endif

/*
 * Purpose...
 * -------
 * Given the three vertices of a triangle, assume a linear
 * distribution and plot at the specified level.
 * <<< Three-Dimensional Version. >>>
 *
 * Input...
 * -----
 * A, B, C  : pointers to the corner vertices
 *            (Assume that they are in a counter-clockwise order.)
 * A_value,
 * B_value,
 * C_value  : Data values at the vertices
 * level    : value of the function at which the contour
 *            is to be drawn
 *
 *                     B -------- A
 *                      \       /
 *                        \   /
 *                          C
 *
 */

{  /* begin ssp_ContourTriangle3D () */
double delA, delB, delC, delmax;
double multAB, multBC, multCA;
double lambdaAB, lambdaBC, lambdaCA;
double xAB, yAB, zAB, xBC, yBC, zBC, xCA, yCA, zCA;

/*
 * Work in terms of the differences to the specified level.
 */
delA = A_value - level;
delB = B_value - level;
delC = C_value - level;

/* Suppresses a "maybe uninitialized" warning -- should fix */
xAB = 0.0; yAB = 0.0; zAB = 0.0;
xBC = 0.0; yBC = 0.0; zBC = 0.0;
xCA = 0.0; yCA = 0.0; zCA = 0.0;

/*
 * Pick the maximum magnitude and scale if it is not zero.
 */
delmax = fabs(delA);
if (fabs(delB) > delmax) delmax = fabs(delB);
if (fabs(delC) > delmax) delmax = fabs(delC);
if (delmax > 0.0) {
   delA /= delmax;
   delB /= delmax;
   delC /= delmax;
}

if (delmax == 0.0) {
   /*
    * All points lie on the level; plot nothing.
    */
   return (0);
}

/*
 * Check for a complete miss.
 */
if (delA > 0.0 && delB > 0.0 && delC > 0.0) return (0);
if (delA < 0.0 && delB < 0.0 && delC < 0.0) return (0);

/*
 * Special cases where the contour is along one edge.
 */
if (fabs(delA) < CON_TOL && fabs(delB) < CON_TOL) {
   /* The contour lies along AB, plot it and leave. */
   ssp_Move3D (A->x, A->y, A->z);
   ssp_Plot3D (B->x, B->y, B->z);
   return (0);
}

if (fabs(delB) < CON_TOL && fabs(delC) < CON_TOL) {
   /* The contour lies along BC, plot it and leave. */
   ssp_Move3D (B->x, B->y, B->z);
   ssp_Plot3D (C->x, C->y, C->z);
   return (0);
}

if (fabs(delC) < CON_TOL && fabs(delA) < CON_TOL) {
   /* The contour lies along CA, plot it and leave. */
   ssp_Move3D (C->x, C->y, C->z);
   ssp_Plot3D (A->x, A->y, A->z);
   return (0);
}


/*
 * Locate intercepts along each side.
 */

multAB = delA * delB;
lambdaAB = -1.0;
if (multAB < 0.0) {
   lambdaAB = -delA / (delB - delA);
   xAB = A->x + lambdaAB * (B->x - A->x);
   yAB = A->y + lambdaAB * (B->y - A->y);
   zAB = A->z + lambdaAB * (B->z - A->z);
} else if (delA == 0.0) {
   lambdaAB = 0.0;
   xAB = A->x;
   yAB = A->y;
   zAB = A->z;
}

multBC = delB * delC;
lambdaBC = -1.0;
if (multBC < 0.0) {
   lambdaBC = -delB / (delC - delB);
   xBC = B->x + lambdaBC * (C->x - B->x);
   yBC = B->y + lambdaBC * (C->y - B->y);
   zBC = B->z + lambdaBC * (C->z - B->z);
} else if (delB == 0.0) {
   lambdaBC = 0.0;
   xBC = B->x;
   yBC = B->y;
   zBC = B->z;
}

multCA = delC * delA;
lambdaCA = -1.0;
if (multCA < 0.0) {
   lambdaCA = -delC / (delA - delC);
   xCA = C->x + lambdaCA * (A->x - C->x);
   yCA = C->y + lambdaCA * (A->y - C->y);
   zCA = C->z + lambdaCA * (A->z - C->z);
} else if (delC == 0.0) {
   lambdaCA = 0.0;
   xCA = C->x;
   yCA = C->y;
   zCA = C->z;
}

/*
 * Now plot the contours (straight lines) connecting the
 * intercepts.
 */

if ( (lambdaAB >= 0.0 && lambdaAB <= 1.0) &&
     (lambdaBC >= 0.0 && lambdaBC <= 1.0) ) {
   ssp_Move3D (xAB, yAB, zAB);
   ssp_Plot3D (xBC, yBC, zBC);
}

if ( (lambdaBC >= 0.0 && lambdaBC <= 1.0) &&
     (lambdaCA >= 0.0 && lambdaCA <= 1.0) ) {
   ssp_Move3D (xBC, yBC, zBC);
   ssp_Plot3D (xCA, yCA, zCA);
}

if ( (lambdaCA >= 0.0 && lambdaCA <= 1.0) &&
     (lambdaAB >= 0.0 && lambdaAB <= 1.0) ) {
   ssp_Move3D (xCA, yCA, zCA);
   ssp_Plot3D (xAB, yAB, zAB);
}


return (0);
}  /* end of ssp_ContourTriangle3D () */

/*--------------------------------------------------------------*/

#if (PROTO)

int ssp_ContourQuad (struct ssp_data_point *A,
		     struct ssp_data_point *B,
		     struct ssp_data_point *C,
		     struct ssp_data_point *D,
		     double level)

#else

int ssp_ContourQuad (A, B, C, D, level)
struct ssp_data_point *A, *B, *C, *D;
double level;

#endif

/*
 * Purpose...
 * -------
 * Given the 4 corners of a quadrilateral, split it into
 * 8 triangles and contour each of the triangles.
 * This arrangement follows the scheme described by
 * James Quirk in his PhD thesis.
 *
 * Input...
 * -----
 * A, B, C, D : pointers to the corner data structures
 *              These must be specified in a counter-clockwise
 *              order.
 *              C--------B
 *              |        |
 *              |        |
 *              D--------A
 *
 * level     : value of the function at which the contour is to
 *             be drawn.
 *
 * Intermediate Points ...
 * -------------------
 *
 *            C --- MBC --- B
 *            |  \   |  /   |
 *           MCD -- Ctr -- MAB
 *            |  /   |  \   |
 *            D --- MDA --- A
 *
 * Revisions...
 * ---------
 * 30-Sep-92 : check to see whether all data points miss the
 * 	present contour before contouring each triangle.
 *	(Keith Weinman's suggestion)
 *
 */

{  /* begin ssp_ContourQuad () */
struct ssp_data_point MAB, MBC, MCD, MDA, Ctr;

/*
 * Set up midpoints along each side
 */
MAB.x = 0.5 * (A->x + B->x);
MAB.y = 0.5 * (A->y + B->y);
MAB.value = 0.5 * (A->value + B->value);

MBC.x = 0.5 * (B->x + C->x);
MBC.y = 0.5 * (B->y + C->y);
MBC.value = 0.5 * (B->value + C->value);

MCD.x = 0.5 * (C->x + D->x);
MCD.y = 0.5 * (C->y + D->y);
MCD.value = 0.5 * (C->value + D->value);

MDA.x = 0.5 * (D->x + A->x);
MDA.y = 0.5 * (D->y + A->y);
MDA.value = 0.5 * (D->value + A->value);

/*
 * Centre point.
 */
Ctr.x = 0.5 * (MBC.x + MDA.x);
Ctr.y = 0.5 * (MBC.y + MDA.y);
Ctr.value = 0.5 * (MBC.value + MDA.value);

/*
 * Now, contour each of the 8 triangles.
 */
if ( !(A->value > level  && MAB.value > level && MDA.value > level) &&
   !(A->value < level  && MAB.value < level && MDA.value < level)    )
   ssp_ContourTriangle (A, &MAB, &MDA, level);

if ( !(Ctr.value > level  && MAB.value > level && MDA.value > level) &&
   !(Ctr.value < level  && MAB.value < level && MDA.value < level)    )
   ssp_ContourTriangle (&MAB, &Ctr, &MDA, level);

if ( !(B->value > level  && MAB.value > level && MBC.value > level) &&
   !(B->value < level  && MAB.value < level && MBC.value < level)    )
   ssp_ContourTriangle (&MAB, B, &MBC, level);

if ( !(Ctr.value > level  && MAB.value > level && MBC.value > level) &&
   !(Ctr.value < level  && MAB.value < level && MBC.value < level)    )
   ssp_ContourTriangle (&MAB, &MBC, &Ctr, level);

if ( !(C->value > level  && MCD.value > level && MBC.value > level) &&
   !(C->value < level  && MCD.value < level && MBC.value < level)    )
   ssp_ContourTriangle (&MBC, C, &MCD, level);

if ( !(Ctr.value > level  && MCD.value > level && MBC.value > level) &&
   !(Ctr.value < level  && MCD.value < level && MBC.value < level)    )
   ssp_ContourTriangle (&MBC, &MCD, &Ctr, level);

if ( !(Ctr.value > level  && MCD.value > level && MDA.value > level) &&
   !(Ctr.value < level  && MCD.value < level && MDA.value < level)    )
   ssp_ContourTriangle (&MCD, &MDA, &Ctr, level);

if ( !(D->value > level  && MCD.value > level && MDA.value > level) &&
   !(D->value < level  && MCD.value < level && MDA.value < level)    )
   ssp_ContourTriangle (&MCD, D, &MDA, level);

return (0);
}  /* end of ssp_ContourQuad () */

/*--------------------------------------------------------------*/

#if (PROTO)

int ssp_ContourQuad3D (struct point_3D *A,
		     struct point_3D *B,
		     struct point_3D *C,
		     struct point_3D *D,
		     double A_value,
		     double B_value,
		     double C_value,
		     double D_value,
		     double level)

#else

int ssp_ContourQuad3D (A, B, C, D,
                     A_value, B_value, C_value, D_value,
                     level)
struct point_3D *A, *B, *C, *D;
double A_value, B_value, C_value, D_value;
double level;

#endif

/*
 * Purpose...
 * -------
 * Given the 4 corners of a quadrilateral, split it into
 * 8 triangles and contour each of the triangles.
 * This arrangement follows the scheme described by
 * James Quirk in his PhD thesis.
 *
 * Input...
 * -----
 * A, B, C, D : pointers to the corner points
 *              These must be specified in a counter-clockwise
 *              order.
 *              C--------B
 *              |        |
 *              |        |
 *              D--------A
 *
 * A_value, B_value,
 * C_value, D_value : data values at the corner points
 *
 * level     : value of the function at which the contour is to
 *             be drawn.
 *
 * Intermediate Points ...
 * -------------------
 *
 *            C --- MBC --- B
 *            |  \   |  /   |
 *           MCD -- Ctr -- MAB
 *            |  /   |  \   |
 *            D --- MDA --- A
 *
 * Revisions...
 * ---------
 * 08-Jun-93 : adapted from ssp_ContourQuad()
 *
 */

{  /* begin ssp_ContourQuad3D () */
struct point_3D MAB, MBC, MCD, MDA, Ctr;
double MAB_value, MBC_value, MCD_value, MDA_value, Ctr_value;

/*
 * Set up midpoints along each side
 */
MAB.x = 0.5 * (A->x + B->x);
MAB.y = 0.5 * (A->y + B->y);
MAB.z = 0.5 * (A->z + B->z);
MAB_value = 0.5 * (A_value + B_value);

MBC.x = 0.5 * (B->x + C->x);
MBC.y = 0.5 * (B->y + C->y);
MBC.z = 0.5 * (B->z + C->z);
MBC_value = 0.5 * (B_value + C_value);

MCD.x = 0.5 * (C->x + D->x);
MCD.y = 0.5 * (C->y + D->y);
MCD.z = 0.5 * (C->z + D->z);
MCD_value = 0.5 * (C_value + D_value);

MDA.x = 0.5 * (D->x + A->x);
MDA.y = 0.5 * (D->y + A->y);
MDA.z = 0.5 * (D->z + A->z);
MDA_value = 0.5 * (D_value + A_value);

/*
 * Centre point.
 */
Ctr.x = 0.5 * (MBC.x + MDA.x);
Ctr.y = 0.5 * (MBC.y + MDA.y);
Ctr.z = 0.5 * (MBC.z + MDA.z);
Ctr_value = 0.5 * (MBC_value + MDA_value);

/*
 * Now, contour each of the 8 triangles.
 */
if ( !(A_value > level  && MAB_value > level && MDA_value > level) &&
   !(A_value < level  && MAB_value < level && MDA_value < level)    )
   ssp_ContourTriangle3D
      (A, &MAB, &MDA, A_value, MAB_value, MDA_value, level);

if ( !(Ctr_value > level  && MAB_value > level && MDA_value > level) &&
   !(Ctr_value < level  && MAB_value < level && MDA_value < level)    )
   ssp_ContourTriangle3D
      (&MAB, &Ctr, &MDA, MAB_value, Ctr_value, MDA_value, level);

if ( !(B_value > level  && MAB_value > level && MBC_value > level) &&
   !(B_value < level  && MAB_value < level && MBC_value < level)    )
   ssp_ContourTriangle3D
      (&MAB, B, &MBC, MAB_value, B_value, MBC_value, level);

if ( !(Ctr_value > level  && MAB_value > level && MBC_value > level) &&
   !(Ctr_value < level  && MAB_value < level && MBC_value < level)    )
   ssp_ContourTriangle3D
      (&MAB, &MBC, &Ctr, MAB_value, MBC_value, Ctr_value, level);

if ( !(C_value > level  && MCD_value > level && MBC_value > level) &&
   !(C_value < level  && MCD_value < level && MBC_value < level)    )
   ssp_ContourTriangle3D
      (&MBC, C, &MCD, MBC_value, C_value, MCD_value, level);

if ( !(Ctr_value > level  && MCD_value > level && MBC_value > level) &&
   !(Ctr_value < level  && MCD_value < level && MBC_value < level)    )
   ssp_ContourTriangle3D
      (&MBC, &MCD, &Ctr, MBC_value, MCD_value, Ctr_value, level);

if ( !(Ctr_value > level  && MCD_value > level && MDA_value > level) &&
   !(Ctr_value < level  && MCD_value < level && MDA_value < level)    )
   ssp_ContourTriangle3D
      (&MCD, &MDA, &Ctr, MCD_value, MDA_value, Ctr_value, level);

if ( !(D_value > level  && MCD_value > level && MDA_value > level) &&
   !(D_value < level  && MCD_value < level && MDA_value < level)    )
   ssp_ContourTriangle3D
      (&MCD, D, &MDA, MCD_value, D_value, MDA_value, level);

return (0);
}  /* end of ssp_ContourQuad3D () */

/*---------------------------------------------------------------*/

#if (PROTO)

int ssp_ContourArray (struct ssp_data_array *F,
		      double level)

#else

int ssp_ContourArray (F, level)
struct ssp_data_array *F;
double level;

#endif

/*
 * Purpose...
 * -------
 * Given the two-dimensional array of positions and values,
 * plot the level contour one quadrilateral at a time.
 *
 * Input...
 * -----
 * F   : pointer to the array data structure
 *       This structure contains the arrays of x, y positions and
 *       the function values as well as the minimum and maximum
 *       indices for array addressing.
 * level : contour value
 *
 */

{  /* begin ssp_ContourArray() */
int    ix, iy;
struct ssp_data_point A, B, C, D;

for (ix = F->ixmin; ix < F->ixmax; ++ix)
   for (iy = F->iymin; iy < F->iymax; ++iy)
      {
      /*
       * Set up the data for each of the four corners.
       */
      A.x = F->x[ix+1][iy];
      A.y = F->y[ix+1][iy];
      A.value = F->value[ix+1][iy];

      B.x = F->x[ix+1][iy+1];
      B.y = F->y[ix+1][iy+1];
      B.value = F->value[ix+1][iy+1];

      C.x = F->x[ix][iy+1];
      C.y = F->y[ix][iy+1];
      C.value = F->value[ix][iy+1];

      D.x = F->x[ix][iy];
      D.y = F->y[ix][iy];
      D.value = F->value[ix][iy];

      if ( !(A.value > level  &&  B.value > level &&
	     C.value > level  &&  D.value > level) &&
	   !(A.value < level  &&  B.value < level &&
	     C.value < level  &&  D.value < level)  )
	 ssp_ContourQuad (&A, &B, &C, &D, level);
      }

return (0);
}  /* end of ssp_ContourArray() */

/*---------------------------------------------------------------*/

int ssp_PlotMesh (struct ssp_data_array *F, int ixskip, int iyskip)

/*
 * Purpose...
 * -------
 * Given the two-dimensional array of positions and values,
 * plot the (x, y) mesh.
 *
 * Input...
 * -----
 * F   : pointer to the array data structure
 *       This structure contains the arrays of x, y positions and
 *       the function values as well as the minimum and maximum
 *       indices for array addressing.
 * ixskip : increment for the ix direction (usually 1)
 * iyskip : increment for the iy direction (usually 1)
 *
 */

{  /* begin ssp_PlotMesh() */
int    ix, iy;
double x, y;

/*
 * Plot the vertical lines.
 */
for (ix = F->ixmin; ix <= F->ixmax; ix += ixskip)
   {
   iy = F->iymin;
   x = F->x[ix][iy];
   y = F->y[ix][iy];
   ssp_Move (x, y);

   for (iy = F->iymin + 1; iy <= F->iymax; ++iy)
      {
      x = F->x[ix][iy];
      y = F->y[ix][iy];
      ssp_Plot (x, y);
      }
   }

/*
 * Plot the horizontal lines.
 */
for (iy = F->iymin; iy <= F->iymax; iy += iyskip)
   {
   ix = F->ixmin;
   x = F->x[ix][iy];
   y = F->y[ix][iy];
   ssp_Move (x, y);

   for (ix = F->ixmin + 1; ix <= F->ixmax; ++ix)
      {
      x = F->x[ix][iy];
      y = F->y[ix][iy];
      ssp_Plot (x, y);
      }
   }

return (0);
}  /* end of ssp_PlotMesh() */

/*-----------------------------------------------------------------*/

#if (PROTO)

int ssp_PlotEdge (struct ssp_data_array *F)

#else

int ssp_PlotEdge (F)
struct ssp_data_array *F;

#endif

/*
 * Purpose...
 * -------
 * Given the two-dimensional array of positions and values,
 * plot the edge of the (x, y) mesh.
 *
 * Input...
 * -----
 * F   : pointer to the array data structure
 *       This structure contains the arrays of x, y positions and
 *       the function values as well as the minimum and maximum
 *       indices for array addressing.
 *
 */

{  /* begin ssp_PlotEdge() */
int    ix, iy;
double x, y;

/*
 * "West" boundary.
 */
ix = F->ixmin;
iy = F->iymin;
x = F->x[ix][iy];
y = F->y[ix][iy];
ssp_Move (x, y);

for (iy = F->iymin + 1; iy <= F->iymax; ++iy)
   {
   x = F->x[ix][iy];
   y = F->y[ix][iy];
   ssp_Plot (x, y);
   }

/*
 * "East" boundary.
 */
ix = F->ixmax;
iy = F->iymin;
x = F->x[ix][iy];
y = F->y[ix][iy];
ssp_Move (x, y);

for (iy = F->iymin + 1; iy <= F->iymax; ++iy)
   {
   x = F->x[ix][iy];
   y = F->y[ix][iy];
   ssp_Plot (x, y);
   }

/*
 * "South" boundary.
 */
iy = F->iymin;
ix = F->ixmin;
x = F->x[ix][iy];
y = F->y[ix][iy];
ssp_Move (x, y);

for (ix = F->ixmin + 1; ix <= F->ixmax; ++ix)
   {
   x = F->x[ix][iy];
   y = F->y[ix][iy];
   ssp_Plot (x, y);
   }

/*
 * "North" boundary
 */
iy = F->iymax;
ix = F->ixmin;
x = F->x[ix][iy];
y = F->y[ix][iy];
ssp_Move (x, y);

for (ix = F->ixmin + 1; ix <= F->ixmax; ++ix)
   {
   x = F->x[ix][iy];
   y = F->y[ix][iy];
   ssp_Plot (x, y);
   }

return (0);
}  /* end of ssp_PlotEdge() */

/*-------------------------------------------------------------*/

#if (PROTO)

int ssp_ColourFillTriangle (struct ssp_data_point *A,
		            struct ssp_data_point *B,
  		            struct ssp_data_point *C,
		            double level_min,
		            double level_max,
		            int grey_fill )

#else

int ssp_ColourFillTriangle (A, B, C, level_min, level_max, grey_fill)
struct ssp_data_point *A, *B, *C;
double level_max, level_min;
int    grey_fill;

#endif

/*
 * Purpose...
 * -------
 * Given the three vertices of a triangle, assume a linear
 * distribution and colour-code or grey-shade the levels.
 * <<< Two-Dimensional Version. >>>
 *
 * Input...
 * -----
 * A, B, C  : pointers to the corner vertices
 *            (Assume that they are in a counter-clockwise order.)
 *
 *                     B -------- A
 *                      \       /
 *                        \   /
 *                          C
 *
 * level_max  : maximum level expected (mapped to red)
 * level_min  : minimum level expected (mapped to blue)
 * grey_fill  : == 1: use grey-shading
 *              == 0: use colour filling
 *
 */

{  /* begin ssp_ContourTriangle () */
double xAB, yAB, xBC, yBC, xCA, yCA;
double mid_value;
double xvtx[3], yvtx[3];
double hue, bright, sat;

/*
 * Split the triangle into four pieces using the mid-point of
 * each side as a new vertex.
 */
mid_value = (A->value + B->value + C->value) / 3.0;

xAB = 0.5 * (A->x + B->x);
yAB = 0.5 * (A->y + B->y);

xCA = 0.5 * (A->x + C->x);
yCA = 0.5 * (A->y + C->y);

xBC = 0.5 * (B->x + C->x);
yBC = 0.5 * (B->y + C->y);

/*
 * Now colour each of the sub-triangles.
 */

/* Corner nearest A */
xvtx[0] = A->x; yvtx[0] = A->y;
xvtx[1] = xAB;  yvtx[1] = yAB;
xvtx[2] = xCA;  yvtx[2] = yCA;
if ( grey_fill == 1 )
   {
   bright = (A->value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon (3, xvtx, yvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (A->value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon (3, xvtx, yvtx, 0, 1, 0);
   }

/* Corner nearest B */
xvtx[0] = B->x; yvtx[0] = B->y;
xvtx[1] = xBC;  yvtx[1] = yBC;
xvtx[2] = xAB;  yvtx[2] = yAB;
if ( grey_fill == 1 )
   {
   bright = (B->value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon (3, xvtx, yvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (B->value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon (3, xvtx, yvtx, 0, 1, 0);
   }

/* Corner nearest C */
xvtx[0] = C->x; yvtx[0] = C->y;
xvtx[1] = xCA;  yvtx[1] = yCA;
xvtx[2] = xBC;  yvtx[2] = yBC;
if ( grey_fill == 1 )
   {
   bright = (C->value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon (3, xvtx, yvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (C->value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon (3, xvtx, yvtx, 0, 1, 0);
   }

/* Middle sub-triangle */
xvtx[0] = xBC;  yvtx[0] = yBC;
xvtx[1] = xCA;  yvtx[1] = yCA;
xvtx[2] = xAB;  yvtx[2] = yAB;
if ( grey_fill == 1 )
   {
   bright = (mid_value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon (3, xvtx, yvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (mid_value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon (3, xvtx, yvtx, 0, 1, 0);
   }

return (0);
}  /* end of ssp_ColourFillTriangle () */


/*--------------------------------------------------------------*/

#if (PROTO)

int ssp_ColourFillQuad (struct ssp_data_point *A,
		        struct ssp_data_point *B,
		        struct ssp_data_point *C,
  		        struct ssp_data_point *D,
		        double level_min,
		        double level_max,
		        int    grey_fill )

#else

int ssp_ColourFillQuad (A, B, C, D, level_min, level_max, grey_fill )
struct ssp_data_point *A, *B, *C, *D;
double level_min, level_max;
int    grey_fill;

#endif

/*
 * Purpose...
 * -------
 * Given the 4 corners of a quadrilateral, split it into
 * 8 triangles and colour-fill each of the triangles.
 * This arrangement follows the scheme described by
 * James Quirk in his PhD thesis.
 *
 * Input...
 * -----
 * A, B, C, D : pointers to the corner data structures
 *              These must be specified in a counter-clockwise
 *              order.
 *              C--------B
 *              |        |
 *              |        |
 *              D--------A
 *
 * level_min : lowest expected value (mapped to blue)
 * level_max : highest expected value (mapped to red)
 * grey_fill  : == 1: use grey-shading
 *              == 0: use colour filling
 *
 * Intermediate Points ...
 * -------------------
 *
 *            C --- MBC --- B
 *            |  \   |  /   |
 *           MCD -- Ctr -- MAB
 *            |  /   |  \   |
 *            D --- MDA --- A
 *
 * Revisions...
 * ---------
 * 1.0 : 01-Jan-94 : Adapted from ssp_ContourQuad().
 * 1.1 : 12-Jun-94 : grey-filling added
 *
 */

{  /* begin ssp_ColourFillQuad () */
struct ssp_data_point MAB, MBC, MCD, MDA, Ctr;

/*
 * Set up midpoints along each side
 */
MAB.x = 0.5 * (A->x + B->x);
MAB.y = 0.5 * (A->y + B->y);
MAB.value = 0.5 * (A->value + B->value);

MBC.x = 0.5 * (B->x + C->x);
MBC.y = 0.5 * (B->y + C->y);
MBC.value = 0.5 * (B->value + C->value);

MCD.x = 0.5 * (C->x + D->x);
MCD.y = 0.5 * (C->y + D->y);
MCD.value = 0.5 * (C->value + D->value);

MDA.x = 0.5 * (D->x + A->x);
MDA.y = 0.5 * (D->y + A->y);
MDA.value = 0.5 * (D->value + A->value);

/*
 * Centre point.
 */
Ctr.x = 0.5 * (MBC.x + MDA.x);
Ctr.y = 0.5 * (MBC.y + MDA.y);
Ctr.value = 0.5 * (MBC.value + MDA.value);

/*
 * Now, colour each of the 8 triangles.
 */
ssp_ColourFillTriangle (A, &MAB, &MDA, level_min, level_max, grey_fill);

ssp_ColourFillTriangle (&MAB, &Ctr, &MDA, level_min, level_max, grey_fill);

ssp_ColourFillTriangle (&MAB, B, &MBC, level_min, level_max, grey_fill);

ssp_ColourFillTriangle (&MAB, &MBC, &Ctr, level_min, level_max, grey_fill);

ssp_ColourFillTriangle (&MBC, C, &MCD, level_min, level_max, grey_fill);

ssp_ColourFillTriangle (&MBC, &MCD, &Ctr, level_min, level_max, grey_fill);

ssp_ColourFillTriangle (&MCD, &MDA, &Ctr, level_min, level_max, grey_fill);

ssp_ColourFillTriangle (&MCD, D, &MDA, level_min, level_max, grey_fill);

return (0);
}  /* end of ssp_ColourFillQuad () */


/*---------------------------------------------------------------*/

#if (PROTO)

int ssp_ColourFillArray (struct ssp_data_array *F,
		         double level_min, double level_max,
		         int grey_fill )

#else

int ssp_ColourFillArray (F, level_min, level_max, grey_fill)
struct ssp_data_array *F;
double level_min, level_max;
int    grey_fill;

#endif

/*
 * Purpose...
 * -------
 * Given the two-dimensional array of positions and values,
 * colour the quadrilaterals, one at a time.
 *
 * Input...
 * -----
 * F   : pointer to the array data structure
 *       This structure contains the arrays of x, y positions and
 *       the function values as well as the minimum and maximum
 *       indices for array addressing.
 * level_min : lowest expected value (mapped to blue)
 * level_max : highest expected value (mapped to red)
 * grey_fill  : == 1: use grey-shading
 *              == 0: use colour filling
 *
 */

{  /* begin ssp_ColourFillArray() */
int    ix, iy;
struct ssp_data_point A, B, C, D;

for (ix = F->ixmin; ix < F->ixmax; ++ix)
   for (iy = F->iymin; iy < F->iymax; ++iy)
      {
      /*
       * Set up the data for each of the four corners.
       */
      A.x = F->x[ix+1][iy];
      A.y = F->y[ix+1][iy];
      A.value = F->value[ix+1][iy];

      B.x = F->x[ix+1][iy+1];
      B.y = F->y[ix+1][iy+1];
      B.value = F->value[ix+1][iy+1];

      C.x = F->x[ix][iy+1];
      C.y = F->y[ix][iy+1];
      C.value = F->value[ix][iy+1];

      D.x = F->x[ix][iy];
      D.y = F->y[ix][iy];
      D.value = F->value[ix][iy];

      ssp_ColourFillQuad (&A, &B, &C, &D, level_min, level_max, grey_fill);
      }

return (0);
}  /* end of ssp_ColourFillArray() */

/*-------------------------------------------------------------*/

#if (PROTO)

int ssp_ColourFillTriangle3D (struct point_3D *A,
		              struct point_3D *B,
  		              struct point_3D *C,
		              double A_value,
		              double B_value,
		              double C_value,
		              double level_min,
		              double level_max,
		              int    grey_fill )

#else

int ssp_ColourFillTriangle3D (A, B, C, A_value, B_value, C_value,
                              level_min, level_max, grey_fill )
struct point_3D *A, *B, *C;
double A_value, B_value, C_value;
double level_max, level_min;
int    grey_fill;

#endif

/*
 * Purpose...
 * -------
 * Given the three vertices of a triangle, assume a linear
 * distribution and colour-code the levels.
 * <<< Three-Dimensional Version. >>>
 *
 * Input...
 * -----
 * A, B, C  : pointers to the corner vertices
 *            (Assume that they are in a counter-clockwise order.)
 * A_value,
 * B_value,
 * C_value  : Data values at the vertices
 *
 *                     B -------- A
 *                      \       /
 *                        \   /
 *                          C
 *
 * level_max  : maximum level expected (mapped to red)
 * level_min  : minimum level expected (mapped to blue)
 * grey_fill  : == 1: use grey-shading
 *              == 0: use colour filling
 *
 */

{  /* begin ssp_ContourTriangle3D () */
double xAB, yAB, zAB, xBC, yBC, zBC, xCA, yCA, zCA;
double mid_value;
struct point_3D pvtx[3];
double hue, bright, sat;

/*
 * Split the triangle into four pieces using the mid-point of
 * each side as a new vertex.
 */
mid_value = (A_value + B_value + C_value) / 3.0;

xAB = 0.5 * (A->x + B->x);
yAB = 0.5 * (A->y + B->y);
zAB = 0.5 * (A->z + B->z);

xCA = 0.5 * (A->x + C->x);
yCA = 0.5 * (A->y + C->y);
zCA = 0.5 * (A->z + C->z);

xBC = 0.5 * (B->x + C->x);
yBC = 0.5 * (B->y + C->y);
zBC = 0.5 * (B->z + C->z);

/*
 * Now colour each of the sub-triangles.
 */

/* Corner nearest A */
pvtx[0].x = A->x; pvtx[0].y = A->y; pvtx[0].z = A->z;
pvtx[1].x = xAB;  pvtx[1].y = yAB;  pvtx[1].z = zAB;
pvtx[2].x = xCA;  pvtx[2].y = yCA;  pvtx[2].z = zCA;
if ( grey_fill == 1 )
   {
   bright = (A_value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon3D (3, pvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (A_value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon3D (3, pvtx, 0, 1, 0);
   }

/* Corner nearest B */
pvtx[0].x = B->x; pvtx[0].y = B->y; pvtx[0].z = B->z;
pvtx[1].x = xBC;  pvtx[1].y = yBC;  pvtx[1].z = zBC;
pvtx[2].x = xAB;  pvtx[2].y = yAB;  pvtx[2].z = zAB;
if ( grey_fill == 1 )
   {
   bright = (B_value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon3D (3, pvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (B_value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon3D (3, pvtx, 0, 1, 0);
   }

/* Corner nearest C */
pvtx[0].x = C->x; pvtx[0].y = C->y; pvtx[0].z = C->z;
pvtx[1].x = xCA;  pvtx[1].y = yCA;  pvtx[1].z = zCA;
pvtx[2].x = xBC;  pvtx[2].y = yBC;  pvtx[2].z = zBC;
if ( grey_fill == 1 )
   {
   bright = (C_value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon3D (3, pvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (C_value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon3D (3, pvtx, 0, 1, 0);
   }

/* Middle sub-triangle */
pvtx[0].x = xBC;  pvtx[0].y = yBC;  pvtx[0].z = zBC;
pvtx[1].x = xCA;  pvtx[1].y = yCA;  pvtx[1].z = zCA;
pvtx[2].x = xAB;  pvtx[2].y = yAB;  pvtx[2].z = zAB;
if ( grey_fill == 1 )
   {
   bright = (mid_value - level_min) / (level_max - level_min);
   ssp_SetFillGrey (bright);
   ssp_Polygon3D (3, pvtx, 1, 0, 0);
   }
else
   {
   hue = ssp_MapToColour (mid_value, level_min, level_max);
   bright = 1.0;
   sat = 1.0;
   ssp_SetHSBFillColour (hue, sat, bright);
   ssp_Polygon3D (3, pvtx, 0, 1, 0);
   }

return (0);
}  /* end of ssp_ColourFillTriangle3D () */


/*--------------------------------------------------------------*/

#if (PROTO)

int ssp_ColourFillQuad3D (struct point_3D *A,
		          struct point_3D *B,
		          struct point_3D *C,
		          struct point_3D *D,
		          double A_value,
		          double B_value,
		          double C_value,
		          double D_value,
		          double level_min,
		          double level_max,
		          int    grey_fill )

#else

int ssp_ColourFillQuad3D (A, B, C, D,
                          A_value, B_value, C_value, D_value,
                          level_min, level_max,
                          grey_fill )
struct point_3D *A, *B, *C, *D;
double A_value, B_value, C_value, D_value;
double level_min, level_max;
int    grey_fill;

#endif

/*
 * Purpose...
 * -------
 * Given the 4 corners of a quadrilateral, split it into
 * 8 triangles and contour each of the triangles.
 * This arrangement follows the scheme described by
 * James Quirk in his PhD thesis.
 *
 * Input...
 * -----
 * A, B, C, D : pointers to the corner points
 *              These must be specified in a counter-clockwise
 *              order.
 *              C--------B
 *              |        |
 *              |        |
 *              D--------A
 *
 * A_value, B_value,
 * C_value, D_value : data values at the corner points
 *
 * level     : value of the function at which the contour is to
 *             be drawn.
 * level_max  : maximum level expected (mapped to red)
 * level_min  : minimum level expected (mapped to blue)
 * grey_fill  : == 1: use grey-shading
 *              == 0: use colour filling
 *
 * Intermediate Points ...
 * -------------------
 *
 *            C --- MBC --- B
 *            |  \   |  /   |
 *           MCD -- Ctr -- MAB
 *            |  /   |  \   |
 *            D --- MDA --- A
 *
 * Revisions...
 * ---------
 * 01-Jan-94 : adapted from ssp_ContourQuad3D()
 * 12-Jun-94 : grey-filling included
 *
 */

{  /* begin ssp_ColourFillQuad3D () */
struct point_3D MAB, MBC, MCD, MDA, Ctr;
double MAB_value, MBC_value, MCD_value, MDA_value, Ctr_value;

/*
 * Set up midpoints along each side
 */
MAB.x = 0.5 * (A->x + B->x);
MAB.y = 0.5 * (A->y + B->y);
MAB.z = 0.5 * (A->z + B->z);
MAB_value = 0.5 * (A_value + B_value);

MBC.x = 0.5 * (B->x + C->x);
MBC.y = 0.5 * (B->y + C->y);
MBC.z = 0.5 * (B->z + C->z);
MBC_value = 0.5 * (B_value + C_value);

MCD.x = 0.5 * (C->x + D->x);
MCD.y = 0.5 * (C->y + D->y);
MCD.z = 0.5 * (C->z + D->z);
MCD_value = 0.5 * (C_value + D_value);

MDA.x = 0.5 * (D->x + A->x);
MDA.y = 0.5 * (D->y + A->y);
MDA.z = 0.5 * (D->z + A->z);
MDA_value = 0.5 * (D_value + A_value);

/*
 * Centre point.
 */
Ctr.x = 0.5 * (MBC.x + MDA.x);
Ctr.y = 0.5 * (MBC.y + MDA.y);
Ctr.z = 0.5 * (MBC.z + MDA.z);
Ctr_value = 0.5 * (MBC_value + MDA_value);

/*
 * Now, colour-fill each of the 8 triangles.
 */
ssp_ColourFillTriangle3D
   (A, &MAB, &MDA, A_value, MAB_value, MDA_value,
    level_min, level_max, grey_fill);

ssp_ColourFillTriangle3D
   (&MAB, &Ctr, &MDA, MAB_value, Ctr_value, MDA_value,
    level_min, level_max, grey_fill);

ssp_ColourFillTriangle3D
   (&MAB, B, &MBC, MAB_value, B_value, MBC_value,
    level_min, level_max, grey_fill);

ssp_ColourFillTriangle3D
   (&MAB, &MBC, &Ctr, MAB_value, MBC_value, Ctr_value,
    level_min, level_max, grey_fill);

ssp_ColourFillTriangle3D
   (&MBC, C, &MCD, MBC_value, C_value, MCD_value,
    level_min, level_max, grey_fill);

ssp_ColourFillTriangle3D
   (&MBC, &MCD, &Ctr, MBC_value, MCD_value, Ctr_value,
    level_min, level_max, grey_fill);

ssp_ColourFillTriangle3D
   (&MCD, &MDA, &Ctr, MCD_value, MDA_value, Ctr_value,
    level_min, level_max, grey_fill);

ssp_ColourFillTriangle3D
   (&MCD, D, &MDA, MCD_value, D_value, MDA_value,
    level_min, level_max, grey_fill);

return (0);
}  /* end of ssp_ColourFillQuad3D () */

/*-----------------------------------------------------------------*/
