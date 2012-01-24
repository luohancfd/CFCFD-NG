/** \file moc_bezier.c
 * \ingroup imoc
 * \brief  Functions for handling Bezier curves of third order.
 * 
 * These routines provide support for Bezier curves in three
 * dimensions.  Common uses are CAD and font outlines.
 * The IMOC user should not need to call these procedures directly.
 *
 * \author PA Jacobs
 * 38 Ellington St
 * Ekibin, Qld 4121
 * AUSTRALIA
 *
 * \verbatim
 * Version...
 * -------
 * 1.0  : basic routines for third-order curves
 * 1.1  : 06-nov-92 : full function prototypes
 * 1.2  : 15-nov-92 : trap underflow in eval_bezier_3_poly()
 * 1.3  : 18-nov-92 : fix the polyline evaluation routine so that the
 *                    t_star array is never addressed incorrectly
 *                    (when we have only one segment in the polyline)
 * 1.4  : 22-Jan-93 : add slope calculation deriv_bezier_3()
 * 1.5  : 25-May-93 : Changed location_3D to point_3D.
 * 2.0  : 21-Feb-95 : Added N-degree Bezier curves
 * 2.1  : 11-Mar-95 : Compute the length of the Bezier segments by evaluating
 *                    points along the segment (rather than summing the control
 *                    point segments).
 * 2.2  : 05-Jul-98 : Cubic spline added.
 * 3.0  : 04-Jan-2000 : adapted for use in I-MOC program.
 *                      Added n_max to the Bezier polyline and allowed the
 *                      number of segments to be less than that originally
 *                      allocated.
 *
 * References...
 * ----------
 * Rogers and Adams (1990)
 * Mathematical Elements for Computer Graphics. 2nd ed
 * McGraw Hill
 *
 * G. Farin (1990)
 * Curves and surfaces for computer aided Geometric design. 2nd ed
 * Academic Press
 *
 * Fujio Yamaguchi (1988)
 * Curves and surfaces in computer aided geometric design.
 * Springer-Verlag
 * \endverbatim
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "moc_bezier.h"

/*--------------------------------------------------------------*/

int alloc_bezier_3_poly (struct bezier_3_poly_data *bp, int n)
/* Purpose...
 * -------
 * Allocate memory for the Bezier polyline data.
 *
 * Input...
 * -----
 * bp	: pointer to the polyline data structure
 * n	: number of Bezier segments in the polyline
 *
 * Output...
 * ------
 * The function returns -1 if the allocation fails, 0 otherwise.
 *
 */

{   /* begin alloc_bezier_3_poly() */
int flag, i;

flag = 0;

if (n <= 0)
   {
   flag = -1;
   return (flag);
   }

bp->b3_seg = (struct bezier_3_data *)
	     malloc (n * sizeof(struct bezier_3_data) );
if (bp->b3_seg == NULL)
   {
   flag = -1;
   return (flag);
   }

bp->t_star = (double *) malloc ( n * sizeof(double) );
if (bp->t_star == NULL)
   {
   flag = -1;
   return (flag);
   }

for (i = 0; i < n; ++i)
   {
   /* Set all of the orders to 3 */
   bp->b3_seg[i].n = 3;
   /* Set the normalized distance along the polyline to
    * nominal values */
   bp->t_star[i] = (double) (i + 1.0) / n;
   }

/* This is the assumed number of Bezier elements and
 * the maximum allowed in the allocated arrays. 
 */
bp->n     = n;
bp->n_max = n;

return (flag);
}


/*---------------------------------------------------------------*/

int init_bezier_3 (struct bezier_3_data *b3,
		   struct point_3D *loc0,
		   struct point_3D *loc1,
		   struct point_3D *loc2,
		   struct point_3D *loc3)
/* Purpose...
 * -------
 * Initialize the coordinates of the Bezier control points.
 *
 * Input...
 * -----
 * b3	: (pointer to) the Bezier curve data
 * loc0 : (pointer to) the location of control point 0
 * loc1
 * loc2
 * loc3
 *
 */

{
b3->B[0].x = loc0->x;
b3->B[0].y = loc0->y;
b3->B[0].z = loc0->z;

b3->B[1].x = loc1->x;
b3->B[1].y = loc1->y;
b3->B[1].z = loc1->z;

b3->B[2].x = loc2->x;
b3->B[2].y = loc2->y;
b3->B[2].z = loc2->z;

b3->B[3].x = loc3->x;
b3->B[3].y = loc3->y;
b3->B[3].z = loc3->z;

b3->n = 3;

return (0);
}


/*---------------------------------------------------------------*/

int init_bezier_n ( int n,
                    struct bezier_n_data *bn,
		    struct point_3D locarray[] )
/* Purpose...
 * -------
 * Initialize the coordinates of the Bezier control points.
 *
 * Input...
 * -----
 * n        : degree of the Bezier curve
 *            (must be less than MAX_BEX_DEGREE else it will truncate)
 *            Also, note that there are n+1 polygon points for an
 *            n-degree Bezier curve.
 * bn       : (pointer to) the Bezier curve data
 * locarray : (pointer to) the array of control points
 *
 */

{
int i;

if (n > MAX_BEZ_DEGREE) n = MAX_BEZ_DEGREE;

bn->n = n;

for (i = 0; i <= n; ++i)
   {
   bn->B[i].x = locarray[i].x;
   bn->B[i].y = locarray[i].y;
   bn->B[i].z = locarray[i].z;
   }

return (0);
}


/*---------------------------------------------------------------*/

int dump_bezier_3( struct bezier_3_data *b3 ) {
   /* Write out all of the data. */
   printf( "B[0] %g %g %g B[1] %g %g %g\n",
      b3->B[0].x, b3->B[0].y, b3->B[0].z,
      b3->B[1].x, b3->B[1].y, b3->B[1].z );
   printf( "B[2] %g %g %g B[3] %g %g %g\n",
      b3->B[2].x, b3->B[2].y, b3->B[2].z,
      b3->B[3].x, b3->B[3].y, b3->B[3].z );
   return 0;
} /* end function dump_bezier_3 */


/*---------------------------------------------------------------*/

int eval_bezier_3 (struct bezier_3_data *b3,
		   double t,
		   struct point_3D *loc)
/* Purpose...
 * -------
 * Evaluate the third-order Bezier curve position with
 * parameter t.
 *
 * Input...
 * -----
 * b3	: (pointer to) the Bezier curve data
 * t	: position parameter, 0 <= t <= 1.0
 *
 * Output...
 * ------
 * loc	: (pointer to) the computed location in 3D
 *
 * The function returns 0 if everything was OK, -1 if there was
 * a problem.
 *
 */

{  /* begin eval_bezier_3() */
double one_m_t, J30, J31, J32, J33;

if (t < -0.01 || t > 1.01)
   {
   return (-1);
   }

one_m_t = 1.0 - t;

J30 = one_m_t * one_m_t * one_m_t;
J31 = 3.0 * t * one_m_t * one_m_t;
J32 = 3.0 * t * t * one_m_t;
J33 = t * t * t;

loc->x = b3->B[0].x * J30 + b3->B[1].x * J31 +
	 b3->B[2].x * J32 + b3->B[3].x * J33;
loc->y = b3->B[0].y * J30 + b3->B[1].y * J31 +
	 b3->B[2].y * J32 + b3->B[3].y * J33;
loc->z = b3->B[0].z * J30 + b3->B[1].z * J31 +
	 b3->B[2].z * J32 + b3->B[3].z * J33;

return (0);
}


/*---------------------------------------------------------------*/

int eval_bezier_n (struct bezier_n_data *bn,
		   double t,
		   struct point_3D *loc)
/* Purpose...
 * -------
 * Evaluate the n-degree Bezier curve position with
 * parameter t.
 * The de Casteljau algorithm is used for simplicity not efficiency.
 *
 * Input...
 * -----
 * bn	: (pointer to) the Bezier curve data
 * t	: position parameter, 0 <= t <= 1.0
 *
 * Output...
 * ------
 * loc	: (pointer to) the computed location in 3D
 *
 * The function returns 0 if everything was OK, -1 if there was
 * a problem.
 *
 */

{  /* begin eval_bezier_n() */
double one_m_t;
struct point_3D Q[MAX_BEZ_DEGREE+1];
int    n, k, i;


if (t < -0.01 || t > 1.01)
   {
   return (-1);
   }

n = bn->n;              /* order of the curve */
one_m_t = 1.0 - t;

/*
 * Copy the control points into the work array.
 */
for (i = 0; i <= n; ++i)
   {
   Q[i].x = bn->B[i].x;
   Q[i].y = bn->B[i].y;
   Q[i].z = bn->B[i].z;
   }

/*
 * Now, generate one new level at a time;
 * over-writing the work array at each level.
 */
for (k = 1; k <= n; ++k)
   {
   for (i = 0; i <= n - k; ++i)
      {
      Q[i].x = one_m_t * Q[i].x + t * Q[i+1].x;
      Q[i].y = one_m_t * Q[i].y + t * Q[i+1].y;
      Q[i].z = one_m_t * Q[i].z + t * Q[i+1].z;
      }
   }  /* for (k... */

/*
 * The point on the curve is now in Q[0]
 */
loc->x = Q[0].x;
loc->y = Q[0].y;
loc->z = Q[0].z;

return (0);
}


/*---------------------------------------------------------------*/

int deriv_bezier_3 (struct bezier_3_data *b3,
		    double t,
		    struct point_3D *dxyzdt)
/* Purpose...
 * -------
 * Evaluate the derivative of the third-order Bezier curve
 * with respect to parameter t.
 *
 * Input...
 * -----
 * b3	: (pointer to) the Bezier curve data
 * t	: position parameter, 0 <= t <= 1.0
 *
 * Output...
 * ------
 * dxyzdt  : (pointer to) the computed derivative in 3D
 *
 * The function returns 0 if everything was OK, -1 if there was
 * a problem.
 *
 */

{  /* begin deriv_bezier_3() */
double one_m_t, dJ30dt, dJ31dt, dJ32dt, dJ33dt;

if (t < -0.01 || t > 1.01)
   {
   return (-1);
   }

one_m_t = 1.0 - t;

if (t < 1.0e-6) {
   dJ30dt = -3.0;
   dJ31dt =  3.0;
   dJ32dt =  0.0;
   dJ33dt =  0.0;
} else if (one_m_t < 1.0e-6) {
   dJ30dt =  0.0;
   dJ31dt =  0.0;
   dJ32dt = -3.0;
   dJ33dt =  3.0;
} else {
   dJ30dt = -3.0 * one_m_t * one_m_t;
   dJ31dt = 3.0 * one_m_t * one_m_t - 6.0 * t * one_m_t;
   dJ32dt = 6.0 * t * one_m_t - 3.0 * t * t;
   dJ33dt = 3.0 * t * t;
} /* end if */

dxyzdt->x = b3->B[0].x * dJ30dt + b3->B[1].x * dJ31dt +
	    b3->B[2].x * dJ32dt + b3->B[3].x * dJ33dt;
dxyzdt->y = b3->B[0].y * dJ30dt + b3->B[1].y * dJ31dt +
	    b3->B[2].y * dJ32dt + b3->B[3].y * dJ33dt;
dxyzdt->z = b3->B[0].z * dJ30dt + b3->B[1].z * dJ31dt +
	    b3->B[2].z * dJ32dt + b3->B[3].z * dJ33dt;

return (0);
}


/*---------------------------------------------------------------*/

int segment_bezier_3_poly (struct bezier_3_poly_data *bp,
			   struct point_3D *loc0,
			   struct point_3D *loc1,
			   struct point_3D *loc2,
			   struct point_3D *loc3,
			   int i)
/* Purpose...
 * -------
 * Initialize the coordinates of the Bezier control points for
 * Bezier segment i.
 *
 * Input...
 * -----
 * p3	: (pointer to) the Bezier polyline data
 * loc0 : (pointer to) the location of control point 0
 * loc1
 * loc2
 * loc3
 * i    : index of the desired segment, 0 <= i < n
 *
 */

{
if (i < 0 || i >= bp->n)
   {
   return (-1);
   }

bp->b3_seg[i].B[0].x = loc0->x;
bp->b3_seg[i].B[0].y = loc0->y;
bp->b3_seg[i].B[0].z = loc0->z;

bp->b3_seg[i].B[1].x = loc1->x;
bp->b3_seg[i].B[1].y = loc1->y;
bp->b3_seg[i].B[1].z = loc1->z;

bp->b3_seg[i].B[2].x = loc2->x;
bp->b3_seg[i].B[2].y = loc2->y;
bp->b3_seg[i].B[2].z = loc2->z;

bp->b3_seg[i].B[3].x = loc3->x;
bp->b3_seg[i].B[3].y = loc3->y;
bp->b3_seg[i].B[3].z = loc3->z;

bp->b3_seg[i].n = 3;

return (0);
}

/*---------------------------------------------------------------*/

int normalize_bezier_3_poly (struct bezier_3_poly_data *bp)
/* Purpose...
 * -------
 * Normalize the Bezier polyline so that the parameter t_star
 * varies from 0.0 to 1.0 along the entire curve.
 * Arc length along the Bezier segments is approximated as the
 * sum of the distances between a number of points along each
 * curve.
 *
 * Input...
 * -----
 * bp	: pointer to the polyline data structure
 *
 * Output...
 * ------
 * The function returns -1 if the normalization fails, 0 otherwise.
 *
 */

{
int    i;
double dx, dy, dz, arc_length;
double t, xa, ya, za, xb, yb, zb, seg_length;
struct point_3D loc, *p;
struct bezier_3_data *B3;

/*
 * Compute the approximate arc length along the polyline.
 */
arc_length = 0.0;
for (i = 0; i < bp->n; ++i)
   {
   B3 = &(bp->b3_seg[i]);
   seg_length = 0.0;

   p = &(B3->B[0]);  /* location of first control point on the segment */
   xa = p->x; ya = p->y; za = p->z;

   /* Step along the segment, summing the arc length as we go. */
   for ( t = 0.1; t <= 1.001; t += 0.1)
      {
      eval_bezier_3 ( B3, t, &loc );
      xb = loc.x; yb = loc.y; zb = loc.z;
      dx = xb - xa; dy = yb - ya; dz = zb - za;
      seg_length += sqrt( dx * dx + dy * dy + dz * dz );
      xa = xb; ya = yb; za = zb;  /* move the left point along */
      } /* endfor */

   /* Approximate arc length of this Bezier segment. */
   arc_length += seg_length;
   bp->t_star[i] = arc_length;
   }

if (arc_length <= 0.0)
   {
   return (-1);
   }

/*
 * Normalize the distance along the polyline so that
 * t_star varies from 0.0 to 1.0 along the entire polyline.
 * If there is only one segment, then t_star[0] will be 1.0
 * and no other segments will be valid.
 */
for (i = 0; i < bp->n; ++i) bp->t_star[i] /= arc_length;


return (0);
}


/*---------------------------------------------------------------*/

int dump_bezier_3_poly ( struct bezier_3_poly_data *bp ) {
   /* Write out all details. */
   int iseg;

   printf( "------- Begin Dump ---------\n");
   printf( "nseg = %d, t_star = ", bp->n );

   for (iseg = 0; iseg < bp->n; ++iseg) {
      printf( "%g ", bp->t_star[iseg] );
   } /* end for */
   printf( "\n" );
   
   for (iseg = 0; iseg < bp->n; ++iseg) {
      printf( "Segment %d\n", iseg );
      dump_bezier_3( &(bp->b3_seg[iseg]) );
   } /* end for */
   
   printf( "------- End Dump ---------\n");
   return 0;
} /* end function dump_bezier_3_poly */

/*---------------------------------------------------------------*/

int eval_bezier_3_poly (struct bezier_3_poly_data *bp,
			double t_star,
			struct point_3D *loc)
/* Purpose...
 * -------
 * Evaluate the position on the Bezier polyline at given
 * parameter t_star.
 *
 * Input...
 * -----
 * bp	  : (pointer to) the Bezier polyline data
 * t_star : position parameter, 0 <= t_star <= 1.0
 *
 * Output...
 * ------
 * loc	: (pointer to) the computed location in 3D
 *
 * The function returns 0 if everything was OK, -1 if there was
 * a problem.
 *
 */

{  /* begin eval_bezier_3_poly() */
double t, one_m_t, J30, J31, J32, J33;
int    i;

/*
 * This default behaviour may be unsuitable for some applications.
 */
if (t_star < 0.0) t_star = 0.0;
if (t_star > 1.0) t_star = 1.0;

/*
 * Find the relevant segment.
 */
i = 0;
while (t_star > bp->t_star[i] && i < (bp->n - 1) ) ++i;

/*
 * Compute the local parameter.
 */
if (i == 0)
   {
   t = t_star / bp->t_star[i];
   }
else
   {
   t = (t_star - bp->t_star[i-1]) /
       (bp->t_star[i] - bp->t_star[i-1]);
   }

/*
 * Now compute the position on the identified Bezier segment.
 */
one_m_t = 1.0 - t;

J30 = one_m_t * one_m_t * one_m_t;
J31 = 3.0 * t * one_m_t * one_m_t;
J32 = 3.0 * t * t * one_m_t;
J33 = t * t * t;

if (t < 1.0e-6)
   {
   loc->x = bp->b3_seg[i].B[0].x;
   loc->y = bp->b3_seg[i].B[0].y;
   loc->z = bp->b3_seg[i].B[0].z;
   }
else if (one_m_t < 1.0e-6)
   {
   loc->x = bp->b3_seg[i].B[3].x;
   loc->y = bp->b3_seg[i].B[3].y;
   loc->z = bp->b3_seg[i].B[3].z;
   }
else
   {
   loc->x = bp->b3_seg[i].B[0].x * J30 + bp->b3_seg[i].B[1].x * J31 +
	    bp->b3_seg[i].B[2].x * J32 + bp->b3_seg[i].B[3].x * J33;
   loc->y = bp->b3_seg[i].B[0].y * J30 + bp->b3_seg[i].B[1].y * J31 +
	    bp->b3_seg[i].B[2].y * J32 + bp->b3_seg[i].B[3].y * J33;
   loc->z = bp->b3_seg[i].B[0].z * J30 + bp->b3_seg[i].B[1].z * J31 +
	    bp->b3_seg[i].B[2].z * J32 + bp->b3_seg[i].B[3].z * J33;
   }

return (0);
} /* end function eval_bezier_3_poly */


/*---------------------------------------------------------------*/

int deriv_bezier_3_poly (struct bezier_3_poly_data *bp,
			double t_star,
			struct point_3D *dxyzdt) {
   /* Purpose...
    * -------
    * Evaluate the derivatives of the Bezier polyline at given
    * parameter t_star.
    *
    * Input...
    * -----
    * bp	  : (pointer to) the Bezier polyline data
    * t_star : position parameter, 0 <= t_star <= 1.0
    *
    * Output...
    * ------
    * dxyzdt : (pointer to) the computed derivatives in 3D
    *
    * The function returns 0 if everything was OK, -1 if there was
    * a problem.
    *
    */

   double t, dt_segment, one_m_t, dJ30dt, dJ31dt, dJ32dt, dJ33dt;
   int    i;

   /*
    * This default behaviour may be unsuitable for some applications.
    */
   if (t_star < 0.0) t_star = 0.0;
   if (t_star > 1.0) t_star = 1.0;

   /*
    * Find the relevant segment.
    */
   i = 0;
   while (t_star > bp->t_star[i] && i < (bp->n - 1) ) ++i;

   /*
    * Compute the local parameter.
    */
   if (i == 0) {
      dt_segment = bp->t_star[i];
      if (dt_segment <= 1.0e-12) return -1;
      t = t_star / dt_segment;
   } else {
      dt_segment = bp->t_star[i] - bp->t_star[i-1];
      if (dt_segment <= 1.0e-12) return -1;
      t = (t_star - bp->t_star[i-1]) / dt_segment;
   } /* end if */

   /*
    * Now compute the derivatives of the identified Bezier segment
    * with respect to the t parameter that ranges 0.0 .. 1.0 .
    */
   one_m_t = 1.0 - t;

   if (t < 1.0e-6) {
      dJ30dt = -3.0;
      dJ31dt =  3.0;
      dJ32dt =  0.0;
      dJ33dt =  0.0;
   } else if (one_m_t < 1.0e-6) {
      dJ30dt =  0.0;
      dJ31dt =  0.0;
      dJ32dt = -3.0;
      dJ33dt =  3.0;
   } else {
      dJ30dt = -3.0 * one_m_t * one_m_t;
      dJ31dt = 3.0 * one_m_t * one_m_t - 6.0 * t * one_m_t;
      dJ32dt = 6.0 * t * one_m_t - 3.0 * t * t;
      dJ33dt = 3.0 * t * t;
   }

   dxyzdt->x = bp->b3_seg[i].B[0].x * dJ30dt + bp->b3_seg[i].B[1].x * dJ31dt +
       bp->b3_seg[i].B[2].x * dJ32dt + bp->b3_seg[i].B[3].x * dJ33dt;
   dxyzdt->y = bp->b3_seg[i].B[0].y * dJ30dt + bp->b3_seg[i].B[1].y * dJ31dt +
       bp->b3_seg[i].B[2].y * dJ32dt + bp->b3_seg[i].B[3].y * dJ33dt;
   dxyzdt->z = bp->b3_seg[i].B[0].z * dJ30dt + bp->b3_seg[i].B[1].z * dJ31dt +
       bp->b3_seg[i].B[2].z * dJ32dt + bp->b3_seg[i].B[3].z * dJ33dt;

   /*
    * Rescale to get derivatives wrt to the t parameter that ranges 
    * 0.0 to 1.0 over the whole spline.
    */
   dxyzdt->x /= dt_segment;
   dxyzdt->y /= dt_segment;
   dxyzdt->z /= dt_segment;

   return (0);
} /* end function deriv_bezier_3_poly */

/*---------------------------------------------------------------*/

int bezier_3_spline( struct bezier_3_poly_data *bp,
                     int m,
                     struct point_3D p[] ) {
   /* Purpose...
    * -------
    * Given m+1 interpolation points, determine the m-segment
    * Bezier polyline that interpolates these points as a spline. 
    * This is done by first determining the array of weight points
    * which define the spline and then evaluating the cubic 
    * Bezier segments.
    *
    * Input...
    * -----
    * bp  : (pointer to) the Bezier polyline data
    *       bp needs to be allocated before calling this function.
    * m   : number of segments
    * p[] : array of m+1 data points (i=0...m)
    *
    * Output...
    * ------
    * If successful, data is written directly to the Bezier polyline 
    * data structure.
    * function returns an integer status flag:
    *  0 == normal return
    * -1 == NULL pointer to Bezier polyline
    * -2 == invalid value for m
    * -3 == NULL pointer for interpolation points
    * -4 == could not allocate weight-point array
    * -5 == m larger than n_max originally allocated
    *
    * Reference:
    * ---------
    * G. Engelin & F. Uhlig (1996)
    * Numerical Algorithms with C
    * Springer, Berlin
    * Section 12.3.1
    *
    */

   struct point_3D b0, b1, b2, b3;
   struct point_3D *d;  /* pointer to weight-point array */
   int             i, j, flag;
   double          tolerance;
   double          old_x, old_y, old_z;
   double          diff, max_diff, dx, dy, dz;

   flag = 0; /* Assume OK */

   if ( bp == NULL ) {
      flag = -1;
      return flag;
   } /* end if */

   if ( m <= 0 ) {
      flag = -2;
      return flag;
   } /* end if */

   if ( p == NULL ) {
      flag = -3;
      return flag;
   } /* end if */

   if ( m > bp->n_max ) {
      flag = -5;
      return flag;
   } /* end if */

   /*
    * Set a nominal tolerance.
    * Will be used later to terminate iterations.
    */
   tolerance = 1.0e-10;

   /*
    * Allocate memory for weight points.
    */
   d = NULL;
   d = malloc( (m+1) * sizeof(struct point_3D) );
   if ( d == NULL ) {
      flag = -4;
      return flag;
   } /* end if */

   /*
    * For a natural spline, the first and last weight points
    * are also the first and last interpolation points.
    */
   d[0].x = p[0].x;
   d[0].y = p[0].y;
   d[0].z = p[0].z;
   d[m].x = p[m].x;
   d[m].y = p[m].y;
   d[m].z = p[m].z;

   /*
    * For the initial guess at the weight points,
    * just use the supplied data points.
    */
   for (i = 1; i < m; i++) {
      d[i].x = p[i].x;
      d[i].y = p[i].y;
      d[i].z = p[i].z;
   } /* end for */

   /*
    * Apply Gauss-Seidel iteration until
    * the internal weight points converge.
    */
   for (j = 1; j < 50; j++) {
      max_diff = 0.0;
      for (i = 1; i < m; i++) {
         old_x = d[i].x;
         old_y = d[i].y;
         old_z = d[i].z;

         d[i].x = 0.25 * (6.0 * p[i].x - d[i-1].x - d[i+1].x);
         d[i].y = 0.25 * (6.0 * p[i].y - d[i-1].y - d[i+1].y);
         d[i].z = 0.25 * (6.0 * p[i].z - d[i-1].z - d[i+1].z);

         dx = fabs(d[i].x - old_x);
         dy = fabs(d[i].y - old_y);
         dz = fabs(d[i].z - old_z);
         diff = dx + dy + dz;
         if ( diff > max_diff ) max_diff = diff;
      } /* end for i */

      if ( max_diff < tolerance ) break;
   } /* end for j */

   #if 0
      printf( "Loop terminated at j=%d\n", j);
   #endif

   /*
    * Final stage; calculate the Bezier segments and pack them away.
    * Remember to set the number of segments in case it is smaller
    * than the spaces originally allocated.
    */
   bp->n = m;
   for (i = 0; i < m; i++) {
      b0.x = p[i].x;
      b0.y = p[i].y;
      b0.z = p[i].z;

      b1.x = (2.0 * d[i].x + d[i+1].x) / 3.0;
      b1.y = (2.0 * d[i].y + d[i+1].y) / 3.0;
      b1.z = (2.0 * d[i].z + d[i+1].z) / 3.0;

      b2.x = (d[i].x + 2.0 * d[i+1].x) / 3.0;
      b2.y = (d[i].y + 2.0 * d[i+1].y) / 3.0;
      b2.z = (d[i].z + 2.0 * d[i+1].z) / 3.0;

      b3.x = p[i+1].x;
      b3.y = p[i+1].y;
      b3.z = p[i+1].z;

      segment_bezier_3_poly(bp, &b0, &b1, &b2, &b3, i);
   } /* end for i */

   /*
    * Adjust the parameter break-points so that progression
    * along the spline is smooth.
    */
   normalize_bezier_3_poly (bp);

   /*
    * Clean up.
    */
   if ( d != NULL ) {
      free( d );
      d = NULL;
   } /* end if */

   return flag;
} /* end function bezier_3_spline */

/*---------------------------------------------------------------*/

int coons_bezier_3 (struct bezier_3_poly_data *c1,
		    struct bezier_3_poly_data *c2,
		    struct bezier_3_poly_data *c3,
		    struct bezier_3_poly_data *c4,
		    double r, double s,
		    struct point_3D *d)
/* Purpose...
 * -------
 * Compute a location on the Coons patch at parameter position
 * (r,s).
 *
 * Input...
 * -----
 * c1, c2, c3, c4 : pointers to the bounding curves
 *	These are normalized Bezier polylines.
 *	c1 and c2 are the South	and North boundaries respectively.
 *	The parameter 0 <= r <= 1 traverses them West to East.
 *	c3 and c4 are the West and East boundaries respectively.
 *	Parameter 0 <= s <= 1 traverses them South to North.
 * r  : West to East parameter, 0 <= r <= 1.
 * s  : South to North parameter, 0 <= s <= 1.
 *
 * Output...
 * ------
 * d  : pointer to the computed location in (x,y,z)-space
 * This function returns 0 for a normal return, -1 otherwise.
 *
 */

{  /* begin coons_bezier_3()  */

struct point_3D c1_r, c2_r, c3_s, c4_s;
struct point_3D c1_0, c1_1, c2_0, c2_1;

/* <<<<<<<< need to put error checking in >>>>>>>> */

if (r < 0.0) r = 0.0;
if (r > 1.0) r = 1.0;
if (s < 0.0) s = 0.0;
if (s > 1.0) s = 1.0;

eval_bezier_3_poly (c1, 0.0, &c1_0);
eval_bezier_3_poly (c1, 1.0, &c1_1);

eval_bezier_3_poly (c2, 0.0, &c2_0);
eval_bezier_3_poly (c2, 1.0, &c2_1);

eval_bezier_3_poly (c1, r, &c1_r);
eval_bezier_3_poly (c2, r, &c2_r);

eval_bezier_3_poly (c3, s, &c3_s);
eval_bezier_3_poly (c4, s, &c4_s);

/*
 * Now interpolate to get the location of the new point.
 */
d->x = (1.0 - s) * c1_r.x + s * c2_r.x
     + (1.0 - r) * c3_s.x + r * c4_s.x
     - (1.0 - s) * (1.0 - r) * c1_0.x
     - (1.0 - s) * r * c1_1.x
     - s * (1.0 - r) * c2_0.x
     - s * r * c2_1.x ;

d->y = (1.0 - s) * c1_r.y + s * c2_r.y
     + (1.0 - r) * c3_s.y + r * c4_s.y
     - (1.0 - s) * (1.0 - r) * c1_0.y
     - (1.0 - s) * r * c1_1.y
     - s * (1.0 - r) * c2_0.y
     - s * r * c2_1.y ;

d->z = (1.0 - s) * c1_r.z + s * c2_r.z
     + (1.0 - r) * c3_s.z + r * c4_s.z
     - (1.0 - s) * (1.0 - r) * c1_0.z
     - (1.0 - s) * r * c1_1.z
     - s * (1.0 - r) * c2_0.z
     - s * r * c2_1.z ;

return (0);
}


/*============== end of file moc_bezier.c ================*/

