/** \file bezier_polyline.c
 * \ingroup geom
 * \brief Redundant functions for Bezier polylines -- kept, just in case...
 *
 * \author PA Jacobs
 */

#include <math.h>
#include "geom.h"
#include "bezier.h"
#include "bezier_polyline.h"


/** \brief Allocate memory for the Bezier polyline data.
 *
 * \param bp : pointer to the polyline data structure
 * \param n  : number of Bezier segments in the polyline
 *
 * Returns -1 if the allocation fails, 0 otherwise.
 */
int alloc_bezier_3_poly(struct bezier_3_poly_data *bp, int n)
{
    int i;
    if (n <= 0) {
	return -1;
    }
    bp->b3_seg = (struct bezier_3_data *)
	malloc (n * sizeof(struct bezier_3_data) );
    if (bp->b3_seg == NULL) {
	return -1;
    }
    bp->t_star = (double *) malloc ( n * sizeof(double) );
    if (bp->t_star == NULL) {
	return -1;
    }
    for (i = 0; i < n; ++i) {
	/* Set all of the orders to 3 */
	bp->b3_seg[i].n = 3;
	/* Set the normalized distance along the polyline to
	 * nominal values */
	bp->t_star[i] = (double) (i + 1.0) / n;
    }
    /* This is the assumed number of Bezier elements. */
    bp->n = n;
    return 0;
}

/** \brief Release the memory that had been allocated to the polyline. */
int free_bezier_3_poly(struct bezier_3_poly_data *bp )
{
    if ( bp->b3_seg != NULL ) free( bp->b3_seg );
    if ( bp->t_star != NULL ) free( bp->t_star );
    return 0;
}


/** \brief Initialize the coordinates of the Bezier control points for
 *         Bezier segment i.
 *
 * \param p3     : (pointer to) the Bezier polyline data
 * \param loc0-3 : (pointer to) the location of control points 0-3
 * \param i      : index of the desired segment, 0 <= i < n
 */
int segment_bezier_3_poly(struct bezier_3_poly_data *bp,
			  struct point_3D *loc0,
			  struct point_3D *loc1,
			  struct point_3D *loc2,
			  struct point_3D *loc3,
			  int i)
{
    if (i < 0 || i >= bp->n) {
	return -1;
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
    return 0;
} /* end segment_bezier_3_poly() */


/** \brief Normalize the Bezier polyline so that the parameter t_star
 *         varies from 0.0 to 1.0 along the entire curve.
 *
 * Returns -1 if the normalization fails, 0 otherwise.
 *
 * Arc length along the Bezier segments is approximated as the
 * sum of the distances between a number of points along each
 * curve.
 *
 * \param bp : pointer to the polyline data structure
 */
int normalize_bezier_3_poly(struct bezier_3_poly_data *bp)
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
    for (i = 0; i < bp->n; ++i) {
	B3 = &(bp->b3_seg[i]);
	seg_length = 0.0;
	p = &(B3->B[0]);  /* location of first control point on the segment */
	xa = p->x; ya = p->y; za = p->z;
	/* Step along the segment, summing the arc length as we go. */
	for ( t = 0.1; t <= 1.001; t += 0.1) {
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
    if (arc_length <= 0.0) {
	return -1;
    }
    /*
     * Normalize the distance along the polyline so that
     * t_star varies from 0.0 to 1.0 along the entire polyline.
     * If there is only one segment, then t_star[0] will be 1.0
     * and no other segments will be valid.
     */
    for (i = 0; i < bp->n; ++i) bp->t_star[i] /= arc_length;
    return 0;
} /* end normalize_bezier_3_poly() */


/** \brief Evaluate the position on the Bezier polyline at given
 *         parameter t_star.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param bp     : (pointer to) the Bezier polyline data
 * \param t_star : position parameter, 0 <= t_star <= 1.0
 * \param loc    : (pointer to) the computed location in 3D
 */
int eval_bezier_3_poly(struct bezier_3_poly_data *bp,
		       double t_star,
		       struct point_3D *loc)
{
    double t, one_m_t, J30, J31, J32, J33;
    int    i;

    /* Warning: The following default behaviour may be unsuitable
     *  for some applications. */
    if (t_star < 0.0) t_star = 0.0;
    if (t_star > 1.0) t_star = 1.0;

    /* Find the relevant segment. */
    i = 0;
    while (t_star > bp->t_star[i] && i < (bp->n - 1) ) ++i;

    /* Compute the local parameter. */
    if (i == 0) {
	t = t_star / bp->t_star[i];
    } else {
	t = (t_star - bp->t_star[i-1]) /
	    (bp->t_star[i] - bp->t_star[i-1]);
    }

    /* Now compute the position on the identified Bezier segment. */
    one_m_t = 1.0 - t;
    J30 = one_m_t * one_m_t * one_m_t;
    J31 = 3.0 * t * one_m_t * one_m_t;
    J32 = 3.0 * t * t * one_m_t;
    J33 = t * t * t;
    if (t < 1.0e-6) {
	loc->x = bp->b3_seg[i].B[0].x;
	loc->y = bp->b3_seg[i].B[0].y;
	loc->z = bp->b3_seg[i].B[0].z;
    } else if (one_m_t < 1.0e-6) {
	loc->x = bp->b3_seg[i].B[3].x;
	loc->y = bp->b3_seg[i].B[3].y;
	loc->z = bp->b3_seg[i].B[3].z;
    } else {
	loc->x = bp->b3_seg[i].B[0].x * J30 + bp->b3_seg[i].B[1].x * J31 +
	    bp->b3_seg[i].B[2].x * J32 + bp->b3_seg[i].B[3].x * J33;
	loc->y = bp->b3_seg[i].B[0].y * J30 + bp->b3_seg[i].B[1].y * J31 +
	    bp->b3_seg[i].B[2].y * J32 + bp->b3_seg[i].B[3].y * J33;
	loc->z = bp->b3_seg[i].B[0].z * J30 + bp->b3_seg[i].B[1].z * J31 +
	    bp->b3_seg[i].B[2].z * J32 + bp->b3_seg[i].B[3].z * J33;
    }
    return 0;
} /* end eval_bezier_3_poly() */


/** \brief Given m+1 interpolation points, determine the m-segment
 * Bezier polyline that interpolates these points as a spline. 
 *
 * This is done by first determining the array of weight points
 * which define the spline and then evaluating the cubic 
 * Bezier segments.
 *
 * If successful, data is written directly to the Bezier polyline 
 * data structure.
 * function returns an integer status flag:
 *  0 == normal return
 * -1 == NULL pointer to Bezier polyline
 * -2 == invalid value for m
 * -3 == NULL pointer for interpolation points
 * -4 == could not allocate weight-point array
 *
 * \param bp  : (pointer to) the Bezier polyline data
 *              bp needs to be allocated before calling this function.
 * \param m   : number of segments
 * \param p[] : array of m+1 data points (i=0...m)
 *
 * Reference:
 *
 * G. Engelin & F. Uhlig (1996)
 * Numerical Algorithms with C
 * Springer, Berlin
 * Section 12.3.1
 */
int bezier_3_spline( struct bezier_3_poly_data *bp,
                     int m,
                     struct point_3D p[] ) {
   struct point_3D b0, b1, b2, b3;
   struct point_3D *d;  /* pointer to weight-point array */
   int             i, j;
   double          tolerance;
   double          old_x, old_y, old_z;
   double          diff, max_diff, dx, dy, dz;
   if ( bp == NULL ) {
      return -1;
   }
   if ( m <= 0 ) {
      return -2;
   }
   if ( p == NULL ) {
      return -3;
   }
   /* Allocate memory for weight points. */
   d = NULL;
   d = malloc( (m+1) * sizeof(struct point_3D) );
   if ( d == NULL ) {
      return -4;
   }
   /* Set a nominal tolerance.
    * Will be used later to terminate iterations. */
   tolerance = 1.0e-10;
   /* For a natural spline, the first and last weight points
    * are also the first and last interpolation points. */
   d[0].x = p[0].x;
   d[0].y = p[0].y;
   d[0].z = p[0].z;
   d[m].x = p[m].x;
   d[m].y = p[m].y;
   d[m].z = p[m].z;
   /* For the initial guess at the weight points,
    * just use the supplied data points. */
   for (i = 1; i < m; i++) {
      d[i].x = p[i].x;
      d[i].y = p[i].y;
      d[i].z = p[i].z;
   }
   /* Apply Gauss-Seidel iteration until
    * the internal weight points converge. */
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
    */
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

   /* Adjust the parameter break-points so that progression
    * along the spline is smooth. */
   normalize_bezier_3_poly (bp);

   /* Clean up. */
   if ( d != NULL ) {
      free( d );
      d = NULL;
   }
   return 0;
} /* end function bezier_3_spline() */


/** \brief Compute a location on the Coons patch at parameter position (r,s).
 *
 * Returns 0 for a normal return, -1 otherwise.
 *
 * \param c1, c2, c3, c4 : pointers to the bounding curves
 *	These are normalized Bezier polylines.
 *	c1 and c2 are the South	and North boundaries respectively.
 *	The parameter 0 <= r <= 1 traverses them West to East.
 *	c3 and c4 are the West and East boundaries respectively.
 *	Parameter 0 <= s <= 1 traverses them South to North.
 * \param r  : West to East parameter, 0 <= r <= 1.
 * \param s  : South to North parameter, 0 <= s <= 1.
 * \param d  : pointer to the computed location in (x,y,z)-space
 */
int coons_bezier_3(struct bezier_3_poly_data *c1,
		   struct bezier_3_poly_data *c2,
		   struct bezier_3_poly_data *c3,
		   struct bezier_3_poly_data *c4,
		   double r, double s,
		   struct point_3D *d)
{
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
    /* Now interpolate to get the location of the new point.  */
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
    return 0;
} /* end coons_bezier_3() */

/*-----------------------------------------------------------------*/
