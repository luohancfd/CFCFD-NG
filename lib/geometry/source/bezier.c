/** \file bezier.c
 * \ingroup geom
 * \brief Bezier curves of third order
 *
 * These routines provide support for Bezier curves in three
 * dimensions.  Common uses are CAD and font outlines.
 *
 * \author PA Jacobs
 *
 * \version 1.0  : basic routines for third-order curves
 * \version 1.1  : 06-nov-92 : full function prototypes
 * \version 1.2  : 15-nov-92 : trap underflow in eval_bezier_3_poly()
 * \version 1.3  : 18-nov-92 : fix the polyline evaluation routine so that the
 *                    t_star array is never addressed incorrectly
 *                    (when we have only one segment in the polyline)
 * \version 1.4  : 22-Jan-93 : add slope calculation deriv_bezier_3()
 * \version 1.5  : 25-May-93 : Changed location_3D to point_3D.
 * \version 2.0  : 21-Feb-95 : Added N-degree Bezier curves
 * \version 2.1  : 11-Mar-95 : Compute the length of the Bezier segments by evaluating
 *                    points along the segment (rather than summing the control
 *                    point segments).
 * \version 2.2  : 05-Jul-98 : Cubic spline added.
 * \version 3.0  : 03-Sep-04 : Go back to simple arrays of point_3D to hold the
 *                             control points.
 *
 * References...
 *
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
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "geom.h"
#include "bezier.h"

/*--------------------------------------------------------------*/

/** \brief Evaluate the third-order Bezier curve position with
 *         parameter t.
 *
 * \param B  : (pointer to) the array of control points
 * \param t  : position parameter, 0 <= t <= 1.0
 * \param loc : (pointer to) the computed location in 3D
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 */
int bezier_3_eval(struct point_3D B[],
		  double t,
		  struct point_3D *loc)
{
    double one_m_t, J30, J31, J32, J33;
    if (t < -0.01 || t > 1.01) {
	return -1;
    }
    one_m_t = 1.0 - t;
    J30 = one_m_t * one_m_t * one_m_t;
    J31 = 3.0 * t * one_m_t * one_m_t;
    J32 = 3.0 * t * t * one_m_t;
    J33 = t * t * t;
    loc->x = B[0].x * J30 + B[1].x * J31 + B[2].x * J32 + B[3].x * J33;
    loc->y = B[0].y * J30 + B[1].y * J31 + B[2].y * J32 + B[3].y * J33;
    loc->z = B[0].z * J30 + B[1].z * J31 + B[2].z * J32 + B[3].z * J33;
    return 0;
}


/** \brief Evaluate the derivative of the third-order Bezier curve
 *         with respect to parameter t.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param B  : (pointer to) the array of control points
 * \param t  : position parameter, 0 <= t <= 1.0
 * \param dxyzdt : (pointer to) the computed derivative in 3D
 */
int bezier_3_deriv(struct point_3D B[],
		   double t,
		   struct point_3D *dxyzdt)
{
    double one_m_t, dJ30dt, dJ31dt, dJ32dt, dJ33dt;

    if (t < -0.01 || t > 1.01) {
	return -1;
    }
    one_m_t = 1.0 - t;
    dJ30dt = -3.0 * one_m_t * one_m_t;
    dJ31dt = 3.0 * one_m_t * one_m_t - 6.0 * t * one_m_t;
    dJ32dt = 6.0 * t * one_m_t - 3.0 * t * t;
    dJ33dt = 3.0 * t * t;
    dxyzdt->x = B[0].x * dJ30dt + B[1].x * dJ31dt +
	        B[2].x * dJ32dt + B[3].x * dJ33dt;
    dxyzdt->y = B[0].y * dJ30dt + B[1].y * dJ31dt +
	        B[2].y * dJ32dt + B[3].y * dJ33dt;
    dxyzdt->z = B[0].z * dJ30dt + B[1].z * dJ31dt +
	        B[2].z * dJ32dt + B[3].z * dJ33dt;
    return 0;
} /* end bezier_3_deriv() */


/** \brief Evaluate the n-degree Bezier curve position with
 *         parameter t.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * If n==3, work is delegated to the bezier_3_eval(), otherwise
 * the de Casteljau algorithm is used.
 *
 * \param n : degree of the Bezier curve
 *            (must be less than MAX_BEZ_DEGREE else it will truncate)
 *            Also, note that there are n+1 polygon points for an
 *            n-degree Bezier curve.
 * \param B   : (pointer to) the array of control points
 * \param t   : position parameter, 0 <= t <= 1.0
 * \param loc : (pointer to) the computed location in 3D
 */
int bezier_eval(int n,
		struct point_3D B[],
		double t,
		struct point_3D *loc)
{
    double one_m_t;
    struct point_3D Q[MAX_BEZ_DEGREE+1];
    int    k, i;

    if (t < -0.01 || t > 1.01) {
	return -1;
    }

    if ( n == 3 ) {
	/* Delegate the work. */
	return bezier_3_eval( B, t, loc );
    } else {
	/* Do the real work here. */
	one_m_t = 1.0 - t;
	/* Copy the control points into the work array. */
	for (i = 0; i <= n; ++i) {
	    Q[i].x = B[i].x;
	    Q[i].y = B[i].y;
	    Q[i].z = B[i].z;
	}
	/* Now, generate one new level at a time;
	 * over-writing the work array at each level. */
	for (k = 1; k <= n; ++k) {
	    for (i = 0; i <= n - k; ++i) {
		Q[i].x = one_m_t * Q[i].x + t * Q[i+1].x;
		Q[i].y = one_m_t * Q[i].y + t * Q[i+1].y;
		Q[i].z = one_m_t * Q[i].z + t * Q[i+1].z;
	    }
	}  /* for (k... */
	/* The point on the curve is now in Q[0] */
	loc->x = Q[0].x;
	loc->y = Q[0].y;
	loc->z = Q[0].z;
	return 0;
    }
} /* end bezier_eval() */


/** \brief Translates a Bezier curve by (dx, dy, dz) in Caresian space. 
 *
 * \param n : order of the Bezier curve (there should be n+1 points in B)
 * \param B : array of control points
 * \param dx, dy, dz : displacement in Cartesian coordinates
 */
int bezier_translate( int n, struct point_3D B[],
		      double dx, double dy, double dz )
{
    int i;
    if ( B == NULL ) return -1;
    for ( i = 0; i <= n; ++i ) {
	point_3D_translate( &(B[i]), dx, dy, dz );
    }
    return 0;
}


/** \brief Returns an estimate of the length of a Bezier curve. 
 *
 * \param B : array of control points
 * \param n : order of the Bezier curve (there should be n+1 points in B)
 */
double bezier_length( int n, struct point_3D B[] )
{
    double length, t, dx, dy, dz;
    double xa, ya, za, xb, yb, zb;
    struct point_3D loc;
    length = 0.0;
    xa = B[0].x; ya = B[0].y; za = B[0].z;
    /* Step along the segment, summing the arc length as we go. */
    for ( t = 0.1; t <= 1.001; t += 0.1) {
	bezier_eval( n, B, t, &loc );
	xb = loc.x; yb = loc.y; zb = loc.z;
	dx = xb - xa; dy = yb - ya; dz = zb - za;
	length += sqrt( dx * dx + dy * dy + dz * dz );
	xa = xb; ya = yb; za = zb;  /* move the left point along */
    } /* endfor */
    return length;
}


/** \brief Given m+1 interpolation points, determine the m-segment
 *         Bezier polyline that interpolates these points as a spline. 
 *
 * This is done by first determining the array of weight points
 * which define the spline and then evaluating the cubic 
 * Bezier segments.
 *
 * If successful, data is written back into the array for the Bezier
 * control points. 
 *
 * function returns an integer status flag:
 *  0 == normal return
 * -1 == NULL pointer to the array for the Bezier control points
 * -2 == invalid value for m
 * -3 == NULL pointer for interpolation points
 * -4 == could not allocate weight-point array
 *
 * \param bcp[] : the Bezier control points array
 *                bp needs to be a valid place to store m+1
 *                point_3D structures calling this function.
 * \param m   : number of Bezier-3 segments
 * \param p[] : array of m+1 data points (i=0...m)
 *
 * Reference:
 *
 * G. Engelin & F. Uhlig (1996)
 * Numerical Algorithms with C
 * Springer, Berlin
 * Section 12.3.1
 */
int bezier_3_spline( struct point_3D bcp[],
                     int m,
                     struct point_3D p[] ) 
{
    struct point_3D *d;  /* pointer to weight-point array */
    int             i, j;
    double          tolerance;
    double          old_x, old_y, old_z;
    double          diff, max_diff, dx, dy, dz;
    if ( bcp == NULL ) {
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
    d = (struct point_3D*) malloc( (m+1) * sizeof(struct point_3D) );
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
#   if 0
    printf( "Loop terminated at j=%d\n", j);
#   endif
    /*
     * Final stage; calculate the Bezier segments and pack them away.
     */
    for (i = 0; i < m; i++) {
	bcp[i*4 + 0].x = p[i].x;
	bcp[i*4 + 0].y = p[i].y;
	bcp[i*4 + 0].z = p[i].z;
	bcp[i*4 + 1].x = (2.0 * d[i].x + d[i+1].x) / 3.0;
	bcp[i*4 + 1].y = (2.0 * d[i].y + d[i+1].y) / 3.0;
	bcp[i*4 + 1].z = (2.0 * d[i].z + d[i+1].z) / 3.0;
	bcp[i*4 + 2].x = (d[i].x + 2.0 * d[i+1].x) / 3.0;
	bcp[i*4 + 2].y = (d[i].y + 2.0 * d[i+1].y) / 3.0;
	bcp[i*4 + 2].z = (d[i].z + 2.0 * d[i+1].z) / 3.0;
	bcp[i*4 + 3].x = p[i+1].x;
	bcp[i*4 + 3].y = p[i+1].y;
	bcp[i*4 + 3].z = p[i+1].z;
    } /* end for i */
    
    /* Clean up. */
    if ( d != NULL ) {
	free( d );
	d = NULL;
    }
    return 0;
} /* end function bezier_3_spline() */


/*--------------------- end of file ----------------------*/

