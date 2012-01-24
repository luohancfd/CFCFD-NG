/** \file bezier.h
 * \ingroup geom
 * \brief Bezier curves of third order and arbitrary order -- header file.
 *
 */

#ifndef BEZIER_H
#define BEZIER_H  1

#include "geom.h"

int bezier_3_eval(struct point_3D B[],
		  double t,
		  struct point_3D *loc);
int bezier_3_deriv(struct point_3D B[],
		   double t,
		   struct point_3D *dxyzdt);

#define  MAX_BEZ_DEGREE  20

int bezier_eval(int n,
		struct point_3D *B,
		double t,
		struct point_3D *loc);

int bezier_translate( int n, struct point_3D B[],
		      double dx, double dy, double dz );

double bezier_length( int n, struct point_3D B[] );

int bezier_3_spline( struct point_3D bcp[],
                     int m,
                     struct point_3D p[] );

#endif


