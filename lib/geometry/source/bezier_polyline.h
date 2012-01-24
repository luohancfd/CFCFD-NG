/** \file bezier_polyline.h
 * \ingroup geom
 * \brief Redundant functions for Bezier polylines -- kept, just in case...
 *
 */

#include "stdlib.h"
#include "geom.h"
#include "bezier.h"


/** \brief Data for a polyline constructed from a number of
 *         third-order Bezier curves.
 */
struct bezier_3_poly_data
{
   int n;                          /* number of Bezier segments     */
   struct bezier_3_data *b3_seg;   /* data for each of the segments */
   double *t_star;                 /* Normalized distance along     */
				   /* the polyline                  */
   int closed;                     /* =0, open curve                */
				   /* =1, closed curve              */
};

/* Function Declarations */

int alloc_bezier_3_poly(struct bezier_3_poly_data *bp, int n);
int free_bezier_3_poly(struct bezier_3_poly_data *bp );

int segment_bezier_3_poly(struct bezier_3_poly_data *bp,
			  struct point_3D *loc0,
			  struct point_3D *loc1,
			  struct point_3D *loc2,
			  struct point_3D *loc3,
			  int i);

int normalize_bezier_3_poly(struct bezier_3_poly_data *bp);

int eval_bezier_3_poly(struct bezier_3_poly_data *bp,
			double t_star,
			struct point_3D *loc);

int bezier_3_spline(struct bezier_3_poly_data *bp,
		    int m,
		    struct point_3D p[] );

int coons_bezier_3(struct bezier_3_poly_data *c1,
		   struct bezier_3_poly_data *c2,
		   struct bezier_3_poly_data *c3,
		   struct bezier_3_poly_data *c4,
		   double r, double s,
		   struct point_3D *d);

