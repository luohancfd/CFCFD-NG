/** \file gpath.h
 * \ingroup geom
 * \brief Functions defining the various types of curves used for
 *        specifying flow geometry.
 *
 */

#ifndef GPATH_H_ALREADY_INCLUDED

#include "stdlib.h"
#include "geom.h"
#include "bezier.h"

#define GPATH_NONE   0
#define GPATH_LINE   1
#define GPATH_ARC    2
#define GPATH_BEZIER 3

/** \brief Path elements may be lines, arcs or Bezier curves. */
struct GPathElement
{
    int type; /**< one of GPATH_LINE, GPATH_ARC or GPATH_BEZIER */
    int np; /**< number of points in the array */
    struct Vector3D *p; /**< pointer to the array of defining points */
};

#define GPATH_MAX_ELEMENTS 100

/** \brief Data for a polyline constructed from a number of elements. */
struct GPathPolyLine
{
    int ne; /**< number of elements */
    struct GPathElement pe[GPATH_MAX_ELEMENTS]; /**< element array */
    double el[GPATH_MAX_ELEMENTS]; /**< element length array */
    double length; /**< current length, updated after adding an element */
    double t_star[GPATH_MAX_ELEMENTS]; /**< t-value at end of element */
    double t0, t1; /**< evaluation subrange, dafaults to 0.0, 1.0  */
    int closed; /**< 0= open curve; 1= closed curve */
};


struct GPathPolyLine *gpath_init(struct GPathPolyLine *plp);

int gpath_add_element(struct GPathPolyLine *plp,
		      int type, int np, struct point_3D p[],
		      int direction);
double gpath_element_length(struct GPathElement *pe);
int gpath_write_element_to_file( FILE *fp, struct GPathElement *pe );
int gpath_write_all_elements_to_file( FILE *fp, struct GPathPolyLine *plp );
int gpath_scan_from_string_and_add_element( char str[], struct GPathPolyLine *plp );

int gpath_append_polyline(struct GPathPolyLine *plp1,
			  struct GPathPolyLine *plp2,
			  int direction);

struct point_3D *gpath_get_first_point_on_path( struct GPathPolyLine *plp );
struct point_3D *gpath_get_last_point_on_path( struct GPathPolyLine *plp );

int gpath_polyline_set_subrange(struct GPathPolyLine *plp,
			double t0, double t1);
int gpath_polyline_eval(struct GPathPolyLine *plp,
			double t_in_subrange,
			struct point_3D *loc);
int gpath_polyline_translate(struct GPathPolyLine *plp,
			     double dx, double dy, double dz);

int coons_patch(struct GPathPolyLine *c1,
		   struct GPathPolyLine *c2,
		   struct GPathPolyLine *c3,
		   struct GPathPolyLine *c4,
		   double r, double s,
		   struct point_3D *d);
int TFI_3D_no_array(struct GPathPolyLine *c0,
		    struct GPathPolyLine *c1,
		    struct GPathPolyLine *c2,
		    struct GPathPolyLine *c3,
		    struct GPathPolyLine *c4,
		    struct GPathPolyLine *c5,
		    struct GPathPolyLine *c6,
		    struct GPathPolyLine *c7,
		    struct GPathPolyLine *c8,
		    struct GPathPolyLine *c9,
		    struct GPathPolyLine *c10,
		    struct GPathPolyLine *c11,
		    double r, double s, double t,
		    struct point_3D *d);
int TFI_3D(struct GPathPolyLine *c[],
	   double r, double s, double t,
	   struct point_3D *d);

int line_translate(struct point_3D B[], double dx, double dy, double dz);
int line_eval(struct point_3D B[], double t, struct point_3D *loc);
double line_length(struct point_3D B[]);

int arc_translate(struct point_3D B[], double dx, double dy, double dz);
int arc_eval_position_and_length(struct point_3D B[], 
				 double t, 
				 struct point_3D *loc, 
				 double *length);
int arc_eval(struct point_3D B[], double t, struct point_3D *loc);
double arc_length(struct point_3D B[]);


#define GPATH_H_ALREADY_INCLUDED  1
#endif
