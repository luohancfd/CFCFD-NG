/** \file gpath.c
 * \ingroup geom
 * \brief Functions defining the various types of curves used for
 *        specifying flow geometry.
 *
 * \author PA Jacobs
 *
 * \version 1.0  : adapted from bezier.c and extended 
 *                 with straight lines and arcs
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geom.h"
#include "bezier.h"
#include "gpath.h"

/** \brief Initializes data structure members and returns a pointer.
 *
 * \param plp : pointer to a GPathPolyline structure.
 */
struct GPathPolyLine *gpath_init(struct GPathPolyLine *plp)
{
    if ( plp == NULL ) {
	/* printf("gpath_init: Creating new polyline\n"); */
	plp = (struct GPathPolyLine *) calloc(sizeof(struct GPathPolyLine), 1);
    }
    if ( plp != NULL ) {
	plp->ne = 0;       /**< the path is initially empty          */
	plp->length = 0.0; /**< length should be always current      */
	plp->t0 = 0.0;     /**< default subrange is the whole path   */
	plp->t1 = 1.0;
	plp->closed = 0;   /**< == 0 indicates that the path is open */
    }
    return plp;
}


/** \brief Adds a path element to a polyline.
 *
 * Returns the new number of elements on success, -1 otherwise.
 *
 * Presently, these polylines may consist of disjoint segments.
 * Maybe we should test that the ends of adjacent segments are
 * coincident; maybe we shouldn't.
 *
 * \param plp  : pointer to the polyline
 * \param type : type of element as listed in the gpath.h
 * \param np   : number of control points in the new element.
 * \param p    : the array of control points
 * \param direction : +1 add in points normal order; -1 add in reverse order
 */
int gpath_add_element(struct GPathPolyLine *plp,
		      int type, int np, struct point_3D p[], 
		      int direction)
{
    int i, jp, je, np_expected;
    double partial_length;

    i = plp->ne;  /* index of the new element */
    if ( plp->ne >= GPATH_MAX_ELEMENTS ) {
	printf( "gpath_add_element(): too many elements; %d already\n", plp->ne );
	return -1;
    }
    /* We have space, add the data for the new element */
    plp->pe[i].type = type;

    /* For lines and arcs, we expect particular numbers of points. */
    if ( type == GPATH_LINE ) {
	np_expected = 2;
    } else if ( type == GPATH_ARC ) {
	np_expected = 3;
    } else {
	np_expected = np;  /* accept the supplied input as correct */
    }
    if ( np_expected != np ) {
	printf( "gpath_add_element(): did not get the expected number of points\n" );
        return -1;
    }

    plp->pe[i].p = (struct point_3D *) calloc(np, sizeof(struct point_3D));
    if ( plp->pe[i].p == NULL ) {
	printf("gpath_add_element(): memory could not be allocated for points\n");
	return -1;
    } else {
	plp->pe[i].np = np;
	if ( direction == 1 ) {
	    /* Copy the control points forward direction. */
	    for ( jp = 0; jp < np; ++jp ) {
		vector_copy_3D( &(p[jp]), &(plp->pe[i].p[jp]) );
	    }
	} else if ( direction == -1 ) {
	    /* Copy the control points, going in reverse. */
	    if ( type == GPATH_ARC ) {
		/* Reverse start and end point, but not centre */
		vector_copy_3D( &(p[0]), &(plp->pe[i].p[1]) );
		vector_copy_3D( &(p[1]), &(plp->pe[i].p[0]) );
		vector_copy_3D( &(p[2]), &(plp->pe[i].p[2]) );
	    } else {
		/* General case of copy the control points, going in reverse. */
		for ( jp = 0; jp < np; ++jp ) {
		    vector_copy_3D( &(p[jp]), &(plp->pe[i].p[np-1 - jp]) );
		}
	    }
	} else {
	    printf("gpath_add_element(): invalid direction %d\n", direction);
	    return -1;
	}
    }

    plp->el[i] = gpath_element_length( &(plp->pe[i]) );
    plp->length += plp->el[i];

    /* if everything has gone well, update number of elements */
    plp->ne += 1;

    /* Compute values for the position parameter at the end of each element. 
     * It is something like a normalised length at the end of an element
     * and can be used later to decide which element we are positioned in. */
    partial_length = 0.0;
    for ( je = 0; je < plp->ne; ++je ) {
	partial_length += plp->el[je];
	plp->t_star[je] = partial_length / plp->length;
    }

    return plp->ne;
} /* end gpath_add_element() */


/** \brief Length of the GPathElement.
 *
 * \param pe : pointer to the GPathElement
 * \returns  : the length of the element
 */
double gpath_element_length(struct GPathElement *pe)
{
    double length;
    if ( pe->type == GPATH_BEZIER ) {
	length = bezier_length( (pe->np)-1, pe->p ); /* bezier order = np-1 */
    } else if ( pe->type == GPATH_ARC ) {
	length = arc_length( pe->p );
    } else if ( pe->type == GPATH_LINE ) {
	length = line_length( pe->p );
    } else {
	printf("gpath_element_length(): invalid type of element: %d\n", pe->type);
	length = 0.0;
    }
    return length;
} /* end gpath_element_length() */


/** \brief Writes the all of the GPathElements of GPathPolyline 
 *         to the specified file.
 *
 * \param fp : pointer to the open file
 *             if NULL, stdout will be used.
 * Returns number of characters written if successful, -1 otherwise.
 */
int gpath_write_all_elements_to_file( FILE *fp, struct GPathPolyLine *plp )
{
    int ie, nc, result;
    result = 0;
    if ( fp == NULL ) fp = stdout;
    for (ie = 0; ie < plp->ne; ++ie) {
	nc = gpath_write_element_to_file(fp, &(plp->pe[ie]) );
	if (nc == -1) {
	    result = -1;
	    break;
	} else {
	    result += nc;
	}
    }
    return result;
}


/** \brief Writes the GPathElement to the specified file.
 *
 * Returns number of characters written if successful, -1 otherwise.
 */
int gpath_write_element_to_file( FILE *fp, struct GPathElement *pe )
{
    int nc, total_chars, ip;
    total_chars = 0;
    if ( fp == NULL ) {
	printf( "gpath_write_element_to_file(): NULL file handle.\n" );
	return -1;
    }
    if ( pe->type == GPATH_LINE ) {
	nc = fprintf(fp, "LINE %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", 
		     pe->p[0].x, pe->p[0].y, pe->p[0].z, 
		     pe->p[1].x, pe->p[1].y, pe->p[1].z );
	if ( nc > 0 ) total_chars += nc; else return nc;
    } else if ( pe->type == GPATH_ARC ) {
	nc = fprintf(fp, "ARC %20.12e %20.12e %20.12e ",
		     pe->p[0].x, pe->p[0].y, pe->p[0].z ); 
	if ( nc > 0 ) total_chars += nc; else return nc;
	nc = fprintf(fp, "%20.12e %20.12e %20.12e ",
		     pe->p[1].x, pe->p[1].y, pe->p[1].z );
	if ( nc > 0 ) total_chars += nc; else return nc;
	nc = fprintf(fp, "%20.12e %20.12e %20.12e\n", 
		     pe->p[2].x, pe->p[2].y, pe->p[2].z );
	if ( nc > 0 ) total_chars += nc; else return nc;
    } else if ( pe->type == GPATH_BEZIER ) {
	nc = fprintf(fp, "BEZIER %d", pe->np );
	if ( nc > 0 ) total_chars += nc; else return nc;
	for ( ip = 0; ip < pe->np; ++ip ) {
	    nc = fprintf(fp, " %20.12e %20.12e %20.12e", pe->p[ip].x, 
			 pe->p[ip].y, pe->p[ip].z );
	    if ( nc > 0 ) total_chars += nc; else return nc;
	}
	nc = fprintf(fp, "\n");
	if ( nc > 0 ) total_chars += nc; else return nc;
    } else {
	printf( "gpath_write_element_to_file(): unknown type %d\n", pe->type );
    }
    return total_chars;
} /* end gpath_write_element_to_file() */


/** \brief Scans the GPathElement from a supplied string and appends it to 
 *         the GPathPolyLine.
 *
 * Returns 0 if successful, -1 otherwise.
 */
int gpath_scan_from_string_and_add_element( char str[], struct GPathPolyLine *plp )
{
    int ip, np, count;
    char type_name[32], token[32];
    struct point_3D p[32];

    if ( str == NULL ) {
	printf( "gpath_scan_from_string_and_add_element(): NULL string.\n" );
	return -1;
    }
    count = sscanf(str, "%s", type_name);
    if ( strcmp(type_name, "LINE") == 0 ) {
	count = sscanf(str, "LINE %lf %lf %lf %lf %lf %lf",
		       &(p[0].x), &(p[0].y), &(p[0].z),
		       &(p[1].x), &(p[1].y), &(p[1].z) );
	gpath_add_element(plp, GPATH_LINE, 2, p, 1);
    } else if ( strcmp(type_name, "ARC") == 0 ) {
	count = sscanf(str, "ARC %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		       &(p[0].x), &(p[0].y), &(p[0].z),
		       &(p[1].x), &(p[1].y), &(p[1].z),
		       &(p[2].x), &(p[2].y), &(p[2].z) );
	gpath_add_element(plp, GPATH_ARC, 3, p, 1);
    } else if ( strcmp(type_name, "BEZIER") == 0 ) {
	/* We need to pull apart the line manually because we don't know np before-hand. */
	strcpy( token, strtok(str, " ") );  /* Should have been BEZIER */
	strcpy( token, strtok(NULL, " ") );
	sscanf( token, "%d", &np );
	for ( ip = 0; ip < np; ++ip ) {
	    strcpy( token, strtok(NULL, " ") );
	    sscanf( token, "%lf", &(p[ip].x) );
	    strcpy( token, strtok(NULL, " ") );
	    sscanf( token, "%lf", &(p[ip].y) );
	    strcpy( token, strtok(NULL, " ") );
	    sscanf( token, "%lf", &(p[ip].z) );
	}
	gpath_add_element(plp, GPATH_BEZIER, np, p, 1);
    } else {
	printf( "gpath_scan_from_string_and_add_element(): unknown type %s\n", type_name );
	printf( "Assuming old format: a Bezier with 4 control points and only x, y specified.\n" );
	count = sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf",
		       &(p[0].x), &(p[0].y), &(p[1].x), &(p[1].y),
		       &(p[2].x), &(p[2].y), &(p[3].x), &(p[3].y) );
	p[0].z = 0.0; p[1].z = 0.0; p[2].z = 0.0; p[3].z = 0.0;
	gpath_add_element(plp, GPATH_BEZIER, 4, p, 1);
    }
    return 0;
} /* end  gpath_scan_from_string_and_add_element() */


/** \brief Appends another path polyline to a first polyline.
 *
 * Returns the new number of elements on success, -1 otherwise.
 *
 * \param plp1 : first polyline to which the other is appended
 * \param plp2 : polyline to append to first
 *               (This one remains unchanged.)
 * \param direction : =1 appends the plp2 data in forward direction (of plp2)
 *                    =-1 appends the plp2 in reverse direction (of plp2)
 */
int gpath_append_polyline(struct GPathPolyLine *plp1,
			  struct GPathPolyLine *plp2,
			  int direction)
{
    int je;
    if ( plp1->ne >= GPATH_MAX_ELEMENTS ) {
	printf( "gpath_append_polyline(): too many elements; %d already\n", plp1->ne );
	return -1;
    }
    if ( direction == 1 ) {
	/* Append elements in the forward direction of plp2. */
	for ( je = 0; je < plp2->ne; ++je ) {
	    gpath_add_element(plp1, plp2->pe[je].type, plp2->pe[je].np, plp2->pe[je].p, 1);
	} /* end for je */
    } else if ( direction == -1 ) {
	/* Append elements in the reverse direction of plp2. */
	for ( je = plp2->ne - 1; je >= 0; --je ) {
	    gpath_add_element(plp1, plp2->pe[je].type, plp2->pe[je].np, plp2->pe[je].p, -1);
	} /* end for je */
    } else {
	printf("gpath_append_polyline(): invalid direction %d\n", direction);
	return -1;
    }

    return plp1->ne;
} /* end gpath_append_polyline() */


/** \brief Returns a pointer to the first point on a path. 
 */
struct point_3D * gpath_get_first_point_on_path( struct GPathPolyLine *plp )
{
    return &(plp->pe[0].p[0]);
}


/** \brief Returns a pointer to the last point on a path. 
 */
struct point_3D * gpath_get_last_point_on_path( struct GPathPolyLine *plp )
{
    struct GPathElement *pe;
    pe = &(plp->pe[plp->ne - 1]);
    if ( pe->type == GPATH_ARC ) {
	return &(pe->p[1]);
    } else {
	return &(pe->p[pe->np - 1]);
    }
}


/** \brief Set the parameter subrange for subsequent evaluations of position.
 *
 * \param plp : pointer to the GPathPolyline structure
 * \param t0  : start of subrange
 * \param t1  : end of subrange 
 *
 * When a polyline is first created, the subrange defaults to 0.0, 1.0.
 */
int gpath_polyline_set_subrange(struct GPathPolyLine *plp,
			double t0, double t1)
{
    if (t0 < 0.0) t0 = 0.0;
    if (t0 > 1.0) t0 = 1.0;
    if (t1 < 0.0) t1 = 0.0;
    if (t1 > 1.0) t1 = 1.0;
    plp->t0 = t0;
    plp->t1 = t1;
    return 0;
}


/** \brief Evaluate the position on the polyline at given parameter t.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param bp  : (pointer to) the Bezier polyline data
 * \param t   : position parameter, 0 <= t <= 1.0
 * \param loc : (pointer to) the computed location in 3D
 */
int gpath_polyline_eval(struct GPathPolyLine *plp,
			double t_in_subrange,
			struct point_3D *loc)
{
    int    i;
    double t, t_local;

    /* The following default behaviour may be unsuitable
     * for some applications. Consider yourself warned. */
    if (t_in_subrange < 0.0) t_in_subrange = 0.0;
    if (t_in_subrange > 1.0) t_in_subrange = 1.0;

    /* Transform from a subrange value of t. */
    t = plp->t0 + t_in_subrange * (plp->t1 - plp->t0);

    /* Find the relevant segment. */
    i = 0;
    while (t > plp->t_star[i] && i < (plp->ne - 1) ) ++i;

    /* Compute the local parameter for the selected element. */
    if (i == 0) {
	t_local = t / plp->t_star[i];
    } else {
	t_local = (t - plp->t_star[i-1]) /
	    (plp->t_star[i] - plp->t_star[i-1]);
    }

    /* Now compute the position on the selected element. */
    if ( plp->pe[i].type == GPATH_BEZIER ) {
	bezier_eval( 3, plp->pe[i].p, t_local, loc );
    } else if ( plp->pe[i].type == GPATH_ARC ) {
	arc_eval( plp->pe[i].p, t_local, loc );
    } else if ( plp->pe[i].type == GPATH_LINE ) {
	line_eval( plp->pe[i].p, t_local, loc );
    } else {
	printf("gpath_polyline_eval(): invalid type of element: %d\n",
	       plp->pe[i].type);
	return -1;
    }
    return 0;
} /* end gpath_polyline_eval() */


/** \brief Translate the polyline position by dx,dy,dz.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param plp : (pointer to) the polyline structure
 * \param dx, dy, dz
 */
int gpath_polyline_translate(struct GPathPolyLine *plp,
			     double dx, double dy, double dz)
{
    int je;
    if ( plp == NULL ) return -1;
    for ( je = 0; je < plp->ne; ++je ) {
	if ( plp->pe[je].type == GPATH_BEZIER ) {
	    bezier_translate( 3, plp->pe[je].p, dx, dy, dz );
	} else if ( plp->pe[je].type == GPATH_ARC ) {
	    arc_translate( plp->pe[je].p, dx, dy, dz );
	} else if ( plp->pe[je].type == GPATH_LINE ) {
	    line_translate( plp->pe[je].p, dx, dy, dz );
	} else {
	    printf("gpath_polyline_translate(): element %d has invalid type %d\n",
		   je, plp->pe[je].type);
	    return -1;
	}
    } /* end for je */
    return 0;
}


/** \brief Compute a location on the Coons patch at parameter position (r,s).
 *
 * Returns 0 for a normal return, -1 otherwise.
 *
 * \param c1, c2, c3, c4 : pointers to the bounding curves
 *	These are normalized GPathPolyLines.
 *	c1 and c2 are the South	and North boundaries respectively.
 *	The parameter 0 <= r <= 1 traverses them West to East.
 *	c3 and c4 are the West and East boundaries respectively.
 *	Parameter 0 <= s <= 1 traverses them South to North.
 * \param r  : West to East parameter, 0 <= r <= 1.
 * \param s  : South to North parameter, 0 <= s <= 1.
 * \param d  : pointer to the computed location in (x,y,z)-space
 */
int coons_patch(struct GPathPolyLine *c1,
		struct GPathPolyLine *c2,
		struct GPathPolyLine *c3,
		struct GPathPolyLine *c4,
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
    gpath_polyline_eval(c1, 0.0, &c1_0);
    gpath_polyline_eval(c1, 1.0, &c1_1);
    gpath_polyline_eval(c2, 0.0, &c2_0);
    gpath_polyline_eval(c2, 1.0, &c2_1);
    gpath_polyline_eval(c1, r, &c1_r);
    gpath_polyline_eval(c2, r, &c2_r);
    gpath_polyline_eval(c3, s, &c3_s);
    gpath_polyline_eval(c4, s, &c4_s);
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
} /* end coons_patch() */


/** \brief Compute a location within the 3D volume at parameter position (r,s,t).
 *
 * This function allows us to specify the edges individually but
 * delegates the actual work to TFI_3D. 
 * The point of this function is to allow easy access from Python code.
 */
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
		    struct point_3D *d)
{
    struct GPathPolyLine *c[12];
    c[0] = c0; c[1] = c1; c[2] = c2; c[3] = c3; 
    c[4] = c4; c[5] = c5; c[6] = c6; c[7] = c7;
    c[8] = c8; c[9] = c9; c[10] = c10; c[11] = c11;
    return TFI_3D(c, r, s, t, d);
}


/** \brief Compute a location within the 3D volume at parameter position (r,s,t).
 *
 * Returns 0 for a normal return, -1 otherwise.
 *
 * \param c[] : array of pointers to the bounding curves
 *	These are normalized GPathPolyLines.
 *	See 3D CFD workbook page 27pp 2004.
 * \param r  : West to East parameter, 0 <= r <= 1.
 * \param s  : South to North parameter, 0 <= s <= 1.
 * \param t  : Bottom to Top parameter, 0 <= t <= 1.
 * \param d  : pointer to the computed location in (x,y,z)-space
 */
int TFI_3D(struct GPathPolyLine *c[],
	   double r, double s, double t,
	   struct point_3D *d)
{
    struct point_3D c01r, c32r, c45r, c76r;
    struct point_3D c03s, c12s, c56s, c47s;
    struct point_3D c04t, c15t, c26t, c37t;
    struct point_3D p000, p001, p010, p011, p100, p101, p110, p111;
    double BigC;
    /* <<<<<<<< need to put error checking in >>>>>>>> */
    if (r < 0.0) r = 0.0;
    if (r > 1.0) r = 1.0;
    if (s < 0.0) s = 0.0;
    if (s > 1.0) s = 1.0;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    /* Corner points. Could check consistency here. */
    gpath_polyline_eval(c[0], 0.0, &p000);
    gpath_polyline_eval(c[0], 1.0, &p100);
    gpath_polyline_eval(c[2], 0.0, &p010);
    gpath_polyline_eval(c[2], 1.0, &p110);

    gpath_polyline_eval(c[4], 0.0, &p001);
    gpath_polyline_eval(c[4], 1.0, &p101);
    gpath_polyline_eval(c[6], 0.0, &p011);
    gpath_polyline_eval(c[6], 1.0, &p111);

    /* Points along edges. */
    gpath_polyline_eval(c[0], r, &c01r);
    gpath_polyline_eval(c[2], r, &c32r);
    gpath_polyline_eval(c[1], s, &c12s);
    gpath_polyline_eval(c[3], s, &c03s);

    gpath_polyline_eval(c[4], r, &c45r);
    gpath_polyline_eval(c[6], r, &c76r);
    gpath_polyline_eval(c[5], s, &c56s);
    gpath_polyline_eval(c[7], s, &c47s);

    gpath_polyline_eval(c[8],  t, &c04t);
    gpath_polyline_eval(c[9],  t, &c15t);
    gpath_polyline_eval(c[10], t, &c26t);
    gpath_polyline_eval(c[11], t, &c37t);

    /* Now interpolate to get the location of the new point.  */
    BigC = (1-r) * (1-s) * (1-t) * p000.x + (1-r) * (1-s) * t * p001.x +
	   (1-r) * s * (1-t) * p010.x     + (1-r) * s * t * p011.x +
	   r * (1-s) * (1-t) * p100.x     + r * (1-s) * t * p101.x +
	   r * s * (1-t) * p110.x         + r * s * t * p111.x;
    d->x = 0.5 * ((1-s) * (1-t) * c01r.x + s * (1-t) * c32r.x + 
	   (1-r) * (1-t) * c03s.x + r * (1-t) * c12s.x +
	   (1-s) * t * c45r.x     + s * t * c76r.x + 
	   (1-r) * t * c47s.x     + r * t * c56s.x +
	   (1-r) * (1-t) * c03s.x + (1-r) * t * c47s.x + 
	   (1-r) * (1-s) * c04t.x + (1-r) * s * c37t.x +
	   r * (1-t) * c12s.x     + r * t * c56s.x + 
	   r * (1-s) * c15t.x     + r * s * c26t.x +
	   (1-s) * (1-t) * c01r.x + (1-s) * t * c45r.x	+
	   (1-r) * (1-s) * c04t.x + r * (1-s) * c15t.x +
	   s * (1-t) * c32r.x     + s * t * c76r.x + 
	   (1-r) * s * c37t.x     + r * s * c26t.x) - 2.0 * BigC;

    BigC = (1-r) * (1-s) * (1-t) * p000.y + (1-r) * (1-s) * t * p001.y +
	   (1-r) * s * (1-t) * p010.y     + (1-r) * s * t * p011.y +
	   r * (1-s) * (1-t) * p100.y     + r * (1-s) * t * p101.y +
	   r * s * (1-t) * p110.y         + r * s * t * p111.y;
    d->y = 0.5 * ((1-s) * (1-t) * c01r.y + s * (1-t) * c32r.y + 
	   (1-r) * (1-t) * c03s.y + r * (1-t) * c12s.y +
	   (1-s) * t * c45r.y     + s * t * c76r.y + 
	   (1-r) * t * c47s.y     + r * t * c56s.y +
	   (1-r) * (1-t) * c03s.y + (1-r) * t * c47s.y + 
	   (1-r) * (1-s) * c04t.y + (1-r) * s * c37t.y +
	   r * (1-t) * c12s.y     + r * t * c56s.y + 
	   r * (1-s) * c15t.y     + r * s * c26t.y +
	   (1-s) * (1-t) * c01r.y + (1-s) * t * c45r.y	+
	   (1-r) * (1-s) * c04t.y + r * (1-s) * c15t.y +
	   s * (1-t) * c32r.y     + s * t * c76r.y + 
	   (1-r) * s * c37t.y     + r * s * c26t.y) - 2.0 * BigC;

    BigC = (1-r) * (1-s) * (1-t) * p000.z + (1-r) * (1-s) * t * p001.z +
	   (1-r) * s * (1-t) * p010.z     + (1-r) * s * t * p011.z +
	   r * (1-s) * (1-t) * p100.z     + r * (1-s) * t * p101.z +
	   r * s * (1-t) * p110.z         + r * s * t * p111.z;
    d->z = 0.5 * ((1-s) * (1-t) * c01r.z + s * (1-t) * c32r.z + 
	   (1-r) * (1-t) * c03s.z + r * (1-t) * c12s.z +
	   (1-s) * t * c45r.z     + s * t * c76r.z + 
	   (1-r) * t * c47s.z     + r * t * c56s.z +
	   (1-r) * (1-t) * c03s.z + (1-r) * t * c47s.z + 
	   (1-r) * (1-s) * c04t.z + (1-r) * s * c37t.z +
	   r * (1-t) * c12s.z     + r * t * c56s.z + 
	   r * (1-s) * c15t.z     + r * s * c26t.z +
	   (1-s) * (1-t) * c01r.z + (1-s) * t * c45r.z	+
	   (1-r) * (1-s) * c04t.z + r * (1-s) * c15t.z +
	   s * (1-t) * c32r.z     + s * t * c76r.z + 
	   (1-r) * s * c37t.z     + r * s * c26t.z) - 2.0 * BigC;

    return 0;
} /* end TFI_3D() */

/*-----------------------------------------------------------------*/

/** \brief Translate the straight-line position by dx,dy,dz.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param B   : (pointer to) the array of 2 control points
 * \param dx, dy, dz
 */
int line_translate(struct point_3D B[], double dx, double dy, double dz)
{
    if ( B == NULL ) return -1;
    point_3D_translate( &(B[0]), dx, dy, dz );
    point_3D_translate( &(B[1]), dx, dy, dz );
    return 0;
}


/** \brief Evaluate the straight-line position with parameter t.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param B   : (pointer to) the array of 2 control points
 * \param t   : position parameter, 0 <= t <= 1.0
 * \param loc : (pointer to) the computed location in 3D
 */
int line_eval(struct point_3D B[], double t, struct point_3D *loc)
{
    struct point_3D a, b;
    if ( B == NULL ) return -1;
    vector_copy_3D( &(B[0]), &a );
    vector_copy_3D( &(B[1]), &b );
    vector_scale_3D( &a, (1.0 - t) );
    vector_scale_3D( &b, t );
    vector_sum_3D( &a, &b, loc );
    return 0;
}


/** \brief Returns an estimate of the length of a straight line. 
 *
 * \param B : array of 2 control points
 */
double line_length(struct point_3D B[])
{
    struct point_3D ds;
    vector_diff_3D( &(B[0]), &(B[1]), &ds );
    return magnitude_3D( &ds );
}

/*-----------------------------------------------------------------*/

/** \brief Translate the arc position by dx,dy,dz.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param B   : (pointer to) the array of 3 control points
 * \param dx, dy, dz
 */
int arc_translate(struct point_3D B[], double dx, double dy, double dz)
{
    if ( B == NULL ) return -1;
    point_3D_translate( &(B[0]), dx, dy, dz );
    point_3D_translate( &(B[1]), dx, dy, dz );
    point_3D_translate( &(B[2]), dx, dy, dz );
    return 0;
}


/** \brief Evaluate the position on an arc for parameter t and the length 
 *         of the original arc. (This function does the real work.)
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 *
 * \param B   : (pointer to) the array of 3 control points
 * \param t   : position parameter, 0 <= t <= 1.0
 * \param loc : (pointer to) the computed location in 3D
 * \param length : (pointer to) the computed arc length
 */
int arc_eval_position_and_length(struct point_3D B[], 
				 double t, 
				 struct point_3D *loc, 
				 double *length)
{
    struct point_3D a, b, c;
    struct point_3D n, t1, t2;
    struct point_3D ca, cb, cb_local;
    double ca_mag, cb_mag, theta;

    if ( B == NULL ) {
	printf( "arc_eval(): NULL array of points\n" );
	return -1;
    }
    vector_copy_3D( &(B[0]), &a );   /* starting point */
    vector_copy_3D( &(B[1]), &b );   /* end point */
    vector_copy_3D( &(B[2]), &c );   /* centre of curvature */

    vector_diff_3D( &a, &c, &ca );
    vector_diff_3D( &b, &c, &cb );
    ca_mag = magnitude_3D( &ca );
    cb_mag = magnitude_3D( &cb );
    if ( fabs(ca_mag - cb_mag) / (1.0 + ca_mag) > 1.0e-6 ) {
	printf( "arc_eval(): radii don't match" );
	return -1;
    }
#   if 0
    printf( "arc_eval(): a=(%g %g %g), b=(%g %g %g)\n",
	    a.x, a.y, a.z, b.x, b.y, b.z );
    printf( "          : c=(%g %g %g), t=%g\n",
	    c.x, c.y, c.z, t );
#   endif

    /* First vector in plane */
    vector_copy_3D( &ca, &t1 );
    normalize_3D ( &t1 );
    /* Compute unit normal to plane of all three points. */
    vector_prod_3D( &ca, &cb, &n );
    if ( magnitude_3D( &n ) > 0.0 ) {
	normalize_3D( &n );
    } else {
	printf( "arc_eval(): cannot find plane of three points.\n" );
	return -1;
    }
    /* Third (orthogonal) vector is in the original plane. */
    vector_prod_3D( &n, &t1, &t2 );
#   if 0
    printf( "arc_eval(): n=(%g %g %g), t1=(%g %g %g), t2=(%g %g %g)\n",
	    n.x, n.y, n.z, t1.x, t1.y, t1.z, t2.x, t2.y, t2.z );
#   endif

    /* Now transform to local coordinates so that we can do 
     * the calculation of the point along the arc in 
     * the local xy-plane, with ca along the x-axis. */
    cb_local.x = cb.x * t1.x + cb.y * t1.y + cb.z * t1.z;
    cb_local.y = cb.x * t2.x + cb.y * t2.y + cb.z * t2.z;
    cb_local.z = cb.x * n.x  + cb.y * n.y  + cb.z * n.z;
#   if 0
    printf( "arc_eval(): cb_local=(%g %g %g)\n", 
	    cb_local.x, cb_local.y, cb_local.z );
#   endif
    if ( fabs(cb_local.z) > 1.0e-6 ) {
	printf( "arc_eval(): problems with transformation cb_local=(%g %g %g)\n",
		cb_local.x, cb_local.y, cb_local.z );
	return -1;
    }
    /* Angle of the final point on the arc is in the range -pi < th <= +pi */
    theta = atan2(cb_local.y, cb_local.x);

    /* The length of the circular arc. */
    *length = theta * cb_mag;

    /* Move the second point around the arc in the local xy-plane. */
    theta *= t;
    cb_local.x = cos(theta) * cb_mag;
    cb_local.y = sin(theta) * cb_mag;

    /* Transform back to global xyz coordinates
     * and remember to add the centre coordinates. */
    loc->x = cb_local.x * t1.x + cb_local.y * t2.x + cb_local.z * n.x + c.x;
    loc->y = cb_local.x * t1.y + cb_local.y * t2.y + cb_local.z * n.y + c.y;
    loc->z = cb_local.x * t1.z + cb_local.y * t2.z + cb_local.z * n.z + c.z;
#   if 0
    printf( "arc_eval(): global coords loc=(%g %g %g)\n",
	    loc->x, loc->y, loc->z );
#   endif
    return 0;
}


/** \brief Evaluate the position on an arc for parameter t.
 *
 * Returns 0 if everything was OK, -1 if there was a problem.
 * 
 * The real work is delegated to a helper function because there
 * is so much common code between the arc_length calculation and
 * the position calculation.
 *
 * \param B   : (pointer to) the array of 3 control points
 * \param t   : position parameter, 0 <= t <= 1.0
 * \param loc : (pointer to) the computed location in 3D
 */
int arc_eval(struct point_3D B[], double t, struct point_3D *loc)
{
    double length;
    return arc_eval_position_and_length(B, t, loc, &length);
}

/** \brief Returns an estimate of the length of an arc. 
 *
 * The real work is delegated.
 *
 * \param B : array of 3 control points
 */
double arc_length(struct point_3D B[])
{
    double length;
    struct point_3D loc;
    if (arc_eval_position_and_length(B, 1.0, &loc, &length) != 0) {
	printf( "arc_length(): problems with computing arc length\n" );
	return 0.0;
    }
    return length;
}

/*-----------------------------------------------------------------*/
