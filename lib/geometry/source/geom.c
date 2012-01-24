/** \file  geom.c
 * \ingroup geom
 * \brief Vector Geometry Routines
 *
 * \author PA Jacobs and the CFCFD members 
 *
 * \version 1.0  : 25-May-93 : first cut, 3D routines only
 * \version 1.01 : 31-May-93 : added triple product
 * \version 1.02 : 01-Jun-93 : added normalize vector
 * \version 1.03 : 23-Jun-93 : includes quadrilateral area and centroid routines
 * \version 1.04 : 16-May-04 : as pointed out by Dave Towsey, the small
 *                             magnitude traps may ne be such a good idea.
 * \version 1.05 : 03-Jul-04 : Major tidy-up and addition of doxygen tags.
 * \version 1.06 : 23-Aug-04 : hex, wedge and tetradedron properties.
 *
 */

/*----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "geom.h"

/** A small number to test for effectively-zero values. */
#define VERY_SMALL_MAGNITUDE (1.0e-200)

/*----------------------------------------------------------------*/

/** \brief Returns a pointer to a single point_3D structure allocated on the heap.
 *
 * This will be useful for interfacing to Python and the like.
 */
struct Vector3D *create_point_3D(double x, double y, double z)
{
    struct Vector3D *p;
    p = (struct Vector3D*) calloc(sizeof(struct Vector3D), 1);
    if (p != NULL) {
	p->x = x;
	p->y = y;
	p->z = z;
    }
    return p;
}

char * point_3D_set_label(struct Vector3D *p, const char *label)
{
    if ( p == NULL ) return NULL;
    strncpy(p->label, label, POINT_3D_LABEL_SIZE-1);
    p->label[POINT_3D_LABEL_SIZE-1] = '\0';
    return p->label;
}

char * point_3D_get_label(struct Vector3D *p)
{
    if ( p == NULL ) return NULL;
    return p->label;
}

/** \brief Returns a pointer to an array of point_3D structures allocated on the heap.
 *
 * This will be useful for interfacing to Python and the like.
 */
struct Vector3D *create_point_3D_array(int np)
{
    struct Vector3D *pa;
    int i;
    pa = (struct Vector3D*) calloc(sizeof(struct Vector3D), np);
    if (pa != NULL) {
	for ( i = 0; i < np; ++i ) {
	    pa[i].x = 0.0;
	    pa[i].y = 0.0;
	    pa[i].z = 0.0;
	}
    }
    return pa;
}


/** \brief Gets a pointer to the i-th point_3D structure
 *         that is in an array pa of such structures.
 *
 * This will be useful for interfacing to Python and the like.
 */
struct Vector3D *get_point_3D_ptr(struct Vector3D *pa, int i)
{
    if (pa != NULL) {
	return &(pa[i]);
    } else {
	return NULL;
    }
}


/** \brief Frees the memory associated with a point_3D structure.
 *
 * This will be useful for interfacing to Python and the like.
 */
void free_point_3D(struct Vector3D *p)
{
    if (p != NULL) free(p);
}


/** \brief Copies the contents of vector a into vector b.
 * \param *a  : pointer to vector a
 * \param *b  : pointer to vector b
 */
int vector_copy_3D( struct Vector3D *a, struct Vector3D *b )
{
    b->x = a->x;
    b->y = a->y;
    b->z = a->z;
    return 0;
}  /* end of vector_copy_3D() */


/** \brief Returns the magnitude of vector a.
 * \param *a  : pointer to vector a
 */
double magnitude_3D( struct Vector3D *a )
{
    double t, xx, yy, zz;
    xx = a->x;
    yy = a->y;
    zz = a->z;
    t = xx * xx + yy * yy + zz * zz;
    if (t > VERY_SMALL_MAGNITUDE)   
	t = sqrt ( t );    /* A reasonable finite length. */
    else
	t = 0.0;           /* Assume zero */
    return t;
}  /* end of magnitude_3D() */


/** \brief Replaces contents of vector a with (scale * a).
 * \param *a    : pointer to vector a
 * \param scale : the scale factor
 */
int vector_scale_3D( struct Vector3D *a, double scale )
{
    a->x = a->x * scale;
    a->y = a->y * scale;
    a->z = a->z * scale;
    return 0;
}  /* end of vector_scale_3D() */


/** \brief Normalizes the vector a (so that its magnitude is 1.0).
 * \param *a    : pointer to vector a
 */
int normalize_3D( struct Vector3D *a )
{
    double scale;
    scale = magnitude_3D (a);
    if (scale > 0.0) {
	/* Assume finite length */
	scale = 1.0 / scale;
	a->x = a->x * scale;
	a->y = a->y * scale;
	a->z = a->z * scale;
    }
    return 0;
}  /* end of normalize_3D() */


/** \brief Returns the scalar (dot) product a.b.
 * \param *a  : pointer to vector a
 * \param *b  : pointer to vector b
 */
double scalar_prod_3D( struct Vector3D *a, struct Vector3D *b )
{
    double t;
    t = b->x * a->x + b->y * a->y + b->z * a->z;
    return t;
}  /* end of scalar_prod_3D() */


/** \brief Computes the vector (cross) product axb.
 * \param *a  : pointer to vector a
 * \param *b  : pointer to vector b
 * \param *c  : the vector product is stored in vector c.
 */
int vector_prod_3D( struct Vector3D *a,
		    struct Vector3D *b,
		    struct Vector3D *c )
{
    c->x = a->y * b->z - a->z * b->y;
    c->y = a->z * b->x - a->x * b->z;
    c->z = a->x * b->y - a->y * b->x;
    return 0;
}  /* end of vector_prod_3D() */


/** \brief Returns the (scalar value of the) triple product a.(bxc).
 * \param *a  : pointer to vector a
 * \param *b  : pointer to vector b
 * \param *c  : pointer to vector c
 */
double triple_prod_3D( struct Vector3D *a,
		       struct Vector3D *b,
		       struct Vector3D *c )
{
    struct Vector3D t1;
    vector_prod_3D(b, c, &t1);
    return scalar_prod_3D(a, &t1);
}  /* end of triple_prod_3D() */


/** \brief Computes the vector sum c=a+b.
 * \param *a  : pointer to vector a
 * \param *b  : pointer to vector b
 * \param *c  : the result is stored in vector c.
 */
int vector_sum_3D( struct Vector3D *a,
		   struct Vector3D *b,
		   struct Vector3D *c )
{
    c->x = a->x + b->x;
    c->y = a->y + b->y;
    c->z = a->z + b->z;
    return 0;
}  /* end of vector_sum_3D() */


/** \brief Translates a point in Cartesian space.
 * \param *a  : pointer to vector (or point) a
 * \param dx, dy, dz : the offest to be applied
 */
int point_3D_translate( struct Vector3D *a,
			double dx, double dy, double dz )
{
    a->x += dx;
    a->y += dy;
    a->z += dz;
    return 0;
}  /* end of translate_3D_point() */


/** \brief Computes the vector difference c=a-b.
 * \param *a  : pointer to vector a
 * \param *b  : pointer to vector b
 * \param *c : the result is stored in vector c.
 *
 */
int vector_diff_3D( struct Vector3D *a,
		    struct Vector3D *b,
		    struct Vector3D *c )
{
    c->x = a->x - b->x;
    c->y = a->y - b->y;
    c->z = a->z - b->z;
    return 0;
}  /* end of vector_diff_3D() */


/** \brief Computes the geometric properties of the quadrilateral p0123.
 *
 * \param p0, p1, p2, p3 : IN : vertices of quadrilateral as shown below
 * \param centroid : OUT : pointer to the centroidal position
 * \param n    : OUT : pointer to the unit normal
 * \param t1   : OUT : pointer to unit-tangent 1 (parallel to p0--p1)
 * \param t2   : OUT : pointer to unit-tangent-2
 * \param area : OUT : pointer to a double to store the computed area.
 *
 * Of course, the vertices may not be coplanar.
 * It is all a bit of a fudge in that case.
 *
 * \verbatim
 *   3-----2
 *   |     |
 *   |     |
 *   0-----1
 * \endverbatim
 */
int quad_properties_3D(struct Vector3D *p0,
		       struct Vector3D *p1,
		       struct Vector3D *p2,
		       struct Vector3D *p3,
		       struct Vector3D *centroid,
		       struct Vector3D *n,
		       struct Vector3D *t1,
		       struct Vector3D *t2,
		       double *area )
{
    struct Vector3D s01, s12, s30, s32, s02, s13;
    struct Vector3D area032, area210;

    /* Work out direction vectors */
    vector_diff_3D( p1, p0, &s01 );
    vector_diff_3D( p2, p1, &s12 );
    vector_diff_3D( p0, p3, &s30 );
    vector_diff_3D( p2, p3, &s32 );

    /* Find the cross products */
    vector_prod_3D( &s30, &s32, &area032 );
    vector_prod_3D( &s01, &s12, &area210 );
    vector_scale_3D( &area032, 0.5 );
    vector_scale_3D( &area210, 0.5 );

    /* unit-normal and area */
    vector_sum_3D( &area210, &area032, n );
    *area = magnitude_3D( n );
    normalize_3D( n );

    /* Tangent unit-vectors: 
     * t1 is parallel to s01, 
     * t2 is normal to n and t1 */
    vector_copy_3D( &s01, t1 );
    normalize_3D( t1 );
    vector_prod_3D( n, t1, t2 );

    /* Centroid: average mid-points of the diagonals. */
    vector_sum_3D( p0, p2, &s02 );
    vector_sum_3D( p1, p3, &s13 );
    vector_sum_3D( &s02, &s13, centroid );
    vector_scale_3D( centroid, 0.25 );

    return 0; 
} /* end quad_properties_3D() */


/** \brief Computes the geometric properties of the tetrahedron p0123.
 *
 * \param p0, p1, p2, p3 : IN : vertices of tetrahedron as shown below
 * \param centroid : OUT : pointer to the centroidal position
 * \param volume   : OUT : pointer to a double to store the computed volume.
 *
 * \verbatim
 * Base of tetrahedron: p012.
 * Peak is at p3.
 *
 *         2
 *        /|
 *       / |
 *      /  |
 *     /   |
 *    /    |
 *   0-----1
 * \endverbatim
 */
double tetrahedron_properties_3D(struct Vector3D *p0,
				 struct Vector3D *p1,
				 struct Vector3D *p2,
				 struct Vector3D *p3,
				 struct Vector3D *centroid,
				 double *volume )
{
    struct Vector3D s01, s02, s03, s13;

    /* Work out direction vectors */
    vector_diff_3D( p1, p0, &s01 );
    vector_diff_3D( p2, p0, &s02 );
    vector_diff_3D( p3, p0, &s03 );

    /* Find the cross products */
    *volume = triple_prod_3D( &s03, &s01, &s02 ) / 6.0;

    /* Centroid: average the vertex locations. */
    vector_sum_3D( p0, p2, &s02 );
    vector_sum_3D( p1, p3, &s13 );
    vector_sum_3D( &s02, &s13, centroid );
    vector_scale_3D( centroid, 0.25 );
#   if 0
    printf("tetrahedron: p0=(%e,%e,%e), p1=(%e,%e,%e)\n",
	   p0->x, p0->y, p0->z, p1->x, p1->y, p1->z);
    printf("           : p2=(%e,%e,%e), p3=(%e,%e,%e)\n",
	   p2->x, p2->y, p2->z, p3->x, p3->y, p3->z);
    printf("    centre : pc=(%e,%e,%e), volume=%e\n",
	   centroid->x, centroid->y, centroid->z, *volume);
#   endif
    return 0; 
} /* end tetrahedron_properties_3D() */


/** \brief Computes the geometric properties of the wedge p012345.
 *
 * Divide the triangular-based wedge into 3 tetrahedra and combine
 * the properties for the tetrahedra.
 *
 * Returns 0 if there were no problems, -1 if volume computes as negative.
 *
 * \param p0, p1, p2, p3, p4, p5
 *                    IN : vertices of wedge as shown below
 * \param centroid : OUT : pointer to the centroidal position
 * \param volume   : OUT : pointer to a double to store the computed volume.
 *
 * \verbatim
 * Base of wedge: p012
 * Top of wedge: p345
 * \endverbatim
 */
double wedge_properties_3D(struct Vector3D *p0,
			   struct Vector3D *p1,
			   struct Vector3D *p2,
			   struct Vector3D *p3,
			   struct Vector3D *p4,
			   struct Vector3D *p5,
			   struct Vector3D *centroid,
			   double *volume )
{
    double v1, v2, v3, v_tot;
    struct Vector3D c1, c2, c3;

    tetrahedron_properties_3D( p0, p4, p5, p3, &c1, &v1 );
    tetrahedron_properties_3D( p0, p5, p4, p1, &c2, &v2 );
    tetrahedron_properties_3D( p0, p1, p2, p5, &c3, &v3 );

    v_tot = v1 + v2 + v3;
    if (v_tot < VERY_SMALL_MAGNITUDE) {
	printf("wedge_properties_3D(): zero or negative volume: %e\n", v_tot);
	*volume = 0.0;
	/* equally-weighted tetrahedral centroids. */
	vector_sum_3D( &c1, &c2, centroid );
	vector_sum_3D( &c3, centroid, centroid );
	vector_scale_3D( centroid, 1.0/3.0 );
	return -1;
    }

    /* Weight the tetrahedral centroids with their volumes. */
    vector_scale_3D( &c1, v1 );
    vector_scale_3D( &c2, v2 );
    vector_scale_3D( &c3, v3 );
    vector_sum_3D( &c1, &c2, centroid );
    vector_sum_3D( &c3, centroid, centroid );
    vector_scale_3D( centroid, 1.0/v_tot );
    *volume = v_tot;
#   if 0
    printf("wedge      : p0=(%e,%e,%e), p1=(%e,%e,%e)\n",
	   p0->x, p0->y, p0->z, p1->x, p1->y, p1->z);
    printf("           : p2=(%e,%e,%e), p3=(%e,%e,%e)\n",
	   p2->x, p2->y, p2->z, p3->x, p3->y, p3->z);
    printf("           : p4=(%e,%e,%e), p5=(%e,%e,%e)\n",
	   p4->x, p4->y, p4->z, p5->x, p5->y, p5->z);
    printf("    centre : pc=(%e,%e,%e), volume=%e\n",
	   centroid->x, centroid->y, centroid->z, *volume);
#   endif
    return 0;
} /* end wedge_properties_3D() */


/** \brief Computes the geometric properties of the hexahedron p01234567.
 *
 * \param p0, p1, p2, p3, p4, p5, p6, p7 : 
 *                    IN : vertices of hexahedron as shown below
 * \param centroid : OUT : pointer to the centroidal position
 * \param volume   : OUT : pointer to a double to store the computed volume.
 *
 * \verbatim
 * Base of hexahedron: p0123; this view looking down.
 *
 *       7-----6
 *      /|    /|  Top
 *     / |   / |
 *    /  4-----5
 *   3--/--2  /
 *   | /   | /
 *   |/    |/ Bottom
 *   0-----1
 * \endverbatim
 */
double hexahedron_properties_3D(struct Vector3D *p0,
				struct Vector3D *p1,
				struct Vector3D *p2,
				struct Vector3D *p3,
				struct Vector3D *p4,
				struct Vector3D *p5,
				struct Vector3D *p6,
				struct Vector3D *p7,
				struct Vector3D *centroid,
				double *volume )
{
    double v1, v2, v_tot;
    struct Vector3D c1, c2;

    wedge_properties_3D( p0, p1, p2, p4, p5, p6, &c1, &v1 );
    wedge_properties_3D( p0, p2, p3, p4, p6, p7, &c2, &v2 );

    v_tot = v1 + v2;
    if (v_tot < VERY_SMALL_MAGNITUDE) {
	printf("hexahedron_properties_3D(): zero or negative volume: %e\n", v_tot);
	*volume = 0.0;
	/* equally-weighted prism centroidal values. */
	vector_sum_3D( &c1, &c2, centroid );
	vector_scale_3D( centroid, 0.5 );
	return -1;
    }

    /* Volume-weight the prism centroidal values. */
    vector_scale_3D( &c1, v1 );
    vector_scale_3D( &c2, v2 );
    vector_sum_3D( &c1, &c2, centroid );
    vector_scale_3D( centroid, 1.0/v_tot );
    *volume = v_tot;
    return 0; 
} /* end hexahedron_properties_3D() */


/** \brief Returns the area of the quadrilateral ABCD
 * \param xA,xB,xC,xD : x-coords of vertices of quadrilateral ABCD
 * \param yA,yB,yC,yD : y-coords of vertices of quadrilateral ABCD
 *
 * \verbatim
 *   C-----B
 *   |     |
 *   |     |
 *   D-----A
 * \endverbatim
 */
double quad_area(double xA, double xB, double xC, double xD,
                 double yA, double yB, double yC, double yD)
{
   struct Vector3D a,b,c,d,ab,ac,ad,cross1,cross2;
   double          area,s1,s2;

   /* Convert points to position vectors */
   a.x = xA;  a.y = yA;  a.z = 0.0;
   b.x = xB;  b.y = yB;  b.z = 0.0;
   c.x = xC;  c.y = yC;  c.z = 0.0;
   d.x = xD;  d.y = yD;  d.z = 0.0;

   /* Work out direction vectors */
   vector_diff_3D(&b, &a, &ab);
   vector_diff_3D(&c, &a, &ac);
   vector_diff_3D(&d, &a, &ad);

   /* Find the cross products */
   vector_prod_3D(&ab, &ac, &cross1);
   vector_prod_3D(&ac, &ad, &cross2);

   /* Get the area */
   s1 = 0.5 * cross1.z;
   s2 = 0.5 * cross2.z;
   area=fabs(s1 + s2);

   return area; 
} /* end quad_area() */


/** \brief Returns the centroid (as a 3D vector) of qualrilateral ABCD.
 * \param xA,xB,xC,xD: x-coords of vertices of quadrilateral ABCD
 * \param yA,yB,yC,yD: y-coords of vertices of quadrilateral ABCD
 *
 * \verbatim
 *   C-----B
 *   |     |
 *   |     |
 *   D-----A
 * \endverbatim
 */
struct Vector3D quad_cent(double xA, double xB, double xC, double xD,
                          double yA, double yB, double yC, double yD)
{
    struct Vector3D a,b,c,d,ab,ac,ad,cross1,cross2;
    struct Vector3D cent1, cent2, comp1, comp2, centroid;
    double          signed_area,s1,s2;
    double   one_third = (1.0/3.0);

    /* Work out the area of each triangular component
     * of the quad as in quad_area
     */
    a.x = xA;  a.y = yA;  a.z = 0.0;
    b.x = xB;  b.y = yB;  b.z = 0.0;
    c.x = xC;  c.y = yC;  c.z = 0.0;
    d.x = xD;  d.y = yD;  d.z = 0.0;

    vector_diff_3D(&b, &a, &ab);
    vector_diff_3D(&c, &a, &ac);
    vector_diff_3D(&d, &a, &ad);

    vector_prod_3D(&ab, &ac, &cross1);
    vector_prod_3D(&ac, &ad, &cross2);

    s1 = 0.5 * cross1.z;
    s2 = 0.5 * cross2.z;
    /* Total area of the quad */
    signed_area = s1 + s2;

    /* Find the triangle centroids */
    vector_sum_3D(&a, &b, &cent1);
    vector_sum_3D(&c, &cent1, &cent1);
    vector_sum_3D(&a, &c, &cent2);
    vector_sum_3D(&d, &cent2, &cent2);
    vector_scale_3D(&cent1, one_third);
    vector_scale_3D(&cent2, one_third);

    /* Work out quad centroid. */
    comp1 = cent1;
    comp2 = cent2;
    vector_scale_3D(&comp1, (s1/signed_area));
    vector_scale_3D(&comp2, (s2/signed_area));
    vector_sum_3D(&comp1, &comp2, &centroid);
               
    return centroid;
} /* end quad_cent() */

/*----------------------------------------------------------------*/

/** \brief Change from a global (X,Y) frame of reference to a frame
 *        with its x-direction normal to the interface.
 *
 * The Riemann solver works in this local frame of reference.
 *
 * \param v         : pointer to a Vector3D in xyz-frame
 *                    The vector is written over.
 * \param n, t1, t2 : pointers to the unit vectors for the interface.
 */
int local_frame_3D(struct Vector3D *v, struct Vector3D *n, 
		   struct Vector3D *t1, struct Vector3D *t2)
{
    struct Vector3D v1;
    v1.x = v->x; v1.y = v->y; v1.z = v->z;
    v->x = scalar_prod_3D(&v1, n);     /* Normal velocity       */
    v->y = scalar_prod_3D(&v1, t1);    /* Tangential velocity 1 */
    v->z = scalar_prod_3D(&v1, t2);    /* Tangential velocity 2 */
    return 0;
}


/** \brief Change from a frame of reference 
 *        with its x-direction normal to the interface.
 *        to a global (X,Y) frame of reference. 
 *
 * \param v         : pointer to a Vector3D in local frame
 *                    The vector is written over.
 * \param n, t1, t2 : pointers to the unit vectors for the interface.
 */
int xyz_frame_3D(struct Vector3D *v, struct Vector3D *n, 
		 struct Vector3D *t1, struct Vector3D *t2)
{
    struct Vector3D v1;
    v1.x = v->x; v1.y = v->y; v1.z = v->z;
    v->x = v1.x * n->x + v1.y * t1->x + v1.z * t2->x; /* global-x component */
    v->y = v1.x * n->y + v1.y * t1->y + v1.z * t2->y; /* global-y component */
    v->z = v1.x * n->z + v1.y * t1->z + v1.z * t2->z; /* global-z component */
    return 0;
}

/* -------------- TEST driver --------------*/

#ifdef TEST_GEOM

int main()
{
    double sc;
    struct point_3D v1, v2, v3;

    printf("Begin test of geom.c...\n");

    v1.x = 2.0; v1.y = 2.0; v1.z = 0.0;
    v2.x = 0.0; v2.y = 0.0; v2.z = 1.0;

    vector_copy_3D(&v1, &v3);
    printf("v3 : (%f, %f, %f) \n", v3.x, v3.y, v3.z);

    sc = magnitude_3D( &v3 );
    printf("mag(v3) = %f \n", sc);

    vector_scale_3D( &v3, 1.0 / sc );
    printf("v3 : (%f, %f, %f) \n", v3.x, v3.y, v3.z);
    sc = magnitude_3D( &v3 );
    printf("mag(v3) = %f \n", sc);

    sc = scalar_prod_3D( &v1, &v2 );
    printf ("v1.v2 = %f \n", sc);

    vector_prod_3D( &v1, &v2, &v3 );
    printf ("v3 : (%f, %f, %f) \n", v3.x, v3.y, v3.z);

    vector_sum_3D( &v1, &v2, &v3 );
    printf ("v3 : (%f, %f, %f) \n", v3.x, v3.y, v3.z);

    vector_diff_3D( &v1, &v2, &v3 );
    printf("v3 : (%f, %f, %f) \n", v3.x, v3.y, v3.z);

    v1.x = 1.414; v1.y =  1.414; v1.z = 0.0;
    v2.x = 1.414; v2.y = -1.414; v2.z = 0.0;
    v3.x = 0.0;   v3.y =  0.0;   v3.z = 2.0;
    normalize_3D(&v1);
    normalize_3D(&v2);
    normalize_3D(&v3);
    printf("v1.(v2xv3) = %f\n", triple_prod_3D(&v1, &v2, &v3) );

    printf("Press <RETURN>...");
    getchar();
    printf("Done.\n");
    return 0;
}

#endif

/* -------------- end file ---------------- */


