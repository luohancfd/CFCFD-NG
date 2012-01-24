/** \file  geom.h
 * \ingroup geom
 * \brief Vector Geometry Routines
 *
 * \author PA Jacobs and the CFCFD members
 *
 * \version 1.0    : 25-May-93 : extracted from ssp.h
 * \version 1.1      23-Jun-95 : added quad_area and quad_cent.  IAJ
 * \version 1.2      09-Jun-96 : merged with fixes from original
 *
 */

#ifndef GEOM_H

/*
 * ---------------
 * Data structures
 * ---------------
 */

#define POINT_3D_LABEL_SIZE 32

/** \brief A point in 2-dimensional vector space. */
struct Vector2D 
{ 
    double x;            /**< Cartesian coordinate x   */
    double y;            /**< Cartesian coordinate y   */
    unsigned char code;  /**< used for clipping in ssp */
};

/** \brief A point in 3-dimensional vector space. */
struct Vector3D 
{ 
    double x;           /**< Cartesian coordinate x   */
    double y;           /**< Cartesian coordinate y   */
    double z;           /**< Cartesian coordinate z   */
    char label[POINT_3D_LABEL_SIZE]; /**< text label  */
    unsigned char code; /**< used for clipping in ssp */
};

/* The following are for historical reasons... */
/** \brief synonom for Vector2D */
#define vector_2D Vector2D 
/** \brief synonom for Vector2D */
#define point_2D  Vector2D
/** \brief synonom for Vector3D */
#define vector_3D Vector3D
/** \brief synonom for Vector3D */
#define point_3D  Vector3D


/*
 * -------------------
 * Function Prototypes...
 * -------------------
 */
struct Vector3D *create_point_3D(double x, double y, double z);
void free_point_3D(struct Vector3D *p);
char * point_3D_set_label(struct Vector3D *p, const char *label);
char * point_3D_get_label(struct Vector3D *p);
struct Vector3D *create_point_3D_array(int np);
struct Vector3D *get_point_3D_ptr(struct Vector3D *pa, int i);


int vector_copy_3D( struct Vector3D *a, struct Vector3D *b );

double magnitude_3D( struct Vector3D *a );

int normalize_3D( struct Vector3D *a );

int vector_scale_3D( struct Vector3D *a, double scale );

double scalar_prod_3D( struct Vector3D *a, struct Vector3D *b );

int vector_prod_3D( struct Vector3D *a,
                     struct Vector3D *b,
                     struct Vector3D *c );

double triple_prod_3D( struct Vector3D *a,
                     struct Vector3D *b,
                     struct Vector3D *c );

int vector_sum_3D( struct Vector3D *a,
                     struct Vector3D *b,
                     struct Vector3D *c );

int point_3D_translate( struct Vector3D *a,
			double dx, double dy, double dz );

int vector_diff_3D( struct Vector3D *a,
                     struct Vector3D *b,
                     struct Vector3D *c );

int quad_properties_3D(struct Vector3D *p0,
		       struct Vector3D *p1,
		       struct Vector3D *p2,
		       struct Vector3D *p3,
		       struct Vector3D *centroid,
		       struct Vector3D *n,
		       struct Vector3D *t1,
		       struct Vector3D *t2,
		       double *area );

double tetrahedron_properties_3D(struct Vector3D *p0,
				 struct Vector3D *p1,
				 struct Vector3D *p2,
				 struct Vector3D *p3,
				 struct Vector3D *centroid,
				 double *volume );

double wedge_properties_3D(struct Vector3D *p0,
			   struct Vector3D *p1,
			   struct Vector3D *p2,
			   struct Vector3D *p3,
			   struct Vector3D *p4,
			   struct Vector3D *p5,
			   struct Vector3D *centroid,
			   double *volume );

double hexahedron_properties_3D(struct Vector3D *p0,
				struct Vector3D *p1,
				struct Vector3D *p2,
				struct Vector3D *p3,
				struct Vector3D *p4,
				struct Vector3D *p5,
				struct Vector3D *p6,
				struct Vector3D *p7,
				struct Vector3D *centroid,
				double *volume );

double quad_area(double xA, double xB, double xC, double xD,
                 double yA, double yB, double yC, double yD);

struct Vector3D quad_cent(double xA, double xB, double xC, double xD,
                          double yA, double yB, double yC, double yD);

int local_frame_3D(struct Vector3D *v, struct Vector3D *n, 
		   struct Vector3D *t1, struct Vector3D *t2);

int xyz_frame_3D(struct Vector3D *v, struct Vector3D *n, 
		 struct Vector3D *t1, struct Vector3D *t2);

/** A tag so that we know that this header file has been included. */
#define  GEOM_H  1
#endif
