/** \file geom.hh
 *  \ingroup libgeom2
 *  \brief Declarations for the C++ Vector3 class.
 *  \author PJ
 *  \version 27-Dec-2005
 *  \version 31-Dec-2005 eliminated label from Vector3 but kept it in Node3.
 *                       disposed of Node3::copy()
 *  \version 20-Jan-2006 Doxygen comments.
 *
 */

#ifndef GEOM_HH
#define GEOM_HH

#include <string>
#include <iostream>
#include <vector>
using namespace std;

/** \brief Basic 3D vectors.
 *
 * These are used as the foundation of the geometric paths, surfaces and volumes
 * used to generate grids for the CFCFD simulation codes.
 */
class Vector3 {
public:
    /// The Cartesian coordinates are accessible.
    double x;
    double y;
    double z;
    /// Constructor accept Cartesian coordinates.
    Vector3( double x=0.0, double y=0.0, double z=0.0 );
    /// Copy constructor.
    Vector3( const Vector3 &v );
    virtual ~Vector3();
    // Clone function (when dealing with pointers to Vector3s)
    Vector3* clone() const;
    /// Returns a string representation like Python __str__()
    virtual string str(int precision=6) const;
    /// Returns Virtual Reality Modelling Language string representation.
    virtual string vrml_str( double radius ) const;
    /// Returns VTK string representation.
    virtual string vtk_str() const;

    // Try inlining the following few transformation functions.
    // They will be used frequently in the flow simulation code Elmer3.
    // Change of coordinate system; rotation with translation.

    /// Transform coordinates from global frame to local (dash) frame.
    /// Local frame is defined by unit vectors (xdash, ydash and zdash) at location c.
    inline Vector3& transform_to_local( const Vector3 &xdash, const Vector3 &ydash, 
					const Vector3 &zdash, const Vector3 &c )
    {
	x -= c.x; y -= c.y; z -= c.z; // shift to local origin
	double new_x = x * xdash.x + y * xdash.y + z * xdash.z;
	double new_y = x * ydash.x + y * ydash.y + z * ydash.z;
	double new_z = x * zdash.x + y * zdash.y + z * zdash.z;
	x = new_x; y = new_y; z = new_z;
	return *this;
    }
    /// Transform coordinates from local (dash) frame to global frame.
    /// Local frame is defined by unit vectors (xdash, ydash and zdash) at location c.
    inline Vector3& transform_to_global( const Vector3 &xdash, const Vector3 &ydash, 
					 const Vector3 &zdash, const Vector3 &c )
    {
	double new_x = x * xdash.x + y * ydash.x + z * zdash.x + c.x;
	double new_y = x * xdash.y + y * ydash.y + z * zdash.y + c.y;
	double new_z = x * xdash.z + y * ydash.z + z * zdash.z + c.z;
	x = new_x; y = new_y; z = new_z;
	return *this;
    }
    // Change of coordinate system; just rotation.
    inline Vector3& transform_to_local( const Vector3 &xdash, const Vector3 &ydash, 
					const Vector3 &zdash )
    {
	double new_x = x * xdash.x + y * xdash.y + z * xdash.z;
	double new_y = x * ydash.x + y * ydash.y + z * ydash.z;
	double new_z = x * zdash.x + y * zdash.y + z * zdash.z;
	x = new_x; y = new_y; z = new_z;
	return *this;
    }
    inline Vector3& transform_to_global( const Vector3 &xdash, const Vector3 &ydash, 
					 const Vector3 &zdash )
    {
	double new_x = x * xdash.x + y * ydash.x + z * zdash.x;
	double new_y = x * xdash.y + y * ydash.y + z * zdash.y;
	double new_z = x * xdash.z + y * ydash.z + z * zdash.z;
	x = new_x; y = new_y; z = new_z;
	return *this;
    }

    /// Mirror image in the plane defined by a point and a normal vector.
    Vector3& mirror_image( const Vector3 &point, const Vector3 &normal );
    /// Rotate the point about the z-axis by the specified angle (in radians).
    Vector3& rotate_about_zaxis( double dtheta );
    // arithmetic operators follow; more are defined as helper functions, later.
    // We rely upon implicit conversions to get double and int operators.
    // Assignment operator
    Vector3& operator=( const Vector3 &v );
    /// Add v to this vector.
    Vector3& operator+=( const Vector3 &v );
    /// Subtract v from this vector.
    Vector3& operator-=( const Vector3 &v );
    /// Scale by v (real number).
    Vector3& operator*=( double v );
    /// Scale by 1/v.
    Vector3& operator/=( double v );
    /// Normalize the vector (so that it is a unit vector.
    Vector3& norm();
};

// Helper functions.
/// Output a string representation of the vector.
ostream& operator<<( ostream &os, const Vector3 &v );
/// Returns the magnitude of the vector.
double vabs( const Vector3 &v );
/// Returns a (new) unit vector with the same direction as v.
Vector3 unit( const Vector3 &v );
/// Unary +
Vector3 operator+( const Vector3 &v );
/// Unary -
Vector3 operator-( const Vector3 &v );
/// Binary + returns new vector (v1 + v2).
Vector3 operator+( const Vector3 &v1, const Vector3 &v2 );
/// Binary - returns new vector (v1 - v2).
Vector3 operator-( const Vector3 &v1, const Vector3 &v2 );
/// Returns True if v1 and v2 are within tolerance of each other.
bool equal( const Vector3 &v1, const Vector3 &v2, double tolerance=1.0e-12);
/// Binary * returns a (new) vector (v1 * v2).
Vector3 operator*( const Vector3 &v1, double v2 );
/// Binary * returns a (new) vector (v1 * v2).
Vector3 operator*( double v1, const Vector3 &v2 );
/// Binary / returns a (new) vector (v1 / v2).
Vector3 operator/( const Vector3 &v1, double v2 );
/// Vector dot product returns a real number.
double dot( const Vector3 &v1, const Vector3 &v2 );
/// Vector cross product returns a (new) vector value.
Vector3 cross( const Vector3 &v1, const Vector3 &v2 );

//--------------------------------------------------------------------------
// Geometry functions that use Vector3 objects.

/// Projects q along qr onto the plane defined by the triangle abc.
/// Returns 0 for successful projection, 1 if qr is parallel and q is on plane,
/// 2 if qr is parallel and q is not on plane.
int project_onto_plane( Vector3 &q, const Vector3 &qr,
			const Vector3 &a, const Vector3 &b, const Vector3 &c );
/// Returns 1 if p in inside triangle abc or on its boundary, 0 otherwise.
int inside_triangle( const Vector3 &p, 
		     const Vector3 &a, const Vector3 &b, const Vector3 &c);

/// For Hannes and Paul's turbine-blade grids,
/// map space so that a neutral plane wraps onto a cylinder of radius H.
int map_neutral_plane_to_cylinder( Vector3 &p, double H );

// -------------------------------------------------------------------- 

/** \brief Computes the geometric properties of the quadrilateral p0123.
 *
 * \param p0, p1, p2, p3 : IN : vertices of quadrilateral as shown below
 * \param centroid : OUT : reference to the centroidal position
 * \param n    : OUT : reference to the unit normal
 * \param t1   : OUT : reference to unit-tangent 1 (parallel to p0--p1)
 * \param t2   : OUT : reference to unit-tangent-2
 * \param area : OUT : reference to a double to store the computed area.
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
int quad_properties( const Vector3 &p0, const Vector3 &p1, 
		     const Vector3 &p2, const Vector3 &p3,
		     Vector3 &centroid,
		     Vector3 &n, Vector3 &t1, Vector3 &t2,
		     double &area );

Vector3 quad_centroid( const Vector3 &p0, const Vector3 &p1, 
		       const Vector3 &p2, const Vector3 &p3 );

Vector3 quad_normal( const Vector3 &p0, const Vector3 &p1, 
		     const Vector3 &p2, const Vector3 &p3 );

Vector3 quad_tangent1( const Vector3 &p0, const Vector3 &p1, 
		       const Vector3 &p2, const Vector3 &p3 );

Vector3 quad_tangent2( const Vector3 &p0, const Vector3 &p1, 
		       const Vector3 &p2, const Vector3 &p3 );

double quad_area( const Vector3 &p0, const Vector3 &p1, 
		  const Vector3 &p2, const Vector3 &p3 );

/** \brief Computes the geometric properties of the tetrahedron p0123.
 *
 * \param p0, p1, p2, p3 : IN : vertices of tetrahedron as shown below
 * \param centroid : OUT : reference to the centroidal position
 * \param volume   : OUT : reference to a double to store the computed volume.
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
int tetrahedron_properties( const Vector3 &p0, const Vector3 &p1,
			    const Vector3 &p2, const Vector3 &p3,
			    Vector3 &centroid, double &volume );
Vector3 tetrahedron_centroid( const Vector3 &p0, const Vector3 &p1,
			      const Vector3 &p2, const Vector3 &p3 );
double tetrahedron_volume( const Vector3 &p0, const Vector3 &p1,
			   const Vector3 &p2, const Vector3 &p3 );

/** \brief Computes the geometric properties of the wedge p012345.
 *
 * Divide the triangular-based wedge into 3 tetrahedra and combine
 * the properties for the tetrahedra.
 *
 * Returns 0 if there were no problems, -1 if volume computes as negative.
 *
 * \param p0, p1, p2, p3, p4, p5
 *                    IN : vertices of wedge as shown below
 * \param centroid : OUT : reference to the centroidal position
 * \param volume   : OUT : reference to a double to store the computed volume.
 *
 * \verbatim
 * Base of wedge: p012
 * Top of wedge: p345
 * \endverbatim
 */
int wedge_properties( const Vector3 &p0, const Vector3 &p1,
		      const Vector3 &p2, const Vector3 &p3,
		      const Vector3 &p4, const Vector3 &p5,
		      Vector3 &centroid, double &volume );
Vector3 wedge_centroid( const Vector3 &p0, const Vector3 &p1,
			const Vector3 &p2, const Vector3 &p3,
			const Vector3 &p4, const Vector3 &p5);
double wedge_volume( const Vector3 &p0, const Vector3 &p1,
		     const Vector3 &p2, const Vector3 &p3,
		     const Vector3 &p4, const Vector3 &p5);

/** \brief Computes the geometric properties of the hexahedron p01234567.
 *
 * \param p0, p1, p2, p3, p4, p5, p6, p7 : 
 *                    IN : vertices of hexahedron as shown below
 * \param centroid : OUT : reference to the centroidal position
 * \param volume   : OUT : reference to a double to store the computed volume.
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
int hexahedron_properties( const Vector3 &p0, const Vector3 &p1,
			   const Vector3 &p2, const Vector3 &p3,
			   const Vector3 &p4, const Vector3 &p5,
			   const Vector3 &p6, const Vector3 &p7,
			   Vector3 &centroid, double &volume );

double hexahedron_volume( const Vector3 &p0, const Vector3 &p1,
			  const Vector3 &p2, const Vector3 &p3,
			  const Vector3 &p4, const Vector3 &p5,
			  const Vector3 &p6, const Vector3 &p7 );

Vector3 hexahedron_centroid( const Vector3 &p0, const Vector3 &p1,
			     const Vector3 &p2, const Vector3 &p3,
			     const Vector3 &p4, const Vector3 &p5,
			     const Vector3 &p6, const Vector3 &p7 );

/** \brief Change from a global Cartesian frame of reference to a frame
 *        with its x-direction normal to the interface.
 *
 * The Riemann solver works in this local frame of reference.
 *
 * \param v         : reference to a Vector3D in xyz-frame
 *                    The vector is written over.
 * \param n, t1, t2 : references to the unit vectors for the interface.
 */
int local_frame( Vector3 &v, const Vector3 &n, 
		 const Vector3 &t1, const Vector3 &t2 );

/** \brief Change from a frame of reference 
 *        with its x-direction normal to the interface.
 *        to a global Cartesian frame of reference. 
 *
 * \param v         : reference to a Vector3D in local frame
 *                    The vector is written over.
 * \param n, t1, t2 : references to the unit vectors for the interface.
 */
int xyz_frame( Vector3 &v, const Vector3 &n, 
	       const Vector3 &t1, const Vector3 &t2 );

#endif
