/** \file geom.cxx
 *  \ingroup libgeom2
 *  \brief Implementation of the Vector3 class, 
 *         together with basic 3D geometry functions.
 *  \author PJ
 *  \version 27-Dec-2005 initial coding
 *  \version 31-Dec-2005 eliminated label from Vector3 but kept it in Node3.
 *                       disposed of Node3::copy()
 *  \version 17-Jan-2006 remove Node3 class
 *  \version 25-Jan-2006 Projection of a point onto the plane of the triangle.
 *
 * This module is a C++ replacement for the combined C+Python geom module that
 * was built for mb_cns and Elmer.
 */
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <math.h>
#include "geom.hh"
using namespace std;


/// A small number to test for effectively-zero values.
#define VERY_SMALL_MAGNITUDE (1.0e-200)
#define SMALL_BUT_SIGNIFICANT (1.0e-20)

// Constructors accept Cartesian coordinates or another vector.
Vector3::Vector3( double x, double y, double z )
    : x(x), y(y), z(z) {} 
Vector3::Vector3( const Vector3 &v )
    : x(v.x), y(v.y), z(v.z) {}
Vector3::~Vector3() {}

Vector3*
Vector3::clone() const
{
    return new Vector3(*this);
}

// string representation
string Vector3::str(int precision) const
{
    ostringstream ost;
    ost << setprecision(precision);
    ost << "Vector3(" << x << ", " << y << ", " << z << ")";
    return ost.str();
}
string Vector3::vrml_str( double radius ) const
{
    ostringstream ost;
    ost << "\nTransform {";
    ost << "\n  # Vector3 point";
    ost << "\n  translation " << x << " " << y << " " << z;
    ost << "\n  children [";
    ost << "\n    Shape {";
    ost << "\n      appearance Appearance {";
    ost << "\n        material Material { diffuseColor 1 0 0 }";
    ost << "\n      }";
    ost << "\n      geometry Sphere { radius " << radius << " }";
    ost << "\n    }";
    ost << "\n  ]";
    ost << "\n}";
    return ost.str();
}
string Vector3::vtk_str() const
{
    ostringstream ost;
    ost << x << ' ' << y << ' ' << z;
    return ost.str();
}

Vector3& Vector3::mirror_image( const Vector3 &point, const Vector3 &normal )
{
    Vector3 n = unit(normal);
    // Construct tangents to the plane.
    Vector3 different = n + Vector3(1.0, 1.0, 1.0);
    Vector3 t1 = cross(n, different).norm();
    Vector3 t2 = cross(n, t1).norm();
    // Mirror image the vector in a frame local to the plane.
    transform_to_local(n, t1, t2, point);
    x = -x;
    transform_to_global(n, t1, t2, point);
    return *this;
}

/// Rotate the point about the z-axis by the specified angle (in radians).
Vector3& Vector3::rotate_about_zaxis( double dtheta )
{
    double theta = atan2(y,x) + dtheta;
    double r = sqrt(x*x + y*y);
    x = r * cos(theta);
    y = r * sin(theta);
    return *this;
}

// some combined assignment, arithmetic operators.
Vector3& Vector3::operator+=( const Vector3 &v )
{
    x += v.x; y += v.y; z += v.z;
    return *this;
}
Vector3& Vector3::operator-=( const Vector3 &v )
{
    x -= v.x; y -= v.y; z -= v.z;
    return *this;
}
Vector3& Vector3::operator*=( double v )
{
    x *= v; y *= v; z *= v;
    return *this;
}
Vector3& Vector3::operator/=( double v )
{
    x /= v; y /= v; z /= v;
    return *this;
}
Vector3& Vector3::operator=( const Vector3 &v )
{
    x = v.x; y = v.y; z = v.z;
    return *this;
}
Vector3& Vector3::norm()
{
    double mag = vabs(*this);
    if ( mag > 0.0 ) {
	*this /= mag;
    }
    return *this;
}

// Helper functions for the Vector3 class.
//
// Overload stream output for Vector3 objects
ostream& operator<<( ostream &os, const Vector3 &v )
{
    os << v.str();
    return os;
}
double vabs( const Vector3 &v )
{
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}
Vector3 unit( const Vector3 &v )
{
    return v/vabs(v);
}
Vector3 operator+( const Vector3 &v )
{
    return Vector3(v.x, v.y, v.z);
}
Vector3 operator-( const Vector3 &v )
{
    return Vector3(-v.x, -v.y, -v.z);
}
Vector3 operator+( const Vector3 &v1, const Vector3 &v2 )
{
    return Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}
Vector3 operator-( const Vector3 &v1, const Vector3 &v2 )
{
    return Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}
bool equal( const Vector3 &v1, const Vector3 &v2, double tolerance)
{
    return vabs(v1 - v2) < tolerance;
}
Vector3 operator*( const Vector3 &v1, double v2 )
{
    return Vector3(v1.x*v2, v1.y*v2, v1.z*v2);
}
Vector3 operator*( double v1, const Vector3 &v2 )
{
    return Vector3(v2.x*v1, v2.y*v1, v2.z*v1);
}
Vector3 operator/( const Vector3 &v1, double v2 )
{
    return Vector3(v1.x/v2, v1.y/v2, v1.z/v2);
}
double dot( const Vector3 &v1, const Vector3 &v2 )
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
Vector3 cross( const Vector3 &v1, const Vector3 &v2 )
{
    Vector3 v3;
    v3.x = v1.y*v2.z - v2.y*v1.z;
    v3.y = v2.x*v1.z - v1.x*v2.z;
    v3.z = v1.x*v2.y - v2.x*v1.y;
    return v3;
}

//------------------------------------------------------------------------
// Geometry functions that make use of Vector3 objects...

int project_onto_plane( Vector3 &q, const Vector3 &qr,
			const Vector3 &a, const Vector3 &b, const Vector3 &c )
{
    // See Section 7.2 in
    // J. O'Rourke (1998)
    // Computational Geometry in C (2nd Ed.)
    // Cambridge Uni Press 
    // See 3D CFD workbook p17, 25-Jan-2006

    // Define a plane Ax + By + Cz = D using the corners of the triangle abc.
    Vector3 N = cross(a-c, b-c); // N = Vector3(A, B, C)
    double D = dot(a, N);

    double numer = D - dot(q, N);
    double denom = dot(qr, N);

    double tol = 1.0e-12;  // floating point tolerance
    if ( fabs(denom) < tol ) {
	if ( fabs(numer) < tol ) {
	    return 1;  // qr is parallel to the plane and q is on the plane
	} else {
	    return 2;  // qr is parallel to the plane and q is off the plane
	}
    } else {
	q = q + (numer/denom) * qr;
	return 0;  // point q has been projected onto the plane.
    } 
}

int inside_triangle( const Vector3 &p, 
		     const Vector3 &a, const Vector3 &b, const Vector3 &c)
{
    Vector3 n = unit(cross(a-c, b-c)); // normal to plane of triangle
    double area1 = 0.5 * dot(cross(p-a, p-b), n); // projected area of triangle pab
    double area2 = 0.5 * dot(cross(p-b, p-c), n);
    double area3 = 0.5 * dot(cross(p-c, p-a), n);
    // Only a point inside the triangle will have all areas positive.
    if ( area1 > 0.0 && area2 > 0.0 && area3 > 0.0 ) return 1;
    // However, the point may be damned close to an edge but
    // missing because of floating-point round-off.
    double tol = 1.0e-12;
    if ( (fabs(area1) < tol && area2 > -tol && area3 > -tol) || 
	 (fabs(area2) < tol && area1 > -tol && area3 > -tol) || 
	 (fabs(area3) < tol && area1 > -tol && area2 > -tol) ) return 2;
    // Otherwise, the point is outside the triangle.
    return 0;
}
// ---------------------------------------------------------------------------

/// Map space so that a neutral plane wraps onto a cylinder of radius H.
int map_neutral_plane_to_cylinder( Vector3 &p, double H )
{
    // The axis of the hypothetical cylinder coincides with the x-axis thus
    // H is also the distance of the neutral plane above the x-axis.
    // For Hannes Wojciak and Paul Petrie-Repar's turbomachinery grids.
    if ( H <= 0.0 ) return -1;
    double theta = p.y / H;
    double old_z = p.z;
    p.z = old_z * cos(theta);
    p.y = old_z * sin(theta);
    // x remains the same
    return 0;
}

// ---------------------------------------------------------------------------
// Compute some properties of geometric objects.
// These functions have been copied from geometry/source/geom.c
// and then mangled to use the C++ Vector3 class methods.
// This gives much tidier code, almost like IAJ's code.

int quad_properties( const Vector3 &p0, const Vector3 &p1, 
		     const Vector3 &p2, const Vector3 &p3,
		     Vector3 &centroid,
		     Vector3 &n, Vector3 &t1, Vector3 &t2,
		     double &area )
{
    // Compute areas via the cross products.
    Vector3 vector_area = 0.5 * (cross(p0-p3, p2-p3) + cross(p1-p0, p2-p1));
    // unit-normal and area
    n = unit(vector_area);
    area = vabs(vector_area);
    // Tangent unit-vectors: 
    // t1 is parallel to s01, 
    // t2 is normal to n and t1
    t1 = unit(p1-p0);
    t2 = cross(n, t1);
    // Centroid: average mid-points of the diagonals.
    centroid = 0.25 * (p0 + p1 + p2 + p3);
    return 0; 
} /* end quad_properties() */

Vector3 quad_centroid( const Vector3 &p0, const Vector3 &p1, 
		       const Vector3 &p2, const Vector3 &p3 )
{
    Vector3 centroid, n, t1, t2;
    double area;
    quad_properties( p0, p1, p2, p3, centroid, n, t1, t2, area );
    return centroid;
}

Vector3 quad_normal( const Vector3 &p0, const Vector3 &p1, 
		     const Vector3 &p2, const Vector3 &p3 )
{
    Vector3 centroid, n, t1, t2;
    double area;
    quad_properties( p0, p1, p2, p3, centroid, n, t1, t2, area );
    return n;
}

Vector3 quad_tangent1( const Vector3 &p0, const Vector3 &p1, 
		       const Vector3 &p2, const Vector3 &p3 )
{
    Vector3 centroid, n, t1, t2;
    double area;
    quad_properties( p0, p1, p2, p3, centroid, n, t1, t2, area );
    return t1;
}

Vector3 quad_tangent2( const Vector3 &p0, const Vector3 &p1, 
		       const Vector3 &p2, const Vector3 &p3 )
{
    Vector3 centroid, n, t1, t2;
    double area;
    quad_properties( p0, p1, p2, p3, centroid, n, t1, t2, area );
    return t2;
}

double quad_area( const Vector3 &p0, const Vector3 &p1, 
		  const Vector3 &p2, const Vector3 &p3 )
{
    Vector3 centroid, n, t1, t2;
    double area;
    quad_properties( p0, p1, p2, p3, centroid, n, t1, t2, area );
    return area;
}


int tetrahedron_properties( const Vector3 &p0, const Vector3 &p1,
			    const Vector3 &p2, const Vector3 &p3,
			    Vector3 &centroid, double &volume )
{
    volume = dot(p3-p0, cross(p1-p0, p2-p0)) / 6.0;
    // Centroid: average the vertex locations.
    centroid = 0.25 * (p0 + p1 + p2 + p3);
    return 0; 
}

Vector3 tetrahedron_centroid( const Vector3 &p0, const Vector3 &p1,
			      const Vector3 &p2, const Vector3 &p3 )
{
    double volume;
    Vector3 centroid;
    tetrahedron_properties( p0, p1, p2, p3, centroid, volume );
    return centroid; 
}

double tetrahedron_volume( const Vector3 &p0, const Vector3 &p1,
			   const Vector3 &p2, const Vector3 &p3 )
{
    double volume;
    Vector3 centroid;
    tetrahedron_properties( p0, p1, p2, p3, centroid, volume );
    return volume; 
}


int wedge_properties( const Vector3 &p0, const Vector3 &p1,
		      const Vector3 &p2, const Vector3 &p3,
		      const Vector3 &p4, const Vector3 &p5,
		      Vector3 &centroid, double &volume )
{
    double v1, v2, v3;
    Vector3 c1, c2, c3;
    tetrahedron_properties( p0, p4, p5, p3, c1, v1 );
    tetrahedron_properties( p0, p5, p4, p1, c2, v2 );
    tetrahedron_properties( p0, p1, p2, p5, c3, v3 );
    volume = v1 + v2 + v3;
    if ( (volume < 0.0 && fabs(volume) < SMALL_BUT_SIGNIFICANT) ||
	 (volume >= 0.0 && volume < VERY_SMALL_MAGNITUDE) ) {
	// We assume that we have a collapsed wedge; no real problem.
	volume = 0.0;
	// equally-weighted tetrahedral centroids.
	centroid = (c1 + c2 + c3) / 3.0;
	return 0;
    }
    if ( volume < 0.0 ) {
	// Something has gone wrong with our wedge geometry.
	printf("wedge_properties(): significant negative volume: %e\n", volume);
	printf("                    setting volume to zero\n");
	volume = 0.0;
	// equally-weighted tetrahedral centroids.
	centroid = (c1 + c2 + c3) / 3.0;
	return -1;
    }
    // Weight the tetrahedral centroids with their volumes.
    centroid = (c1*v1 + c2*v2 +c3*v3) / volume;
    return 0;
} /* end wedge_properties() */

Vector3 wedge_centroid( const Vector3 &p0, const Vector3 &p1,
			const Vector3 &p2, const Vector3 &p3,
			const Vector3 &p4, const Vector3 &p5)
{
    double volume;
    Vector3 centroid;
    wedge_properties( p0, p1, p2, p3, p4, p5, centroid, volume);
    return centroid; 
}

double wedge_volume( const Vector3 &p0, const Vector3 &p1,
		     const Vector3 &p2, const Vector3 &p3,
		     const Vector3 &p4, const Vector3 &p5)
{
    double volume;
    Vector3 centroid;
    wedge_properties( p0, p1, p2, p3, p4, p5, centroid, volume );
    return volume; 
}


int hexahedron_properties( const Vector3 &p0, const Vector3 &p1,
			   const Vector3 &p2, const Vector3 &p3,
			   const Vector3 &p4, const Vector3 &p5,
			   const Vector3 &p6, const Vector3 &p7,
			   Vector3 &centroid, double &volume )
{
    double v1, v2;
    struct Vector3 c1, c2;
    wedge_properties( p0, p1, p2, p4, p5, p6, c1, v1 );
    wedge_properties( p0, p2, p3, p4, p6, p7, c2, v2 );
    volume = v1 + v2;
    if (volume < VERY_SMALL_MAGNITUDE) {
	printf("hexahedron_properties(): zero or negative volume: %e\n", volume);
	printf("                         setting volume to zero.\n");
	volume = 0.0;
	// equally-weighted prism centroidal values.
	centroid = 0.5 * (c1 + c2);
	return -1;
    }
    // Volume-weight the prism centroidal values.
    centroid = (c1*v1 + c2*v2) / volume;
    return 0; 
} /* end hexahedron_properties() */


double hexahedron_volume( const Vector3 &p0, const Vector3 &p1,
			  const Vector3 &p2, const Vector3 &p3,
			  const Vector3 &p4, const Vector3 &p5,
			  const Vector3 &p6, const Vector3 &p7 )
{
    double volume;
    Vector3 centroid;
    hexahedron_properties( p0, p1, p2, p3, p4, p5, p6, p7, centroid, volume );
    return volume;
}

Vector3 hexahedron_centroid( const Vector3 &p0, const Vector3 &p1,
			     const Vector3 &p2, const Vector3 &p3,
			     const Vector3 &p4, const Vector3 &p5,
			     const Vector3 &p6, const Vector3 &p7 )
{
    double volume;
    Vector3 centroid;
    hexahedron_properties( p0, p1, p2, p3, p4, p5, p6, p7, centroid, volume );
    return centroid;
}

// ----------------------------------------------------------------

int local_frame( Vector3 &v, const Vector3 &n, 
		 const Vector3 &t1, const Vector3 &t2 )
{
    Vector3 v1 = v;
    v.x = dot( v1, n );     /* Normal velocity       */
    v.y = dot( v1, t1 );    /* Tangential velocity 1 */
    v.z = dot( v1, t2 );    /* Tangential velocity 2 */
    return 0;
}

int xyz_frame( Vector3 &v, const Vector3 &n, 
	       const Vector3 &t1, const Vector3 &t2 )
{
    Vector3 v1 = v;
    v.x = v1.x * n.x + v1.y * t1.x + v1.z * t2.x; /* global-x component */
    v.y = v1.x * n.y + v1.y * t1.y + v1.z * t2.y; /* global-y component */
    v.z = v1.x * n.z + v1.y * t1.z + v1.z * t2.z; /* global-z component */
    return 0;
}


//------------------------------------------------------------------------
