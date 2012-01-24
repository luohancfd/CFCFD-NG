/** \file volume.hh
 *  \ingroup libgeom2
 *  \brief Declarations for the C++ geometric-volume classes.
 *  \author PJ
 *  \version 12-Jan-2006 -- initial code for base class ParametricVolume.
 *  \version 13-Jan-2006 -- derived WireFrameVolume, SimpleBoxVolume classes.
 *  \version 17-Jan-2006 -- get rid of std::vector<T> from WireFrameVolume and SimpleBoxVolume
 *  \version 12-Apr-2006 -- Add subsections.
 */

#ifndef VOLUME_HH
#define VOLUME_HH

#include <string>
#include <iostream>
#include <vector>
#include "geom.hh"
#include "gpath.hh"
#include "surface.hh"
using namespace std;

/** \brief A parametric volume has the topology of a hexahedron. 
 */
class ParametricVolume {
public:
    string label;
    // Order of Surfaces and corner points as found in Elmer.
    ParametricSurface* south;
    ParametricSurface* bottom;
    ParametricSurface* west;
    ParametricSurface* east;
    ParametricSurface* north;
    ParametricSurface* top;
    Vector3 p000, p100, p110, p010; // corners around the bottom surface
    Vector3 p001, p101, p111, p011; // corners around the top surface
    // Ranges of the subsection of each parameter.
    double r0, r1;
    double s0, s1;
    double t0, t1;
    /// Construct the volume from 6 surfaces.
    ParametricVolume( const ParametricSurface* faceN, 
		      const ParametricSurface* faceE,
		      const ParametricSurface* faceS,
		      const ParametricSurface* faceW,
		      const ParametricSurface* faceT,
		      const ParametricSurface* faceB,
		      string label="",
		      double r0 = 0.0, double r1 = 1.0,
		      double s0 = 0.0, double s1 = 1.0,
		      double t0 = 0.0, double t1 = 1.0 );
    /// Construct the volume from a collection of 6 surfaces.
    ParametricVolume( const vector<ParametricSurface*> &face, const string label="",
		      double r0 = 0.0, double r1 = 1.0,
		      double s0 = 0.0, double s1 = 1.0,
		      double t0 = 0.0, double t1 = 1.0 );
    /// Constructor for the case when we don't yet have the surfaces.
    ParametricVolume( const string label="",
		      double r0 = 0.0, double r1 = 1.0,
		      double s0 = 0.0, double s1 = 1.0,
		      double t0 = 0.0, double t1 = 1.0 );
    /// Construct as a copy of another parametric volume. 
    ParametricVolume( const ParametricVolume &vol );
    virtual ~ParametricVolume();
    /// Returns a pointer to a newly cloned volume.
    virtual ParametricVolume* clone() const; 
    /// Also returns a pointer to a newly cloned volume.
    virtual ParametricVolume* copy() const; 
    /// Returns a point at parametric coordinates (r,s,t).
    virtual Vector3 eval( double r, double s, double t ) const;
    /// Returns a string representation of the Path (much like Python's __str__ method).
    virtual string str() const;
    /// Shifts the whole path by Cartesian coordinates.
    virtual ParametricVolume* translate( const Vector3 &v );
    /// Shifts the whole path by Cartesian coordinates.
    virtual ParametricVolume* translate( double vx, double vy, double vz );
    /// Mirror image in the plane defined by a point and a normal vector.
    virtual ParametricVolume* mirror_image( const Vector3 &point, const Vector3 &normal );
    /// Rotate the volume about the z-axis by the specified angle (in radians).
    virtual ParametricVolume* rotate_about_zaxis( double dtheta );
protected:
    int set_and_check_corners();
};

// Helper functions
/// Writes a string representation of the volume to the output stream.
ostream& operator<<( ostream &os, const ParametricVolume &p );

// Special cases.
/// A volume defined by its 12 edges.
class WireFrameVolume : public ParametricVolume {
public:
    /// Construct the volume from 12 (edge) paths.
    WireFrameVolume( const Path &c01, const Path &c12, const Path &c32, const Path &c03, 
		     const Path &c45, const Path &c56, const Path &c76, const Path &c47,
		     const Path &c04, const Path &c15, const Path &c26, const Path &c37,
		     const string label="",
		     double r0 = 0.0, double r1 = 1.0,
		     double s0 = 0.0, double s1 = 1.0,
		     double t0 = 0.0, double t1 = 1.0 );
    /// Construct the volume by extruding a surface along a path.
    ///
    /// The extrude_path forms one of the edges of the volume with the t=0 point
    /// on the path becoming p000 for the volume. The surface is translated so that
    /// its (0,0) point is at p000 and it becomes the WEST, SOUTH or BOTTOM face
    /// of the volume for directions i,j and k, respectively.
    WireFrameVolume( const CoonsPatch &base_surf, 
		     Path &extrude_path, 
		     const string direction="k", const
		     string label="",
		     double r0 = 0.0, double r1 = 1.0,
		     double s0 = 0.0, double s1 = 1.0,
		     double t0 = 0.0, double t1 = 1.0 );
};

/// A volume defined as a box with straight-line edges.
class SimpleBoxVolume : public ParametricVolume {
public:
    /// Construct the volume from the 8 corner points.
    SimpleBoxVolume( const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,
		     const Vector3 &p4, const Vector3 &p5, const Vector3 &p6, const Vector3 &p7, 
		     const string label="",
		     double r0 = 0.0, double r1 = 1.0,
		     double s0 = 0.0, double s1 = 1.0,
		     double t0 = 0.0, double t1 = 1.0  );
};

/// A volume defined by a 3D mesh of hexahedral cells.
class MeshVolume : public ParametricVolume {
public:
    vector<Vector3> p; // vertices of the 3D mesh
    int ni; // number of vertices in the r-coordinate direction
    int nj; // number of vertices in the s-coordinate direction
    int nk; // number of vertices in the t-coordinate direction
    /// Construct the volume from a list of vertices.
    MeshVolume( const vector<Vector3*> _p, int _ni, int _nj, int _nk,
		const string label="",
		double r0 = 0.0, double r1 = 1.0,
		double s0 = 0.0, double s1 = 1.0,
		double t0 = 0.0, double t1 = 1.0  );
    /// Construct the volume from data in a file.
    MeshVolume( const string file_name, int vtk_format = 1,
		const string label="",
		double r0 = 0.0, double r1 = 1.0,
		double s0 = 0.0, double s1 = 1.0,
		double t0 = 0.0, double t1 = 1.0  );
    ~MeshVolume();
    /// Returns a point at parametric coordinates (r,s,t).
    virtual Vector3 eval( double r, double s, double t ) const;
    /// Returns a string representation of the Path (much like Python's __str__ method).
    virtual string str() const;
    /// Shifts the whole path by Cartesian coordinates.
    virtual MeshVolume* translate( const Vector3 &v );
    /// Reverses the path in parameter space. new(t) == old(1.0-t)
    virtual MeshVolume* translate( double vx, double vy, double vz );
    /// Mirror image in the plane defined by a point and a normal vector.
    virtual MeshVolume* mirror_image( const Vector3 &point, const Vector3 &normal );
    /// Rotate the volume about the z-axis by the specified angle (in radians).
    virtual MeshVolume* rotate_about_zaxis( double dtheta );
protected:
    int set_surfaces();
};

#endif
