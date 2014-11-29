/** \file surface.hh
 *  \ingroup libgeom2
 *  \brief Declarations for the C++ geometric-surface classes.
 *  \author PJ
 *  \version 31-Dec-2005 -- Initial code for base class ParametricSurface and
 *                          derived CoonsPatch class.
 *  \version 23-Jan-2006 -- Added AOPatch.
 *  \version 26-Jan-2005 -- Added TrianglePatch class.
 *  \version 02-Mar-2006 -- Extended TrianglePatch class.
 *  \version 12-Apr-2006 -- subsections added to all surfaces
 *
 */

#ifndef SURFACE_HH
#define SURFACE_HH

#include <string>
#include <iostream>
#include <vector>
#include "geom.hh"
#include "gpath.hh"
#include "nurbs.hh"

using namespace std;

/// \brief Base class for parametric surfaces, p = p(r,s).
///
/// \note A subsection of a path can be specified 
///       by setting suitable values for r0, r1, s0, s1.
class ParametricSurface {
public:
    string label;  ///< May be handy for VRML rendering.
    double r0;     ///< Value of r at left side of surface subsection.
    double r1;     ///< Value of r at right side of subsection.
    double s0;     ///< Value of s at bottom of subsection.
    double s1;     ///< Value of s at top of subsection.

    /// Constructor that should not be called directly.
    ParametricSurface( const string label="", 
		       double r0=0.0, double r1=1.0,
		       double s0=0.0, double s1=1.0 );
    /// Copy constructor (also should not be called directly).
    ParametricSurface( const ParametricSurface &p );
    virtual ~ParametricSurface();
    /// Returns a pointer to a (new) clone of the surface object.
    virtual ParametricSurface* clone() const; 
    /// Also returns a pointer to a (new) clone of the surface object.
    virtual ParametricSurface* copy() const; 
    /// Returns a point on the surface at parametric coordinates (r,s).
    virtual Vector3 eval( double r, double s ) const;
    /// Returns the local tangent (gradient) in the r-parameter direction.
    virtual Vector3 dpdr( double r, double s ) const;
    /// Returns the local tangent (gradient) in the s-parameter direction.
    virtual Vector3 dpds( double r, double s ) const;
    /// Returns a string representation of the Path (much like Python's __str__ method).
    virtual string str() const;
    /// Shifts the whole surface by Cartesian coordinates.
    virtual ParametricSurface* translate( const Vector3 &v );
    /// Shifts the whole surface by Cartesian coordinates.
    virtual ParametricSurface* translate( double vx, double vy, double vz );
    /// Mirror image in the plane defined by a point and a normal vector.
    virtual ParametricSurface* mirror_image( const Vector3 &point, const Vector3 &normal );
    /// Rotate the surface about the z-axis by the specified angle (in radians).
    virtual ParametricSurface* rotate_about_zaxis( double dtheta );
};

ostream& operator<<( ostream &os, const ParametricSurface &p );

/** \brief A surface defined as a blend of 4 bounding paths.
 *
 * The topology of the parametric surface is shown here.
 *
 * \verbatim
 *      1   +-----B-----+            p01----B----p11
 *      ^   |           |             |           |
 *      |   |           |             |           |
 *      s   C           D             C           D
 *      |   |           |             |           |
 *      |   |           |             |           |
 *      0   +-----A-----+            p00----A----p10
 *
 *          0-----r---->1
 * \endverbatim
 *
 * The bounding paths are A, B, C and D while 
 * the corner points are labelled with r and s parameter values.
 */
class CoonsPatch : public ParametricSurface {
public:
    /// Pointers to the bounding paths.
    Path* cA; Path* cB; Path* cC; Path* cD;
    /// The corner points.
    Vector3 p00; Vector3 p10; Vector3 p11; Vector3 p01;
    /// Construct a parametric surface from 4 bounding paths.
    CoonsPatch( const Path &_cA, const Path &_cB, 
		const Path &_cC, const Path &_cD, 
		string label="", 
		double r0 = 0.0, double r1 = 1.0,
		double s0 = 0.0, double s1 = 1.0 );
    /// Construct a surface from 4 corner points, assuming straight-line bounding paths.
    CoonsPatch( const Vector3 p00, const Vector3 p10, 
		const Vector3 p11, const Vector3 p01,
		string label="",
		double r0 = 0.0, double r1 = 1.0,
		double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another CoonsPatch surface.
    CoonsPatch( const CoonsPatch &surf );
    virtual ~CoonsPatch();
    virtual CoonsPatch* clone() const; 
    virtual CoonsPatch* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
    virtual CoonsPatch* translate( const Vector3 &v );
    virtual CoonsPatch* translate( double vx, double vy, double vz );
    virtual CoonsPatch* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual CoonsPatch* rotate_about_zaxis( double dtheta );
};

/** \brief A surface defined between two bounding paths.
 *         Bezier3 curves (normal to the defining paths)
 *         or straight lines (for a ruled surface) are used
 *         to bridge the region between the defining paths.
 *
 * The topology of the parametric surface is shown here.
 *
 * \verbatim
 *      1   +-----B-----+            p01----B----p11
 *      ^   |           |             |           |
 *      |   |           |             |           |
 *      s   |           |             |           |
 *      |   |           |             |           |
 *      |   |           |             |           |
 *      0   +-----A-----+            p00----A----p10
 *
 *          0-----r---->1
 * \endverbatim
 *
 * The bounding paths are A and B while 
 * the corner points are labelled with r and s parameter values.
 *
 * 2014-Nov-19: Peter J. inspired by Wilson's expansion-region surface.
 */
class ChannelPatch : public ParametricSurface {
public:
    /// Pointers to the bounding paths.
    Path* cA; Path* cB; 
    bool ruled;
    bool pure2D;
    /// Construct a parametric surface from 2 bounding paths.
    ChannelPatch( const Path &_cA, const Path &_cB,
		  bool ruled=false, bool pure2D=false,
		  string label="", 
		  double r0 = 0.0, double r1 = 1.0,
		  double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another CoonsPatch surface.
    ChannelPatch( const ChannelPatch &surf );
    virtual ~ChannelPatch();
    virtual ChannelPatch* clone() const; 
    virtual ChannelPatch* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    Path* make_bridging_path( double r ) const;
    virtual string str() const;
    virtual ChannelPatch* translate( const Vector3 &v );
    virtual ChannelPatch* translate( double vx, double vy, double vz );
    virtual ChannelPatch* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual ChannelPatch* rotate_about_zaxis( double dtheta );
};

/** \brief Another surface defined as a blend of 4 bounding paths.
 *
 * The topology of the parametric surface is the same as that of the CoonsPatch.
 * The difference, however, is in the evaluation of points on the surface.
 * Here, points are interpolated within a background mesh that has been fitted
 * with the Area-Orthogonality (AO) elliptic grid generator 
 * described by Patrick M. Knupp "A Robust Elliptic Grid Generator"
 * J. Computational Physics Vol.100 pp409-418 (1992)
 *
 */
class AOPatch : public ParametricSurface {
public:
    CoonsPatch tfi_surface;
    int nx, ny;
    /// Construct a parametric surface from 4 bounding paths.
    AOPatch( const Path &_cA, const Path &_cB, 
	     const Path &_cC, const Path &_cD, 
	     string label="", int nx=20, int ny=20,
	     double r0 = 0.0, double r1 = 1.0,
	     double s0 = 0.0, double s1 = 1.0 );
    /// Construct a surface from 4 corner points, assuming straight-line bounding paths.
    AOPatch( const Vector3 p00, const Vector3 p10, 
	     const Vector3 p11, const Vector3 p01,
	     string label="", int nx=20, int ny=20,
	     double r0 = 0.0, double r1 = 1.0,
	     double s0 = 0.0, double s1 = 1.0  );
    /// Construct a surface from a given CoonsPatch surface.
    AOPatch( const CoonsPatch &tfi_surf, string label="", int nx=20, int ny=20,
	     double r0 = 0.0, double r1 = 1.0,
	     double s0 = 0.0, double s1 = 1.0  );
    /// Construct as a copy of another AO surface.
    AOPatch( const AOPatch &surf );
    virtual ~AOPatch();
    virtual AOPatch* clone() const; 
    virtual AOPatch* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
    virtual AOPatch* translate( const Vector3 &v );
    virtual AOPatch* translate( double vx, double vy, double vz );
    virtual AOPatch* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual AOPatch* rotate_about_zaxis( double dtheta );
private:
    Vector3 **bgmesh;  // background mesh
    int mesh_ok; // flag to indicate that the mesh memory is allocated
    void allocate_background_mesh();
    void compute_background_mesh();
    void free_background_mesh();
};


/** \brief A surface defined by a mesh of quadrilateral facets in a regular grid.
 *
 */
class MeshPatch : public ParametricSurface {
public:
    vector<Vector3> p; // vertices of the quads making up the mesh.
    int ni; // number of vertices in the r-coordinate direction.
    int nj; // number of vertices in the s-coordinate direction.
    /// Construct a parametric surface from a collection of points in a 2D array.
    MeshPatch( const vector<Vector3*> &_p, int _ni, int _nj,
	       string label="",
	       double r0 = 0.0, double r1 = 1.0,
	       double s0 = 0.0, double s1 = 1.0  );
#ifndef SWIG
    // We need to hide this constructor from SWIG to avoid shadowing
    // the previous constructor.
    MeshPatch( const vector<Vector3> &_p, int _ni, int _nj,
	       string label="",
	       double r0 = 0.0, double r1 = 1.0,
	       double s0 = 0.0, double s1 = 1.0  );
#endif
    /// Construct the volume from data in a file.
    MeshPatch( const string file_name, int vtk_format = 1,
		const string label="",
		double r0 = 0.0, double r1 = 1.0,
		double s0 = 0.0, double s1 = 1.0);
    /// Construct as a copy of another MeshPatch surface.
    MeshPatch( const MeshPatch &surf );
    virtual ~MeshPatch();
    virtual MeshPatch* clone() const; 
    virtual MeshPatch* copy() const; 
    virtual Vector3 eval( double r, double s) const;
    virtual string str() const;
    virtual MeshPatch* translate( const Vector3 &v );
    virtual MeshPatch* translate( double vx, double vy, double vz );
    virtual MeshPatch* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual MeshPatch* rotate_about_zaxis( double dtheta );
};

/** \brief A surface defined by triangular facets.
 *
 */
class TrianglePatch : public ParametricSurface {
public:
    vector<Vector3> p; // vertices of the triangles and the Polyline edges
    vector<int> itri; // indices of the vertices that form the triangles
    vector<int> iA, iB, iC, iD; // indices of the points forming the boundaries  
    /// Construct a parametric surface from a collection of points,
    /// the indices of those points that from the collection of triangular facets,
    /// and the indices of the points that form the 4 bounding paths.
    TrianglePatch( const vector<Vector3*> &_p, const vector<int> &_itri,
		   const vector<int> &_iA, const vector<int> &_iB, 
		   const vector<int> &_iC, const vector<int> &_iD, 
		   string label="",
		   double r0 = 0.0, double r1 = 1.0,
		   double s0 = 0.0, double s1 = 1.0  );
    /// Construct a TrianglePatch from another ParametricSurface.
    /// @param Nr : number of panels in the r-parameter direction
    /// @param Ns : number of panels in the s-parameter direction
    TrianglePatch( const ParametricSurface *surf, int Nr=2, int Ns=2,
		   string label="",
		   double r0 = 0.0, double r1 = 1.0,
		   double s0 = 0.0, double s1 = 1.0  );
    /// Construct as a copy of another TrianglePatch surface.
    TrianglePatch( const TrianglePatch &surf );
    virtual ~TrianglePatch();
    virtual TrianglePatch* clone() const; 
    virtual TrianglePatch* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
    virtual TrianglePatch* translate( const Vector3 &v );
    virtual TrianglePatch* translate( double vx, double vy, double vz );
    virtual TrianglePatch* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual TrianglePatch* rotate_about_zaxis( double dtheta );
    /// Appends another TrianglePatch to this one.
    /// @returns: 0 if all went OK, 1 otherwise.
    /// For the moment, the only allowable extensions are those
    /// that extend in the direction of one of the (r or s) parameter directions.
    /// That is, the additional patch may join only edge D or B.
    int add( const TrianglePatch &surf );
private:
    CoonsPatch *tfi_surface; // An overlaid parametric surface to get the query points.
    int ntri; // number of triangular facets
    void setup_tfi_surface();
};

// Returns 1 if the specified edges match, point-to-point; returns 0 otherwise.
// Used internally when adding one patch to another.
bool edges_are_matched( const TrianglePatch &surf0, const vector<int> &edge0, 
			const TrianglePatch &surf1, const vector<int> &edge1,
			int opposite_direction );


/** \brief A Tensor-Product Bezier surface.
 *
 */
class BezierPatch : public ParametricSurface {
public:
    // The control net is defined by points Q[i,j] with 0 <= i <= n, 0 <= j <= m.
    vector<Vector3> Q; // Points are stored as Q[i*(m+1) + j].
    int n; // order of the polynomial in the r-parameter direction
    int m; // order of the polynonial in the s-parameter direction 
    /// Construct a parametric surface from a collection of points.
    BezierPatch( const vector<Vector3*> &_Q, int n, int m, string label="",
		 double r0 = 0.0, double r1 = 1.0,
		 double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another BezierPatch surface.
    BezierPatch( const BezierPatch &surf );
    virtual ~BezierPatch();
    virtual BezierPatch* clone() const; 
    virtual BezierPatch* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
    virtual BezierPatch* translate( const Vector3 &v );
    virtual BezierPatch* translate( double vx, double vy, double vz );
    virtual BezierPatch* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual BezierPatch* rotate_about_zaxis( double dtheta );
};


/** \brief A surface defined by revolving a Path about the x-axis.
 *
 */
class RevolvedSurface : public ParametricSurface {
public:
    Path *path;  // Path that is to be revolved.
    RevolvedSurface( const Path* pathptr, string label="",
		     double r0 = 0.0, double r1 = 1.0,
		     double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another RevolvedSurface surface.
    RevolvedSurface( const RevolvedSurface &surf );
    virtual ~RevolvedSurface();
    virtual RevolvedSurface* clone() const; 
    virtual RevolvedSurface* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
};


/** \brief A surface defined by mapping one surface onto another (true) surface.
 *
 */
class MappedSurface : public ParametricSurface {
public:
    ParametricSurface *query_surf;  // Sample points are generated on this surface.
    ParametricSurface *true_surf;   // Final points are on this surface.
    MappedSurface( const ParametricSurface* qsptr, 
		   const ParametricSurface* tsptr,
		   string label="",
		   double r0 = 0.0, double r1 = 1.0,
		   double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another MappedSurface.
    MappedSurface( const MappedSurface &surf );
    virtual ~MappedSurface();
    virtual MappedSurface* clone() const; 
    virtual MappedSurface* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
};


/** \brief A surface transformed from Cartesian space to Polar space.
 *
 */
class PolarSurface : public ParametricSurface {
public:
    ParametricSurface *original_surf;  // in Cartesian space.
    double H; // Height of the neutral plane above the x-axis.
    PolarSurface( const ParametricSurface* surf, 
		  double H,
		  string label="",
		  double r0 = 0.0, double r1 = 1.0,
		  double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another PolarSurface.
    PolarSurface( const PolarSurface &surf );
    virtual ~PolarSurface();
    virtual PolarSurface* clone() const; 
    virtual PolarSurface* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
};


/** \brief A surface defined through a parametric volume.
 *
 * SurfaceThruVolume added for turbomachinery work, 19-Oct-2008, PJ
 * but will be more generally useful.
 */
class ParametricVolume;
class SurfaceThruVolume : public ParametricSurface {
public:
    ParametricVolume *pvol;
    BivariateFunction *fr;
    BivariateFunction *fs;
    BivariateFunction *ft;
    SurfaceThruVolume( const ParametricVolume& _pvol, 
		       const BivariateFunction& _fr,
		       const BivariateFunction& _fs,
		       const BivariateFunction&_ft,
		       string label="",
		       double r0 = 0.0, double r1 = 1.0,
		       double s0 = 0.0, double s1 = 1.0 );
    /// Construct as a copy of another PolarSurface.
    SurfaceThruVolume( const SurfaceThruVolume &surf );
    virtual ~SurfaceThruVolume();
    virtual SurfaceThruVolume* clone() const; 
    virtual SurfaceThruVolume* copy() const; 
    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
};


/** \brief A Tensor-Product Nurbs surface.
 *
 */
class NurbsSurface : public ParametricSurface {
public:
    vector<vector<Vector3> > Q; // Net of control points
    vector<vector<double> > w;  // weights associated with control points
    int p;                      // degree of surface in u-direction
    vector<double> U;           // knot vector for u-direction
    int q;                      // degree of surface in v-direction
    vector<double> V;           // knot vector in v-direction

    NurbsSurface( const vector<vector<Vector3> > &Q, const vector<vector<double> > &w,
		  int p, const vector<double> &U, int q, const vector<double> &V,
		  string label="",
		  double r0=0.0, double r1=1.0,
		  double s0=0.0, double s1=1.0 );
    NurbsSurface( const vector<vector<Mapped_point> > &Qw,
		  int p, const vector<double> &U, int q, const vector<double> &V,
		  string label="",
		  double r0=0.0, double r1=1.0,
		  double s0=0.0, double s1=1.0 );
    NurbsSurface( const NurbsSurface &n );

    virtual ~NurbsSurface();
    virtual NurbsSurface* clone() const; 
    virtual NurbsSurface* copy() const; 

    double map_r2u(double r) const;
    double map_s2v(double s) const;

    virtual Vector3 eval( double r, double s ) const;
    virtual string str() const;
    virtual NurbsSurface* translate( const Vector3 &v );
    virtual NurbsSurface* translate( double vx, double vy, double vz );
    virtual NurbsSurface* mirror_image( const Vector3 &point, const Vector3 &normal );
						
private:
    vector<vector<Mapped_point> > Qw_;
    double umin_, umax_, vmin_, vmax_;
};

int write_STL(const ParametricSurface &s, int nr, int ns, std::string fname);


#endif
