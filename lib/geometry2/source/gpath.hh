/** \file gpath.hh
 *  \ingroup libgeom2
 *  \brief Declarations for the C++ geometric-path classes.
 *  \author PJ
 *  \version 31-Dec-2005 -- initial code for base class Path and
 *                          derived Line class
 *
 */

#ifndef GPATH_HH
#define GPATH_HH

#include <string>
#include <iostream>
#include <vector>
#include "geom.hh"
#include "../../nm/source/fobject.hh"
#include "nurbs.hh"
using namespace std;

/** \brief Base class for the parametric path objects. 
 *
 * Our main interest is being able to evaluate points on the path
 * as parameter t varies from 0.0 to 1.0.
 * 
 * \note A subsection of a path can be specified 
 *       by setting suitable values of t0 and t1.
 * \note We shouldn't ever construct an object of this base class.
 *       This class exists only so that we can deal with collections
 *       of derived-class objects.
 */
class Path {
public:
    string label; //!< A label may be handy when looking at the VRML or SVG rendering.
    double t0;    //!< start of subsection of the full path
    double t1;    //!< end of subsection of the full path
    Path( const string label="", double t0=0.0, double t1=1.0 );
    Path( const Path &p );
    virtual ~Path();
    /// Returns a pointer to a (new) cloned path.
    virtual Path* clone() const; 
    /// Returns a pointer to a (new) cloned path, optionally reversed in direction.
    virtual Path* copy(int direction=1) const; 
    /// Returns the position at parametric value t.
    virtual Vector3 eval( double t ) const;
    /// Returns the parametric value t corresponding to a certain position.
    virtual double locate(const Vector3 &p, int &result_flag,
			  double tolerance=1.0e-6, int max_iterations=1000);
    /// Returns the local tangent (gradient) with respect to the t-parameter.
    virtual Vector3 dpdt( double t ) const;
    /// Returns the geometric length of the path (or subsection).
    virtual double length() const;
    /// Returns the geometric length of a subsection from ta --> tb.
    /// Note that this is not the same subrange as t0->t1 for the Path.
    virtual double partial_length(double ta, double tb) const;
    /// Returns a Vector3 of the location of a length along the path
    virtual Vector3 point_from_length(double length, double &t) const;
    /// Returns a string representation of the Path (much like Python's __str__ method).
    virtual string str() const;
    /// Shifts the whole path by vector displacement v.
    virtual Path* translate( const Vector3 &v );
    /// Shifts the whole path by Cartesian coordinates.
    virtual Path* translate( double vx, double vy, double vz );
    /// Reverses the path in parameter space. new(t) == old(1.0-t)
    virtual Path* reverse();
    /// Mirror image in the plane defined by a point and a normal vector.
    virtual Path* mirror_image( const Vector3 &point, const Vector3 &normal );
    /// Rotate the path about the z-axis by the specified angle (in radians).
    virtual Path* rotate_about_zaxis( double dtheta );
};

ostream& operator<<( ostream &os, const Path &p );

/** \brief Straight-line segments are defined by end points. */
class Line : public Path {
public:
    Vector3 a; // beginning (t = 0)
    Vector3 b; // end (t = 1)
    /// \brief Line constructed from start-point a (t=0.0) to end-point b (t=1.0).
    Line( const Vector3 &a, const Vector3 &b, 
	  const string label="", double t0=0.0, double t1=1.0 );
    /// \brief Line constructed as a copy of another Line.
    Line( const Line &line );
    virtual ~Line();
    virtual Line* clone() const; 
    virtual Line* copy(int direction=1) const; 
    virtual Vector3 eval( double t ) const;
    virtual double length() const;
    virtual string str() const;
    virtual Line* translate( const Vector3 &v );
    virtual Line* translate( double vx, double vy, double vz );
    virtual Line* reverse();
    virtual Line* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual Line* rotate_about_zaxis( double dtheta );
};


/** \brief Circular arc. */
class Arc : public Path {
public:
    Vector3 a; ///< beginning point (t = 0)
    Vector3 b; ///< end point (t = 1)
    Vector3 c; ///< centre of curvature
    /// Arc constructed from start-point a, end-point b and centre of curvature c.
    Arc( const Vector3 &a, const Vector3 &b, const Vector3 &c, 
	 const string label="", double t0=0.0, double t1=1.0 );
    /// Arc constructed as a copy.
    Arc( const Arc &arc );
    virtual ~Arc();
    virtual Arc* clone() const; 
    virtual Arc* copy(int direction=1) const; 
    virtual Vector3 eval( double t ) const;
    virtual double length() const;
    virtual string str() const;
    virtual Arc* translate( const Vector3 &v );
    virtual Arc* translate( double vx, double vy, double vz );
    virtual Arc* reverse();
    virtual Arc* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual Arc* rotate_about_zaxis( double dtheta );
    int evaluate_position_and_length( double t, Vector3& loc, double &L ) const;
};


Vector3 locate_centre( double s );
double error_in_radius( double s );

/** \brief Circular arc constructed through three points.
 *
 * Underneath, it is actually stored as start-point, end-point and centre.
 */
class Arc3 : public Arc {
public:
    /// Arc constructed through three points, starting at a, 
    /// going through b and endin at c.
    Arc3( const Vector3 &a, const Vector3 &b, const Vector3 &c, 
	  const string label="", double t0=0.0, double t1=1.0 );
    /// Arc constructed as a copy of another Arc.
    Arc3( const Arc3 &arc );
    virtual ~Arc3();
};


/** \brief Helix. */
class Helix : public Path {
public:
    Vector3 a0; ///< beginning point on local z-axis (t = 0)
    Vector3 a1; ///< end point on local z-axis (t = 1)
    Vector3 xdsh; ///< local x-axis, unit vector
    Vector3 ydsh; ///< local y-axis, unit vector
    Vector3 zdsh; ///< local z-axis, unit vector
    double r0, r1; ///< starting and ending radii
    double theta01; ///< angle (in radians) from starting point to ending point
    // assuming the right-hand screw convention.
    /// Helix constructed from fundamental quantities.
    Helix( const Vector3 &a0, const Vector3 &a1, 
	   const Vector3 &xlocal, double r0, double r1, double dtheta, 
	   const string label="", double t0=0.0, double t1=1.0 );
    /// Helix constructed from pnt_start to pnt_end
    /// about an axis from axis0 to axis1.
    Helix( const Vector3 &point_start, const Vector3 &point_end, 
	   const Vector3 &axis0, const Vector3 &axis1, 
	   const string label="", double t0=0.0, double t1=1.0 );
    /// Helix constructed as a copy.
    Helix( const Helix &h );
    virtual ~Helix();
    virtual Helix* clone() const; 
    virtual Helix* copy(int direction=1) const; 
    virtual Vector3 eval( double t ) const;
    virtual string str() const;
    virtual Helix* translate( const Vector3 &v );
    virtual Helix* translate( double vx, double vy, double vz );
    virtual Helix* reverse();
    virtual Helix* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual Helix* rotate_about_zaxis( double dtheta );
};


/** \brief Bezier curve of arbitrary order.
 */
class Bezier : public Path {
public:
    vector<Vector3> B; ///< collection of control points
    int arc_length_param_flag;
    /// Construct the curve from a collection of points.
    Bezier( const vector<Vector3> &B, string label="", 
	    double t0=0.0, double t1=1.0, int arc_length_p=0 );
    /// Construct the curve as a copy of another Bezier curve.
    Bezier( const Bezier &bez );
    virtual ~Bezier();
    virtual Bezier* clone() const; 
    virtual Bezier* copy(int direction=1) const; 
    /// Add another point to the end of the collection. (Will increase order of curve.)
    Bezier* add_point( const Vector3 &p );
    Vector3 raw_eval( double t ) const;
    virtual Vector3 eval( double t ) const;
    Vector3 dpdt( double t ) const;
    // virtual double length() const; Get it from Path
    virtual string str() const;
    virtual Bezier* translate( const Vector3 &v );
    virtual Bezier* translate( double vx, double vy, double vz );
    virtual Bezier* reverse();
    virtual Bezier* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual Bezier* rotate_about_zaxis( double dtheta );

    // Added: R.J. Gollan for use with XBezier
    Vector3 d2pdt2( double t ) const;

private:
    // Don't expect users to access these
    void set_deriv_control_points();
    void set_arc_length_vector();
    double t_from_arc_length(double t) const;
    vector<double> arc_length;
    int n_arc_length; // Number of pieces in the arc_length vector, if used.
    vector<Vector3> C; // storage of control points
                       // for quick evaluation of first derivative
    vector<Vector3> D; // storage of control points
                       // for quick evaluation of second derivative
};

/** \brief NURBS curve.
 */
class Nurbs : public Path {
public:
    vector<Vector3> P; // control points
    vector<double> w;  // weights
    int p;             // degree
    vector<double> U;  // knot vector
    vector<Mapped_point> Pw; // control points in 4D space
    vector<int> inf_list; // list of infinite control points
  
    //-- Constructors
    Nurbs( const vector<Vector3> &P, const vector<double> w, 
	   int p, const vector<double> U,
	   string label="", double t0=0.0, double t1=1.0 );
    Nurbs( const vector<Vector3> &P, const vector<double> w,
	   int p, const vector<double> U, const vector<int> inf_list,
	   string label="", double t0=0.0, double t1=1.0 );

    Nurbs( const vector<Mapped_point> &Pw, int p,
	   const vector<double> U,
	   string label="", double t0=0.0, double t1=1.0 );

    Nurbs( const Nurbs &n );

    //-- Destructor
    virtual ~Nurbs();

    //-- Clone
    virtual Nurbs* clone() const;
    
    //-- General services/methods
    Vector3 eval(double t) const;
    string str() const;

    Vector3 dpdu( double t ) const;
    Vector3 d2pdu2( double t ) const;

    double map_t2u(double t) const;

    //-- Modifiers (but leave oringal untouched.)
    Nurbs knot_insertion(double u, int k, int s, int r) const;
    Nurbs knot_refinement(const std::vector<double> &X) const;
    
private:
    // Store these to help calculate linear re-parmeterization
    // in terms of t (0 --> 1).
    double umin_;
    double umax_;
};


/** \brief A Polyline consists of a number of Path elements.
 *
 * \note It is possible that the Path elements are disjoint.
 */
class Polyline : public Path {
public:
    vector<Path*> seg; ///< collection of pointers to Path segments
    // We have to use pointers for the Path segments because they
    // may be of varying type: Line, Arc, Bezier, etc.
    vector<double> t_seg; ///< collection of segment break-points (in parameter t)
    int arc_length_param_flag;
    /// Construct the Polyline from a collection of Path segments.
    Polyline( const vector<Path*> &segments, 
	      string label="",
	      double t0=0.0, double t1=1.0, int arc_length_p=0 );
    /// Construct the Polyline as Path segments between supplied points.
    /// Note that type_of_segments is required, just give a value of 0.
    /// This odd requirement is here to avoid shadowing of the constructor for SWIG.
    Polyline( const vector<Vector3*> &points, int type_of_segments,
	      string label="", int arc_length_p=0 );
    /// A constructor for when we don't yet have the segments (e.g. Spline)
    Polyline( int nseg, string label="", double t0=0.0, double t1=1.0,
	      int arc_length_p=0);
    /// Construct as a copy of another Polyline.
    Polyline( const Polyline &pline );
    virtual ~Polyline();
    virtual Polyline* clone() const;
    virtual Polyline* copy(int direction=1) const; 
    Polyline* add_segment( const Path *segment, int direction=1 );
    Vector3 raw_eval( double t ) const;
    virtual Vector3 eval( double t ) const;
    virtual double length() const;
    virtual string str() const;
    virtual Polyline* translate( const Vector3 &v );
    virtual Polyline* translate( double vx, double vy, double vz );
    virtual Polyline* reverse();
    virtual Polyline* mirror_image( const Vector3 &point, const Vector3 &normal );
    virtual Polyline* rotate_about_zaxis( double dtheta );
protected:
    void reset_breakpoints();
private:
    void set_arc_length_vector();
    double t_from_arc_length(double t) const;
    vector<double> arc_length;
    int n_arc_length; // Number of pieces in the arc_length vector, if used.
};


/// \brief A Spline curve is a specialized Polyline.
class Spline : public Polyline {
public:
    /// Construct the Spline through a collection of points.
    Spline( const vector<Vector3> &p, string label="", 
	    double t0=0.0, double t1=1.0, double tolerance=1.0e-10 );
    /// Construct as a copy of another Spline.
    Spline( const Spline &spl );
    Vector3 eval_from_x(double x) const;
    Vector3 eval_from_y(double y) const;
private:
    double Fx(double t, double x_) const;
    double Fy(double t, double y_) const;
};


/// \brief A path that is part of a surface.
class ParametricSurface;
class PathOnSurface : public Path {
public:
    ParametricSurface* surf; ///< point_position = surf(r,s)
    UnivariateFunction* fr;  ///< r = fr(t)
    UnivariateFunction* fs;  ///< s = fs(t)
    /// Construct the Path from an underlying surface, and two parametric functions.
    PathOnSurface( const ParametricSurface &_surf, 
		   const UnivariateFunction &_fr,
		   const UnivariateFunction &_fs,
		   string label="", double t0=0.0, double t1=1.0 );
    /// Construct the Path as a copy of another PathOnSurface.
    PathOnSurface( const PathOnSurface &p );
    virtual ~PathOnSurface();
    virtual PathOnSurface* clone() const; 
    virtual Vector3 eval( double t ) const;
    virtual string str() const;
    virtual PathOnSurface* translate( const Vector3 &v );
    virtual PathOnSurface* translate( double vx, double vy, double vz );
//     virtual PathOnSurface* reverse();
//     virtual PathOnSurface* mirror_image( const Vector3 &point, const Vector3 &normal );
};

// -------------------------------------------------------------------------

/// \brief Paths in a space that has a neutral-plane wrapped around a cylinder.
///
/// The axis of the hypothetical cylinder coincides with the x-axis thus
/// H is also the distance of the neutral plane above the x-axis.
/// For Hannes Wojciak and Paul Petrie-Repar's turbomachinery grids.
class PolarPath : public Path {
public:
    Path *original_path;
    double H; // height of the neutral-plane above the (x,y)-plane.
    PolarPath( const Path &pth, double H, 
	 const string label="", double t0=0.0, double t1=1.0 );
    PolarPath( const PolarPath &ppth );
    virtual ~PolarPath();
    virtual PolarPath* clone() const; 
    virtual Vector3 eval( double t ) const;
    virtual string str() const;
    virtual PolarPath* translate( const Vector3 &v );
    virtual PolarPath* translate( double vx, double vy, double vz );
    virtual PolarPath* reverse();
    virtual PolarPath* mirror_image( const Vector3 &point, const Vector3 &normal );
};


//---------------------------------------------------------------
// Paths with the form y(x) and linear mapping between x --> t

// --------------------------------------------------------------
// XPath - a base class.
// Objects of this class are not intended for use.
// It is provided for polymorphism purposes.
//
// Added by: Rowan J. Gollan
// Date: Feb-2008

class XPath : public Path {
public:
    double x0;
    double x1;
    /// Default constructor
    XPath();
    /// Normal constructor
    XPath(double x0, double x1, const string label="");
    /// Copy constructor
    XPath(const XPath &x);
    /// Destructor
    virtual ~XPath();
    virtual XPath* clone() const = 0;

    virtual Vector3 eval(double t) const;

    virtual Vector3 dpdt(double t) const;

    virtual Vector3 xeval(double x) const = 0;
    /// first derivative evaluation
    virtual double dydx(double x) const = 0;
    /// second derivative evaulation
    virtual double d2ydx2(double x) const = 0;
    
    /// Mapping t --> x
    virtual double map_t2x(double t) const
    { return (1.0 - t)*x0 + t*x1; }
    
    /// Mapping x --> t
    virtual double map_x2t(double x) const
    { return (x - x0)/(x1 - x0); }

    virtual double locate_x(double y, int &result_flag);
    
    virtual string lua_str(const string name) const = 0;
};

// XPoly
// Added by: Rowan J. Gollan
// Date: Feb-2008

class XPoly : public XPath {
public:
    vector<double> B;
    vector<double> D;
    XPoly();
    XPoly(double x0, double x1, vector<double> &B,
	  const string label="");
    XPoly( const XPoly &x );
    virtual ~XPoly();
    XPoly& operator=(const XPoly &x);
    /// Returns a pointer to a (new) cloned path.
    XPoly* clone() const; 
    /// Returns the position at parametric value t.
    //    Vector3 eval( double t ) const;
    /// Returns the local tangent (gradient) with respect to the t-parameter.
    //    Vector3 dpdt( double t ) const;
    string str() const;

    // Evaluate the path given x as parameter
    Vector3 xeval(double x) const;
    // Evaluate the first derivative w.r.t. x
    double dydx(double x) const;
    // Evaluate the second derivative w.r.t. x
    double d2ydx2(double x) const;

    string lua_str(const string name) const;

};

// XSpline
// Added by: Rowan J. Gollan
// Date: Feb-2008

class XSpline : public XPath {
public:
    vector<XPoly> X;
    vector<double> x_seg;
    XSpline(double x0, double x1, vector<XPoly> &X,
	    const string label="");
    XSpline( const XSpline &x );
    virtual ~XSpline();
    /// Returns a pointer to a (new) cloned path.
    XSpline* clone() const; 
    /// Returns the position at parametric value t.
    //    Vector3 eval( double t ) const;
    /// Returns the local tangent (gradient) with respect to the t-parameter.
    //    Vector3 dpdt( double t ) const;

    // Evaluate the path given x as parameter
    Vector3 xeval(double x) const;
    // Evaluate the first derivative w.r.t. x
    double dydx(double x) const;
    // Evaluate the second derivative w.r.t. x
    double d2ydx2(double x) const;

    // Find the value of x where dydx is a certain value
    Vector3 locate_dydx(double dydx, int &result_flag) const;
    string lua_str(const string name) const;
};

double error_in_dydx(double x);


// XBezier
// This is a wrapper around a Bezier class but
// it is re-parameterized in terms of x.
//
// In other words, this will only work
// with Bezier curves that have single-valued
// functions of x.
//
// Added by: Rowan J. Gollan
// Date: 20-Feb-2008

class XBezier : public XPath {
public:
    Bezier bez;
    XBezier( double x0, double x1, const vector<Vector3> &B,
	     const string label="" );
    XBezier( const XBezier &x );
    virtual ~XBezier();
    /// Returns a pointer to a new cloned path.
    XBezier* clone() const;

    Vector3 eval( double t ) const;
    
    Vector3 dpdt( double t ) const;

    string str() const;

    Vector3 xeval(double x) const;
    
    double dydx(double x) const;
    
    double d2ydx2(double x) const;

    string lua_str(const string name) const;

};


// -------------------------------------------------------------------------
// XPolyline
// A collection of XPath objects.
//
// This echos the behaviour of a Polyline object.
// The main difference is that the XPolyline may be evaluated
// in terms of x.  The segments should be single-valued piece-wise
// continuous segments in terms of x.
//
// Added by: Rowan J. Gollan
// Date: 12-Feb-2008

class XPolyline : public XPath {
public:
    vector<XPath*> seg;
    vector<double> x_seg;
    XPolyline(const vector<XPath*> &segments, const string label="");
    XPolyline(const XPolyline &x );
    virtual ~XPolyline();
    /// Returns a pointer to a (new) cloned path.
    XPolyline* clone() const; 

    // Evaluate the path given x as parameter
    Vector3 xeval(double x) const;
    // Evaluate the first derivative w.r.t. x
    double dydx(double x) const;
    // Evaluate the second derivative w.r.t. x
    double d2ydx2(double x) const;
    string lua_str(const string name) const;
};


#endif
