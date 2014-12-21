/** \file gpath.cxx
 *  \ingroup libgeom2
 *  \brief Implementation of the geometric-path classes. 
 *  \author PJ
 *  \version 31-Dec-2005 initial coding
 *
 * This module is a C++ replacement for the combined C+Python gpath code that
 * was built for mb_cns and Elmer.
 */
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "geom.hh"
#include "gpath.hh"
#include "../../util/source/useful.h"
#include "../../nm/source/secant.hh"
#include "../../nm/source/fobject.hh"
#include "../../nm/source/golden_section_search.hh"
#include "surface.hh"
#include "nurbs_utils.hh"
using namespace std;

//------------------------------------------------------------------------------

// Base class for geometric paths.
// We'll use this later for building Polylines composed of 
// a mix of Path elements.
Path::Path( const string label, double t0, double t1 )
    : label(label), t0(t0), t1(t1) {}
Path::Path( const Path &p )
    : label(p.label), t0(p.t0), t1(p.t1) {}
Path::~Path() {}
Path* Path::clone() const
{
    return new Path(*this);
}
Path* Path::copy(int direction) const
{
    cout << "Path::copy() should not be called." << endl;
    return NULL;
}
Vector3 Path::eval( double t ) const
{ 
    cout << "Path::eval() does nothing." << endl;
    return Vector3(0.0, 0.0, 0.0); 
}
// A global piece of data so that we can
// create a minimum function for the optimiser
// used to locate a point.

// Note this is a minimisation problem because we are
// trying to approach zero.  It is hard to write this
// as a zero-finding problem because we do not cross the
// axis, that is, the scalar length between the guessed
// value and the desired value is *always* a positive number.


static Vector3 p_sought;
static Path *path = 0;
double error_in_point(double t)
{
    return vabs(p_sought - path->eval(t));
}

// Note this function can only find *one* location on
// the path where the path crosses p.  This may not
// be suitable if your path crosses back on itself
// AND you happen to be looking for the particular
// points where the path crosses.
double Path::locate(const Vector3 &p, int &result_flag,
		    double tolerance, int max_iterations) 
{
    // Do this to mimic "closure"
    p_sought = p;
    path = this;

    double t = golden_section_search(error_in_point, 0.0, 1.0, result_flag, tolerance, max_iterations);

    return t;
}



Vector3 Path::dpdt( double t ) const
{
    // Obtain the derivative approximately, via a finite-difference.
    double dt = 0.001;
    Vector3 p0 = eval(t);
    Vector3 derivative = Vector3();
    if ( t+dt > 1.0 ) {
	// t is close to the t=1.0 boundary, use a one-sided difference.
	Vector3 pminus1 = eval(t-dt);
	derivative = (p0 - pminus1) / dt;
    } else if ( t-dt < 0.0 ) {
	// s is close to the s=0 boundary, use a one-sided difference.
	Vector3 pplus1 = eval(t+dt);
	derivative = (pplus1 - p0) / dt;
    } else {
	// Not near a boundary, use central-difference.
	Vector3 pminus1 = eval(t-dt);
	Vector3 pplus1 = eval(t+dt);
	derivative = (pplus1 - pminus1) / (2.0 * dt);
    }
    return derivative;
}
double Path::length() const
{
    double L = 0.0;
    int n = 20;
    double dt = 1.0 / n;
    Vector3 p0 = eval(0.0);
    Vector3 p1;
    for ( int i = 1; i <= n; ++i ) {
	p1 = eval(dt * i);
	L += vabs(p1 - p0);
	p0 = p1;
    }
    return L;
}
double Path::partial_length(double ta, double tb) const
{
    if( tb < ta ) {
	double tmp = ta; ta = tb; tb = tmp;
    }
    double L = 0.0;
    int n = 100;
    double dt = (tb - ta) / n;
    Vector3 p0 = eval(ta);
    Vector3 p1;
    for ( int i = 1; i <= n; ++i ) {
	p1 = eval(ta + dt * i);
	L += vabs(p1 - p0);
	p0 = p1;
    }
    return L;
}
Vector3 Path::point_from_length(double length, double &t) const
{
    double L = 0.0;
    int n = 1000;
    double dt = 1.0 / n;
    Vector3 p0 = eval(0.0);
    Vector3 p1;
    for ( int i = 1; i <= n; ++i ) {
	p1 = eval(dt * i);
	L += vabs(p1 - p0);
	p0 = p1;
	if( L > length ) {
	    t = dt * i;
	    return p1;
	}
    }
    t = dt * n;
    return p1;
	
}
string Path::str() const
{
    ostringstream ost;
    ost << "Path(" << ", \"" << label << "\", "
	<< t0 << ", " << t1 << "\")";
    return ost.str();
}
Path* Path::translate( const Vector3 &v ) 
{
    cout << "Path::translate() does nothing." << endl 
	 << "We should never have called this code." << endl;
    return this;
}
Path* Path::translate( double vx, double vy, double vz ) 
{
    cout << "Path::translate() does nothing." << endl 
	 << "We should never have called this code." << endl;
    return this;
}
Path* Path::reverse() 
{
    double t0_old = t0;
    double t1_old = t1;
    t0 = 1.0 - t1_old;
    t1 = 1.0 - t0_old;
    return this;
}
Path* Path::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    cout << "Path::mirror_image() does nothing." << endl 
	 << "We should never have called this code." << endl;
    return this;
}

Path* Path::rotate_about_zaxis( double dtheta )
{
    cout << "Path::rotate_about_zaxis() does nothing." << endl 
	 << "We should never have called this code." << endl;
    return this;
}

// Overload stream output for Path objects
ostream& operator<<( ostream &os, const Path &p )
{
    os << p.str();
    return os;
}

//------------------------------------------------------------------------------

// Straight line segments.
Line::Line( const Vector3 &a, const Vector3 &b, 
	    const string label, double t0, double t1 )
    : Path(label, t0, t1), a(a), b(b) {}
Line::Line( const Line &line )
    : Path(line.label, line.t0, line.t1), a(line.a), b(line.b) {}
Line::~Line() {}
Line* Line::clone() const
{
    return new Line(*this);
} 
Line* Line::copy(int direction) const
{
    Line* new_path = new Line(*this);
    if (direction==-1) new_path->reverse();
    return new_path;
}
Vector3 Line::eval( double t ) const
{
    t = t0 + t * (t1 - t0); // scale to subrange
    return (1.0 - t) * a + t * b; 
}
double Line::length() const
{
    return fabs(t1 - t0) * vabs(b - a);
}
string Line::str() const
{
    ostringstream ost;
    ost << "Line(" << a << ", " << b << ", \"" 
	<< label << "\", " << t0 << ", " << t1 << ")";
    return ost.str();
}
Line* Line::translate( const Vector3 &v )
{
    a += v; b += v; 
    return this;
}
Line* Line::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
Line* Line::reverse()
{
    Vector3 tmp;
    tmp = a; a = b; b = tmp;
    double t0_old = t0;
    double t1_old = t1;
    t0 = 1.0 - t1_old;
    t1 = 1.0 - t0_old;
    return this;
}
Line* Line::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    a.mirror_image(point, normal);
    b.mirror_image(point, normal);
    return this;
}
Line* Line::rotate_about_zaxis( double dtheta )
{
    a.rotate_about_zaxis(dtheta);
    b.rotate_about_zaxis(dtheta);
    return this;
}

//------------------------------------------------------------------------------

// Circular arcs with start, end and centre points.
Arc::Arc( const Vector3 &a, const Vector3 &b, const Vector3 &c, 
	  const string label, double t0, double t1 )
    : Path(label, t0, t1), a(a), b(b), c(c) {}
Arc::Arc( const Arc &arc )
    : Path(arc.label, arc.t0, arc.t1), a(arc.a), b(arc.b), c(arc.c) {}
Arc::~Arc() {}
Arc* Arc::clone() const
{
    return new Arc(*this);
} 
Arc* Arc::copy(int direction) const
{
    Arc* new_path = new Arc(*this);
    if (direction==-1) new_path->reverse();
    return new_path;
}
Vector3 Arc::eval( double t ) const
{
    Vector3 p;
    t = t0 + t * (t1 - t0); // scale to subrange
    double L;
    if ( evaluate_position_and_length(t, p, L) == 0 ) {
	return p;
    } else {
	return Vector3(0.0, 0.0, 0.0);
    }
}
double Arc::length() const
{
    Vector3 p;
    double L;
    if ( evaluate_position_and_length(1.0, p, L) == 0 ) {
	return fabs(t1 - t0) * L;
    } else {
	return 0.0;
    }
}
string Arc::str() const
{
    ostringstream ost;
    ost << "Arc(" << a << ", " << b << ", " << c 
	<< ", \"" << label << "\", " << t0 << ", " << t1 << ")";
    return ost.str();
}
Arc* Arc::translate( const Vector3 &v )
{
    a += v; b += v; c += v;
    return this;
}
Arc* Arc::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
Arc* Arc::reverse()
{
    Vector3 tmp;
    tmp = a; a = b; b = tmp;
    double t0_old = t0;
    double t1_old = t1;
    t0 = 1.0 - t1_old;
    t1 = 1.0 - t0_old;
    return this;
}
Arc* Arc::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    a.mirror_image(point, normal);
    b.mirror_image(point, normal);
    c.mirror_image(point, normal);
    return this;
}
Arc* Arc::rotate_about_zaxis( double dtheta )
{
    a.rotate_about_zaxis(dtheta);
    b.rotate_about_zaxis(dtheta);
    c.rotate_about_zaxis(dtheta);
    return this;
}
int Arc::evaluate_position_and_length( double t, Vector3& loc, double &L ) const
{
    // Both the position of the point and the length of the full arc are evaluated
    // using mostly the same process of transforming to the plane local to the arc.
    Vector3 ca, cb, tangent1, tangent2, n, cb_local;
    double ca_mag, cb_mag, theta;

    L = 0.0;
    ca = a - c; ca_mag = vabs(ca);
    cb = b - c; cb_mag = vabs(cb);
    // cout << "ca_mag=" << ca_mag << ", cb_mag=" << cb_mag << endl;
    if ( fabs(ca_mag - cb_mag) > 1.0e-5 ) {
	cout << "Arc::eval(): radii do not match ca="
	     << ca << " cb=" << cb << endl;
	return -1;
    }
    // First vector in plane.
    tangent1 = Vector3(ca); tangent1.norm(); 
    // cout << "tangent1=" << tangent1 << endl;
    // Compute unit normal to plane of all three points.
    n = cross(ca, cb); // cout << "n=" << n << endl;
    if ( vabs(n) > 0.0 ) {
	n.norm();
    } else {
	cout << "Arc.eval(): cannot find plane of three points." << endl;
	return -2;
    }
    // Third (orthogonal) vector is in the original plane.
    tangent2 = cross(n, tangent1); 
    // cout << "tangent2=" << tangent2 << endl;
    // Now transform to local coordinates so that we can do 
    // the calculation of the point along the arc in 
    // the local xy-plane, with ca along the x-axis.
    cb_local = cb;
    Vector3 zero = Vector3(0.0,0.0,0.0);
    cb_local.transform_to_local(tangent1, tangent2, n, zero);
    // cout << "cb_local=" << cb_local << endl;
    if ( fabs(cb_local.z) > 1.0e-6 ) {
	cout << "Arc.eval(): problems with transformation cb_local=" 
	     << cb_local << endl;
	return -3;
    }
    // Angle of the final point on the arc is in the range -pi < th <= +pi.
    theta = atan2(cb_local.y, cb_local.x);
    // The length of the circular arc.
    L = theta * cb_mag;
    // cout << "whole arc: theta=" << theta << " L=" << L << endl;
    // Move the second point around the arc in the local xy-plane.
    theta *= t;
    loc.x = cos(theta) * cb_mag;
    loc.y = sin(theta) * cb_mag;
    loc.z = 0.0;
    // cout << "in local plane: loc=" << loc << endl;
    // Transform back to global xyz coordinates
    // and remember to add the centre coordinates.
    loc.transform_to_global(tangent1, tangent2, n, c);
    // cout << "in global coords: loc=" << loc << endl;
    return 0;
}

//------------------------------------------------------------------------------

// Circular arcs with start-, mid- and end-points.

// Helper functions for the Arc3 constructor.
static Vector3 ab_mid;
static Vector3 nab;
static Vector3 global_a;
static Vector3 global_b;
Vector3 locate_centre( double s ) { return ab_mid + (s * nab); }
double error_in_radius( double s )
{
    Vector3 centre = locate_centre(s); // It's somewhere along the bisector.
    return vabs(global_a - centre) - vabs(global_b - centre); // difference in radii
}

Arc3::Arc3( const Vector3 &_a, const Vector3 &_b, const Vector3 &_c, 
	    const string label, double _t0, double _t1 )
    : Arc(_a,_c, Vector3(), label, _t0, _t1) 
{
    // Compute normal vector to the plane of the circle.
    global_a = _a; // start of Arc3
    global_b = _c; // end of Arc3
    Vector3 n = cross(_b - _a, _b - _c);
    if ( vabs(n) <= 1.0e-11 ) {
	cout << "Arc3: Points appear colinear." << endl;
    }
    // The centre of the circle lies along the bisector of ab (and bc).
    ab_mid = 0.5 * (_a + _b);
    nab = cross(_b - _a, n);
    int result_flag;
    double s = secant_solve(error_in_radius, 1.0, 1.01, result_flag);
    c = locate_centre(s);  // member c is ultimately the centre
}
// Already in start-end-centre form.
Arc3::Arc3( const Arc3 &arc )
    : Arc(arc.a, arc.b, arc.c, arc.label, arc.t0, arc.t1) {} 
Arc3::~Arc3() {}

//------------------------------------------------------------------------------
// Helix

Helix::Helix( const Vector3 &point_start, const Vector3 &point_end, 
	      const Vector3 &axis0, const Vector3 &axis1, 
	      const string label, double t0, double t1 )
    : Path(label, t0, t1) 
{
    // Local vectors relative to axis0.
    Vector3 a = axis1 - axis0;
    Vector3 b = point_start - axis0;
    Vector3 c = point_end - axis0;
    zdsh = unit(a);
    Vector3 a0b = b - dot(b,a)*zdsh;
    xdsh = unit(a0b);
    ydsh = cross(zdsh,xdsh);
    a0 = axis0 + dot(b,zdsh)*zdsh;
    a1 = axis0 + dot(c,zdsh)*zdsh;
    r0 = dot(b,xdsh);
    Vector3 a1c = c - dot(c,zdsh)*zdsh;
    r1 = vabs(a1c);
    a1c.transform_to_local(xdsh, ydsh, zdsh, Vector3(0.0, 0.0, 0.0));
    theta01 = atan2(a1c.y, a1c.x);
}
Helix::Helix( const Vector3 &a0, const Vector3 &a1, 
	      const Vector3 &xlocal, double r0, double r1, double dtheta, 
	      const string label, double t0, double t1 )
    : Path(label, t0, t1), a0(a0), a1(a1), xdsh(unit(xlocal)), 
      r0(r0), r1(r1), theta01(dtheta)
{
    // set up local unit vectors at p0
    zdsh = unit(a1 - a0);  // along the axis of the helix
    ydsh = cross(zdsh, xdsh); // normal to both
}

Helix::Helix( const Helix &h )
    : Path(h.label, h.t0, h.t1), a0(h.a0), a1(h.a1), 
      xdsh(h.xdsh), ydsh(h.ydsh), zdsh(h.zdsh),
      r0(h.r0), r1(h.r1), theta01(h.theta01)
{}
Helix::~Helix() {}
Helix* Helix::clone() const
{
    return new Helix(*this);
} 
Helix* Helix::copy(int direction) const
{
    Helix* new_path = new Helix(*this);
    if (direction==-1) new_path->reverse();
    return new_path;
}
Vector3 Helix::eval( double t ) const
{
    t = t0 + t * (t1 - t0); // scale to subrange
    double r = r0 * (1.0-t) + r1 * t;
    double theta = theta01 * t;
    Vector3 p = r * cos(theta) * xdsh + r * sin(theta) * ydsh +
	a0 * (1.0-t) + a1 * t;
    return p;
}
string Helix::str() const
{
    ostringstream ost;
    ost << "Helix(" << a0 << ", " << a1 << ", " << xdsh << ", " 
	<< r0 << ", " << r1 << ", " << theta01 
	<< ", \"" << label << "\", " << t0 << ", " << t1 << ")";
    return ost.str();
}
Helix* Helix::translate( const Vector3 &v )
{
    a0 += v; a1 += v;
    return this;
}
Helix* Helix::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
Helix* Helix::reverse()
{
    // Original end-point relative to original a1.
    Vector3 c = r1 * cos(theta01) * xdsh + r1 * sin(theta01) * ydsh;
    Vector3 vtmp = a0; a0 = a1; a1 = vtmp;
    double tmp = r0; r0 = r1; r1 = tmp;
    // new local axes
    xdsh = unit(c); 
    zdsh = -zdsh;
    ydsh = cross(zdsh, xdsh);
    // Because of the right-hand rule for screws,
    // theta01 remains unchanged.
    double t0_old = t0;
    double t1_old = t1;
    t0 = 1.0 - t1_old;
    t1 = 1.0 - t0_old;
    return this;
}
Helix* Helix::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    cout << "Helix::mirror_image() not implemented." << endl;
    return this;
}
Helix* Helix::rotate_about_zaxis( double dtheta )
{
    cout << "Helix::rotate_about_zaxis() not implemented." << endl;
    return this;
}

//------------------------------------------------------------------------------

// Bezier curves of arbitrary order.
Bezier::Bezier( const vector<Vector3> &B, string label,
		double t0, double t1, int arc_length_p )
    : Path(label, t0, t1), B(B), arc_length_param_flag(arc_length_p)
{
    if ( arc_length_param_flag ) {
	n_arc_length = 100;
    } else {
	n_arc_length = 0;
    }
    set_arc_length_vector();
    set_deriv_control_points();
}
Bezier::Bezier( const Bezier &bez ) 
    : Path(bez.label, bez.t0, bez.t1), B(bez.B),
      arc_length_param_flag(bez.arc_length_param_flag)
{
    n_arc_length = bez.n_arc_length;
    set_arc_length_vector();
    set_deriv_control_points();
}
Bezier::~Bezier() 
{
    arc_length.resize(0);
    C.resize(0);
    D.resize(0);
}
Bezier* Bezier::clone() const
{
    return new Bezier(*this);
} 
Bezier* Bezier::copy(int direction) const
{
    Bezier* new_path = new Bezier(*this);
    if (direction==-1) new_path->reverse();
    return new_path;
}
Bezier* Bezier::add_point( const Vector3 &p )
{
    B.push_back(Vector3(p));
    set_arc_length_vector();
    set_deriv_control_points();
    return this;
}
Vector3 Bezier::raw_eval( double t ) const
// Evaluate B(t) without considering arc_length parameterization flag or subrange.
{
    if ( B.size() == 0 ) return Vector3(0.0, 0.0, 0.0);
    if ( B.size() == 1 ) return B[0];
    size_t n_order = B.size() - 1;
    // Apply de Casteljau's algorithm. 
    vector<Vector3> Q = B; // work array will be overwritten
    for ( size_t k = 0; k < n_order; ++k ) {
	for ( size_t i = 0; i < n_order-k; ++i ) {
	    Q[i] = (1.0 - t) * Q[i] + t * Q[i+1];
	}
    }
    return Q[0];
}
Vector3 Bezier::eval( double t ) const
// Evaluate B(t) considering arc_length parameterization flag and possible subrange of t.
{
    if ( B.size() == 0 ) return Vector3(0.0, 0.0, 0.0);
    if ( B.size() == 1 ) return B[0];
    t = t0 + t * (t1 - t0); // scale to subrange
    if ( arc_length_param_flag ) {
	t = t_from_arc_length(t);
    }
    return raw_eval(t);
}
Vector3 Bezier::dpdt( double t ) const
{
    if ( C.size() == 0 ) return Vector3(0.0, 0.0, 0.0);
    if ( C.size() == 1 ) return C[0];
    t = t0 + t * (t1 - t0); // scale to subrange
    size_t n_order = C.size() - 1;
    // Apply de Casteljau's algorithm. 
    vector<Vector3> Q = C; // work array; will be overwritten
    for ( size_t k = 0; k < n_order; ++k ) {
	for ( size_t i = 0; i < n_order-k; ++i ) {
	    Q[i] = (1.0 - t) * Q[i] + t * Q[i+1];
	}
    }
    // FIX-ME Please Rowan: The subrange selection (if it's not 0->1)
    // and the arc_length parameterization probably messes with your
    // derivative calculations.
    return Q[0];
}
string Bezier::str() const
{
    ostringstream ost;
    ost << "Bezier([";
    for ( size_t i = 0; i < B.size(); ++i ) {
	ost << B[i];
	if (i < B.size()-1 ) ost << ", "; else ost << "], ";
    }
    ost << "\"" << label << "\", " << t0 << ", " << t1 << ", " 
	<< arc_length_param_flag << ")";
    return ost.str();
}
Bezier* Bezier::translate( const Vector3 &v )
{
    for ( size_t i = 0; i < B.size(); ++i ) {
	B[i] += v;
    }
    set_arc_length_vector();
    set_deriv_control_points();
    return this;
}
Bezier* Bezier::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    set_arc_length_vector();
    set_deriv_control_points();
    return this;
}
Bezier* Bezier::reverse()
{
    std::reverse(B.begin(),B.end()); 
    double t0_old = t0;
    double t1_old = t1;
    t0 = 1.0 - t1_old;
    t1 = 1.0 - t0_old;
    set_arc_length_vector();
    set_deriv_control_points();
    return this;
}
Bezier* Bezier::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    for ( size_t i = 0; i < B.size(); ++i ) {
	B[i].mirror_image(point, normal);
    }
    set_arc_length_vector();
    set_deriv_control_points();
    return this;
}
Bezier* Bezier::rotate_about_zaxis( double dtheta ) 
{
    for ( size_t i = 0; i < B.size(); ++i ) {
	B[i].rotate_about_zaxis(dtheta);
    }
    set_arc_length_vector();
    set_deriv_control_points();
    return this;
}
Vector3 Bezier::d2pdt2( double t ) const
{
    if ( D.size() == 0 ) return Vector3(0.0, 0.0, 0.0);
    if ( D.size() == 1 ) return D[0];
    t = t0 + t * (t1 - t0); // scale to subrange
    size_t n_order = D.size() - 1;
    // Apply de Casteljau's algorithm. 
    vector<Vector3> Q = D; // work array will be overwritten
    for ( size_t k = 0; k < n_order; ++k ) {
	for ( size_t i = 0; i < n_order-k; ++i ) {
	    Q[i] = (1.0 - t) * Q[i] + t * Q[i+1];
	}
    }
    // FIX-ME Please Rowan: The subrange selection (if it's not 0->1)
    // and the arc_length parameterization probably messes with your
    // derivative calculations.  Maybe this function should not be
    // called for curves with arc_lenght_param_flag == 1.
    return Q[0];
}
void Bezier::set_deriv_control_points()
{
    size_t n = B.size() - 1;
    if ( n == 0 ) return;
    C.resize(n);
    for (size_t i = 0; i < n; ++i ) C[i] = n*(B[i+1] - B[i]);

    if ( n == 1 ) return;
    D.resize(n-1);
    for (size_t i = 0; i < (n-1); ++i ) D[i] = (n-1)*(C[i+1] - C[i]);
}
void Bezier::set_arc_length_vector()
{
    // Compute the arc_lengths for a number of sample points 
    // so that these can later be used to do a reverse interpolation
    // on the evaluation parameter.
    if ( n_arc_length == 0 ) return;
    double dt = 1.0 / n_arc_length;
    arc_length.resize(0);
    double L = 0.0;
    arc_length.push_back(L);
    Vector3 p0 = raw_eval(0.0);
    Vector3 p1;
    for ( int i = 1; i <= n_arc_length; ++i ) {
	p1 = raw_eval(dt * i);
	L += vabs(p1 - p0);
	arc_length.push_back(L);
	p0 = p1;
    }
    return;
}
double Bezier::t_from_arc_length(double t) const
{
    // The incoming parameter value, t, is proportional to arc_length.
    // Do a reverse look-up from the arc_length to the original t parameter
    // of the Bezier curve.
    double L_target = t * arc_length[n_arc_length];
    // Starting from the right-hand end,
    // let's try to find a point to the left of L_target.
    // If the value is out of range, this should just result in
    // us extrapolating one of the end segments -- that's OK.
    int i = n_arc_length - 1;
    double dt = 1.0 / n_arc_length;
    while ( L_target < arc_length[i] && i > 0 ) i--;
    double frac = (L_target - arc_length[i]) / (arc_length[i+1] - arc_length[i]);
    return (1.0 - frac) * dt*i + frac * dt*(i+1);
 }

//------------------------------------------------------------------------------
Nurbs::Nurbs( const vector<Vector3> &P, const vector<double> w,
	      int p, const vector<double> U,
	      string label, double t0, double t1 )
    : Path(label, t0, t1), P(P), w(w), p(p), U(U)
{
    Pw.resize(P.size());
    umin_ = U.front();
    umax_ = U.back();
    for ( size_t i = 0; i < Pw.size(); ++i ) {
	Pw[i].wx = w[i]*P[i].x;
	Pw[i].wy = w[i]*P[i].y;
	Pw[i].wz = w[i]*P[i].z;
	Pw[i].w = w[i];
    }
}

Nurbs::Nurbs( const vector<Vector3> &P, const vector<double> w,
	      int p, const vector<double> U, const vector<int> inf_list,
	      string label, double t0, double t1 )
    : Path(label, t0, t1), P(P), w(w), p(p), U(U), inf_list(inf_list)
{
    Pw.resize(P.size());
    umin_ = U.front();
    umax_ = U.back();
    for ( size_t i = 0; i < Pw.size(); ++i ) {
	vector<int>::const_iterator result;
	result = find(inf_list.begin(), inf_list.end(), (int) i);
	
	if ( result == inf_list.end() ) {
	    // This point is NOT in the list of infinite control points.
	    // So proceed normally.
	    Pw[i].wx = w[i]*P[i].x;
	    Pw[i].wy = w[i]*P[i].y;
	    Pw[i].wz = w[i]*P[i].z;
	    Pw[i].w = w[i];
	}
	else {
	    // This is in an infinite control point.  Leave values untouched,
	    // but set weight to zero.
	    Pw[i].wx = P[i].x;
	    Pw[i].wy = P[i].y;
	    Pw[i].wz = P[i].z;
	    Pw[i].w = 0.0;
	}
    }
}


Nurbs::Nurbs( const vector<Mapped_point> &Pw, int p,
	      const vector<double> U,
	      string label, double t0, double t1 )
    : Path(label, t0, t1), p(p), U(U), Pw(Pw)
{
    umin_ = U.front();
    umax_ = U.back();
    P.resize(Pw.size());
    w.resize(Pw.size());

    for ( size_t i = 0; i < Pw.size(); ++i ) {
	if ( Pw[i].w == 0.0 ) {
	    // An infinite control point
	    P[i].x = Pw[i].wx;
	    P[i].y = Pw[i].wy;
	    P[i].z = Pw[i].wz;
	    w[i] = 0.0;
	}
	else {
	    w[i] = Pw[i].w;
	    P[i].x = Pw[i].wx/w[i];
	    P[i].y = Pw[i].wy/w[i];
	    P[i].z = Pw[i].wz/w[i];
	}
    }
}

Nurbs::Nurbs( const Nurbs &n )
    : Path(n.label, n.t0, n.t1), P(n.P), w(n.w),
      p(n.p), U(n.U), Pw(n.Pw), inf_list(n.inf_list),
      umin_(n.umin_), umax_(n.umax_) {}

Nurbs::~Nurbs() {}

Nurbs* Nurbs::clone() const
{
    return new Nurbs(*this);
}

Vector3 Nurbs::eval( double t ) const
{
    // Linear re-parameterization from t --> u.
    double u = map_t2u(t);
    return nurbs_curve_point(u, p, U, Pw);
}

string Nurbs::str() const
{
    ostringstream ost;
    ost << "Nurbs([ \\" << endl;
    for ( size_t i = 0; i < P.size(); ++i ) {
	ost << P[i];
	if ( i < P.size()-1 ) ost << ", "; else ost << "], \\ \n";
    }
    for ( size_t i = 0; i < w.size(); ++i ) {
	ost << w[i];
	if ( i < w.size()-1 ) ost << ", "; else ost << "], \\ \n";
    }
    ost << p << ", \\ \n";
    for ( size_t i = 0; i < U.size(); ++i ) {
	ost << U[i];
	if ( i < U.size()-1 ) ost << ", "; else ost << "], \\ \n";
    }
    ost << "\"" << label << "\", " << t0 << ", " << t1 << ")";
    return ost.str();
}

Vector3
Nurbs::
dpdu(double t) const
{
    double u = map_t2u(t);

    int n = Pw.size() - 1;
    int d = 1;
    vector<Mapped_point> CKw;
    CKw.resize(d+1);

    curve_derivs(n, p, U, Pw, u, d, CKw);

    vector<Vector3> CK;
    CK.resize(d+1);

    rat_curve_derivs(CKw, d, CK);

    return CK[1];
}

Vector3
Nurbs::
d2pdu2(double t) const
{
    double u = map_t2u(t);

    int n = Pw.size() - 1;
    int d = 2;
    vector<Mapped_point> CKw;
    CKw.resize(d+1);

    curve_derivs(n, p, U, Pw, u, d, CKw);

    vector<Vector3> CK;
    CK.resize(d+1);

    rat_curve_derivs(CKw, d, CK);

    return CK[2];
}

double
Nurbs::
map_t2u(double t) const
{
    return umin_ + t*(umax_ - umin_);
}

Nurbs
Nurbs::
knot_insertion(double u, int k, int s, int r) const
{
    int np = P.size()-1;
    int nq;
    vector<double> UQ;
    vector<Mapped_point> Qw;

    if ( curve_knot_ins(np, p, U, Pw,
			u, k, s, r, nq,
			UQ, Qw) != SUCCESS ) {
	cout << "Nurbs::knot_insertion():\n";
	cout << "Problem inserting knot value: " << u << endl;
	cout << "Bailing out!\n";
	exit(1);
    }

    return Nurbs(Qw, p, UQ, label);

}

Nurbs
Nurbs::
knot_refinement(const vector<double> &X) const
{
    if ( X.size() == 0 ) {
	return Nurbs(*this);
    }
    int n = P.size()-1;
    int r = X.size()-1;
    vector<double> Ubar;
    vector<Mapped_point> Qw;

    if ( refine_knot_vector_curve(n, p, U, Pw, X,
				  r, Ubar, Qw) != SUCCESS ) {
	cout << "Nurbs::knot_refinement():\n";
	cout << "Problem refining knot vector.\n";
	cout << "Bailing out!\n";
	exit(1);
    }

    return Nurbs(Qw, p, Ubar, label);
}

//------------------------------------------------------------------------------

// Polyline curves consist of one or more Path segments.
Polyline::Polyline( const vector<Path*> &segments, 
		    string label, double t0, double t1, int arc_length_p )
    : Path(label, t0, t1), seg(segments),
      t_seg(vector<double>(segments.size())),
      arc_length_param_flag(arc_length_p)
{
    // Replace the original pointers with pointers to cloned objects
    // so that we don't have to worry about their persistence.
    for ( size_t i = 0; i < seg.size(); ++i ) {
	seg[i] = segments[i]->clone();
    }
    if ( arc_length_param_flag ) {
	n_arc_length = 100;
    } else {
	n_arc_length = 0;
    }
    reset_breakpoints();
    set_arc_length_vector();
}
Polyline::Polyline( const vector<Vector3*> &points, int type_of_segments,
		    string label, int arc_length_p )
    : Path(label, 0.0, 1.0), 
      seg(vector<Path*>(points.size()-1)), 
      t_seg(vector<double>(points.size()-1) ),
      arc_length_param_flag(arc_length_p)
{
    // For the moment, just implement straight-line segments.
    for ( size_t i = 0; i < points.size()-1; ++i ) {
	seg[i] = new Line(*points[i], *points[i+1]);
    }
    reset_breakpoints();
    set_arc_length_vector();
}
Polyline::Polyline( int nseg, string label, double t0, double t1, int arc_length_p )
    : Path(label, t0, t1), seg(vector<Path*>(nseg)),
      t_seg(vector<double>(nseg)),
      arc_length_param_flag(arc_length_p)
{}
Polyline::Polyline( const Polyline &pline ) 
    : Path(pline.label, pline.t0, pline.t1), seg(pline.seg),
      t_seg(pline.t_seg),
      arc_length_param_flag(pline.arc_length_param_flag)
{
    // Replace the original pointers with pointers to cloned objects
    // so that we don't have to worry about their persistence.
    for ( size_t i = 0; i < seg.size(); ++i ) {
	seg[i] = pline.seg[i]->clone();
    }
    n_arc_length = pline.n_arc_length;
    reset_breakpoints();
    set_arc_length_vector();
}
Polyline::~Polyline()
{
    arc_length.resize(0);
    // Destroy the cloned objects.
    for ( size_t i = 0; i < seg.size(); ++i ) {
	delete seg[i];
    }
}
Polyline* Polyline::clone() const
{
    return new Polyline(*this);
} 
Polyline* Polyline::copy(int direction) const
{
    Polyline* new_path = new Polyline(*this);
    if (direction==-1) new_path->reverse();
    return new_path;
}
Polyline* Polyline::add_segment( const Path* segment, int direction )
{
    Path* new_segment = segment->clone();
    if ( direction == -1 ) new_segment->reverse();
    seg.push_back(new_segment);
    reset_breakpoints();
    set_arc_length_vector();
    return this;
}
Vector3 Polyline::raw_eval( double t ) const
// Evaluate B(t) without considering arc_length parameterization flag or subrange.
{
    double t_local;
    size_t n = seg.size();
    if ( n == 0 ) return Vector3(0.0, 0.0, 0.0);
    if ( n == 1 ) return seg[0]->eval(t);
    // Search the segments until we reach the appropriate range of t.
    size_t i;
    for ( i = 0; i < n; ++i ) {
	if ( t <= t_seg[i] ) break;
    }
    if ( i >= n ) i = n - 1;  // last segment
    // At this point, t_seg[i-1] < t <= t_seg[i] (we hope)
    // Have assumed that the t breakpoints are well behaved.
    if ( i == 0 ) {
	t_local = t / t_seg[i];
    } else {
	t_local = (t - t_seg[i-1]) / (t_seg[i] - t_seg[i-1]);
    }
    return seg[i]->eval(t_local);
}
Vector3 Polyline::eval( double t ) const
{
    // First, map to the subrange of the full polyline.
    t = t0 + t * (t1 - t0);
    if ( arc_length_param_flag ) {
	t = t_from_arc_length(t);
    }
    return raw_eval(t);
}
double Polyline::length() const
{
    double L_total = 0.0;
    for ( size_t i = 0; i < seg.size(); ++i ) {
	L_total += seg[i]->length();
    }
    return L_total; 
}
string Polyline::str() const
{
    ostringstream ost;
    ost << "Polyline([";
    for ( size_t i = 0; i < seg.size(); ++i ) {
	ost << *seg[i];
	if (i < seg.size()-1 ) ost << ", "; else ost << "], ";
    }
    ost << "\"" << label << "\", ";
    ost << t0 << ", " << t1 << ", " 
	<< arc_length_param_flag << ")";
    return ost.str();
}
Polyline* Polyline::translate( const Vector3 &v )
{
    for ( size_t i = 0; i < seg.size(); ++i ) {
	seg[i]->translate(v);
    }
    return this;
}
Polyline* Polyline::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
Polyline* Polyline::reverse()
{
    // First, reverse the individual segments, then the whole collection.
    vector<Path*> new_seg;
    for ( size_t i = 0; i < seg.size(); ++i ) {
	Path* seg_i = seg[i]->clone();
	seg_i->reverse();
	new_seg.push_back(seg_i);
    }
    std::reverse(new_seg.begin(),new_seg.end());
    seg = new_seg;
    reset_breakpoints();
    double t0_old = t0;
    double t1_old = t1;
    t0 = 1.0 - t1_old;
    t1 = 1.0 - t0_old;
    set_arc_length_vector();
    return this;
}
Polyline* Polyline::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    for ( size_t i = 0; i < seg.size(); ++i ) {
	seg[i]->mirror_image(point, normal);
    }
    set_arc_length_vector();
    return this;
}
Polyline* Polyline::rotate_about_zaxis( double dtheta ) 
{
    for ( size_t i = 0; i < seg.size(); ++i ) {
	seg[i]->rotate_about_zaxis(dtheta);
    }
    set_arc_length_vector();
    return this;
}
void Polyline::reset_breakpoints()
{
    // Set up the parameter breakpoints based on cumulative length.
    t_seg[0] = seg[0]->length();
    for ( size_t i = 1; i < seg.size(); ++i ) {
	t_seg[i] = t_seg[i-1] + seg[i]->length();
    }
    double L_total = t_seg[t_seg.size()-1];
    for ( size_t i = 0; i < seg.size(); ++i ) {
	t_seg[i] /= L_total;
    }
    return;
}
void Polyline::set_arc_length_vector()
{
    // Compute the arc_lengths for a number of sample points 
    // so that these can later be used to do a reverse interpolation
    // on the evaluation parameter.
    if ( n_arc_length == 0 ) return;
    double dt = 1.0 / n_arc_length;
    arc_length.resize(0);
    double L = 0.0;
    arc_length.push_back(L);
    Vector3 p0 = raw_eval(0.0);
    Vector3 p1;
    for ( int i = 1; i <= n_arc_length; ++i ) {
	p1 = raw_eval(dt * i);
	L += vabs(p1 - p0);
	arc_length.push_back(L);
	p0 = p1;
    }
    return;
}
double Polyline::t_from_arc_length(double t) const
{
    // The incoming parameter value, t, is proportional to arc_length.
    // Do a reverse look-up from the arc_length to the original t parameter
    // of the Bezier curve.
    double L_target = t * arc_length[n_arc_length];
    // Starting from the right-hand end,
    // let's try to find a point to the left of L_target.
    // If the value is out of range, this should just result in
    // us extrapolating one of the end segments -- that's OK.
    int i = n_arc_length - 1;
    double dt = 1.0 / n_arc_length;
    while ( L_target < arc_length[i] && i > 0 ) i--;
    double frac = (L_target - arc_length[i]) / (arc_length[i+1] - arc_length[i]);
    return (1.0 - frac) * dt*i + frac * dt*(i+1);
 }

//------------------------------------------------------------------------------

// A Spline curve is a specialized Polyline.
Spline::Spline( const vector<Vector3> &p, string label, 
		double t0, double t1,
		int arc_length_p, double tolerance )
    : Polyline( p.size()-1, label, t0, t1, arc_length_p )
{
    size_t m = p.size() - 1;
    // Given m+1 interpolation points p, determine the m-segment
    // Bezier polyline that interpolates these points as a spline. 
    // This is done by first determining the array of weight points
    // which define the spline and then evaluating the cubic 
    // Bezier segments.
    // Reference:
    //     G. Engelin & F. Uhlig (1996)
    //     Numerical Algorithms with C
    //     Springer, Berlin
    //     Section 12.3.1
 
    vector<Vector3> d(m+1);  // weight points
    // For a natural spline, the first and last weight points
    // are also the first and last interpolation points.
    d[0] = p[0];
    d[m] = p[m];

    // For the initial guess at the remaining weight points,
    // just use the supplied data points.
    for ( size_t i = 1; i < m; i++ ) {
	d[i] = p[i];
    }
    // Apply Gauss-Seidel iteration until
    // the internal weight points converge.
    int j;
    Vector3 old_p;
    double max_diff;
    for ( j = 1; j < 50; j++ ) {
	max_diff = 0.0;
	for ( size_t i = 1; i < m; i++ ) {
	    old_p = d[i];
	    d[i] = 0.25 * (6.0 * p[i] - d[i-1] - d[i+1]);
	    double diff = vabs(d[i] - old_p);
	    if ( diff > max_diff ) max_diff = diff;
	} /* end for i */
	if ( max_diff < tolerance ) break;
    } /* end for j */

    // Final stage; calculate the Bezier segments and pack them away.
    vector<Vector3> p03(4);
    for ( size_t i = 0; i < m; i++ ) {
	p03[0] = p[i];
	p03[1] = (2.0 * d[i] + d[i+1]) / 3.0;
	p03[2] = (d[i] + 2.0 * d[i+1]) / 3.0;
	p03[3] = p[i+1];
	seg[i] = new Bezier(p03);
    }
    reset_breakpoints();
    if ( arc_length_param_flag ) {
	n_arc_length = 100;
    } else {
	n_arc_length = 0;
    }
     set_arc_length_vector();
}
Spline::Spline( const Spline &spl ) 
    : Polyline(spl.seg, spl.label, spl.t0, spl.t1, spl.arc_length_param_flag) {}

Vector3 Spline::eval_from_x(double x) const
{
    // Implement golden section search inline.
    double eps = 1.0e-8;
    double a = 0.0;
    double c = 1.0;
    double g = (3 - sqrt(5.0))/2.0;
    double b1 = a + g*(c - a);
    double b2 = a + (1 - g)*(c - a);
    double F1 = Fx(b1, x);
    double F2 = Fx(b2, x);

    while( (c - a) > eps*(a+c) ) {
	if( F1 <= F2 ) {
	    c = b2;
	    b2 = b1;
	    F2 = F1;
	    b1 = a + g*(c - a);
	    F1 = Fx(b1, x);
	}
	else {
	    a = b1;
	    b1 = b2;
	    F1 = F2;
	    b2 = a + (1 - g)*(c - a);
	    F2 = Fx(b2, x);
	}
    }

    double t = (a+c)/2.0;
    return eval(t);
}

double Spline::Fx(double t, double x) const
{
    double dx = x - eval(t).x;
    return (dx*dx);
}

Vector3 Spline::eval_from_y(double y) const
{
    // Implement golden section search inline.
    double eps = 1.0e-8;
    double a = 0.0;
    double c = 1.0;
    double g = (3 - sqrt(5.0))/2.0;
    double b1 = a + g*(c - a);
    double b2 = a + (1 - g)*(c - a);
    double F1 = Fy(b1, y);
    double F2 = Fy(b2, y);

    while( (c - a) > eps*(a+c) ) {
	if( F1 <= F2 ) {
	    c = b2;
	    b2 = b1;
	    F2 = F1;
	    b1 = a + g*(c - a);
	    F1 = Fy(b1, y);
	}
	else {
	    a = b1;
	    b1 = b2;
	    F1 = F2;
	    b2 = a + (1 - g)*(c - a);
	    F2 = Fy(b2, y);
	}
    }

    double t = (a+c)/2.0;
    return eval(t);
}

double Spline::Fy(double t, double y) const
{
    double dy = y - eval(t).y;
    return (dy*dy);
}


//------------------------------------------------------------------------------

PathOnSurface::PathOnSurface( const ParametricSurface &_surf, 
			      const UnivariateFunction &_fr,
			      const UnivariateFunction &_fs,
			      string label, double t0, double t1 )
    : Path(label, t0, t1), surf(_surf.clone()), fr(_fr.clone()), fs(_fs.clone()) {}
PathOnSurface::PathOnSurface( const PathOnSurface &p ) 
    : Path(p.label, p.t0, p.t1), surf(p.surf->clone()), 
      fr(p.fr->clone()), fs(p.fs->clone()) {}
PathOnSurface::~PathOnSurface()
{
    delete fs;
    delete fr;
    delete surf;
}
PathOnSurface* PathOnSurface::clone() const
{
    return new PathOnSurface(*this);
} 
Vector3 PathOnSurface::eval( double t ) const
{
    t = t0 + t * (t1 - t0);
    return surf->eval(fr->eval(t), fs->eval(t));
}
string PathOnSurface::str() const
{
    ostringstream ost;
    ost << "PathOnSurface(";
    ost << surf->str() << ", " << fr->str() << ", " << fs->str();
    ost << ", \"" << label << "\", " << t0 << ", " << t1 << ")";
    return ost.str();
}
PathOnSurface* PathOnSurface::translate( const Vector3 &v )
{
    surf->translate(v);
    return this;
}
PathOnSurface* PathOnSurface::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
// PathOnSurface* PathOnSurface::reverse()
// {
//     fr->reverse();
//     fs->reverse();
//     return this;
// }
// PathOnSurface* PathOnSurface::mirror_image( const Vector3 &point, const Vector3 &normal ) 
// {
//     surf->mirror_image(point, normal);
//     return this;
// }

//------------------------------------------------------------------------------
// PolarPath

PolarPath::PolarPath( const Path &pth, double H, 
		      const string label, double t0, double t1 )
    : Path(label, t0, t1), 
      original_path(pth.clone()), H(H) 
{}
PolarPath::PolarPath( const PolarPath &ppth )
    : Path(ppth.label, ppth.t0, ppth.t1), 
      original_path(ppth.original_path->clone()), H(ppth.H) 
{}
PolarPath::~PolarPath() 
{
    delete original_path;
}
PolarPath* PolarPath::clone() const
{
    return new PolarPath(*this);
} 
Vector3 PolarPath::eval( double t ) const
{
    Vector3 p = original_path->eval(t);
    map_neutral_plane_to_cylinder(p, H);
    return p;
}
string PolarPath::str() const
{
    ostringstream ost;
    ost << "PolarPath(" << original_path->str() << ", " << H
	<< ", \"" << label << "\", " << t0 << ", " << t1 << ")";
    return ost.str();
}
PolarPath* PolarPath::translate( const Vector3 &v )
{
    original_path->translate(v);
    return this;
}
PolarPath* PolarPath::translate( double vx, double vy, double vz )
{
    translate(Vector3(vx, vy, vz));
    return this;
}
PolarPath* PolarPath::reverse()
{
    original_path->reverse();
    return this;
}
PolarPath* PolarPath::mirror_image( const Vector3 &point, const Vector3 &normal ) 
{
    // Not really sure that this function will be well behaved!
    original_path->mirror_image(point, normal);
    return this;
}


// ----------------------------------------------------------------------------
// XPath

XPath::
XPath()
    : Path("", 0.0, 1.0), x0(0.0), x1(1.0) {}

XPath::
XPath(double x0, double x1, const string label)
    : Path(label, 0.0, 1.0), x0(x0), x1(x1) {}

XPath::
XPath(const XPath &x)
    : Path(x.label, 0.0, 1.0), x0(x.x0), x1(x.x1) {}

XPath::
~XPath() {}

Vector3
XPath::
eval( double t ) const
{
    double x = map_t2x(t);
    return xeval(x);
}

Vector3
XPath::
dpdt( double t ) const
{
    double x = map_t2x(t);
    double ydash = dydx(x);

    // Scale derivatives appropriately
    double dpdx = x1 - x0;
    double dpdy = (x1 - x0)*ydash;

    return Vector3(dpdx, dpdy, 0.0);

}

// A global piece of data so that we can
// create a function to pass to the secant method.
static XPath *gpath;
static double y_desired;

double error_in_y(double x)
{
    return y_desired - gpath->xeval(x).y;
}

double 
XPath::
locate_x(double y, int &result_flag)
{
    const double tolerance = 1.0e-6;

    // Initialise spline gspline.
    gpath = this;

    // Initialise dydx_desired.
    y_desired = y;

    // Use the secant method to find answer.
    double x = secant_solve(error_in_y, x0, x1, result_flag, tolerance);
    // Finished with gpath.
    return x;
}

//-----------------------------------------------------------------------------
// XPoly

XPoly::
XPoly()
    : XPath(0.0, 1.0, "") {}

XPoly::
XPoly(double x0, double x1,
      vector<double> &B, const string label)
    : XPath(x0, x1, label), B(B)
{
    // Store first derivative as a polynomial
    // for easy evaluation of second derivative
    for( size_t i = 1; i < B.size(); ++i ) {
	D.push_back(i*B[i]);
    }
}

XPoly::
XPoly( const XPoly &x )
    : XPath(x.x0, x.x1, x.label), B(x.B), D(x.D) {}

XPoly::
~XPoly() {}

XPoly&
XPoly::
operator=(const XPoly &x)
{
    if( this == &x ) // avoid aliasing
	return *this;
    else {
	label = x.label;
	t0 = x.t0;
	t1 = x.t1;
	x0 = x.x0;
	x1 = x.x1;
	B = x.B;
	D = x.D;
    }
    return (*this);
}

XPoly*
XPoly::
clone() const
{
    return new XPoly(*this);
}

// Vector3
// XPoly::
// eval( double t ) const
// {
//     double x = map_t2x(t);
//     return xeval(x);
// }

// Vector3
// XPoly::
// dpdt( double t ) const
// {
//     double x = map_t2x(t);
//     double ydash = dydx(x);

//     // Scale derivatives appropriately
//     double dpdx = x1 - x0;
//     double dpdy = (x1 - x0)*ydash;

//     return Vector3(dpdx, dpdy, 0.0);

// }

string
XPoly::
str() const
{
    ostringstream ost;
    ost << "XPoly(" << ", \"" << label << "\", "
	<< x0 << ", " << x1;
    for( size_t i = 0; i < B.size(); ++i ) ost << ", " << B[i];
    ost << "\")";
    return ost.str();
}    

Vector3
XPoly::
xeval(double x) const
{
    // Using Horner's algorithm
    double y = B.back();
    for( int i = int(B.size() - 2); i >= 0; --i ) {
	y = x*y + B[i];
    }
    return Vector3(x, y, 0.0);
}

double
XPoly::
dydx(double x) const
{
    // Using Horner's algorithm
    double y = B.back();
    double dydx = 0.0;

    for( int i = int(B.size() - 2); i >= 0; --i ) {
	dydx = x*dydx + y;
	y = x*y + B[i];
    }
    return dydx;
}

double
XPoly::
d2ydx2(double x) const
{
    if( D.size() == 0.0 )
	return 0.0;

    // Using Horner's algorithm
    double y = D.back();
    double dydx = 0.0;
    for( int i = int(D.size() - 2); i >= 0; --i ) {
	dydx = x*dydx + y;
	y = x*y + D[i];
    }
    return dydx;
}

string
XPoly::
lua_str(const string name) const
{
    ostringstream ost;
    ost << "tmpd = vectord(" << B.size() << ")\n";
    for( size_t i = 0; i < B.size(); ++i ) {
	ost << "tmpd[" << i << "] = " << B[i] << endl;
    }
    ost << name << " = XPoly(" << x0 << ", " << x1 << ", tmpd)\n";
    return ost.str();
}

//------------------------------------------------------------------------------
// XSpline

XSpline::
XSpline(double x0, double x1, vector<XPoly> &X,
	const string label)
    : XPath(x0, x1, label), X(X)
{
    for(size_t i = 0; i < X.size(); ++i ) {
	x_seg.push_back(X[i].eval(1.0).x);
    }
}

XSpline::
XSpline(const XSpline &x)
    : XPath(x.x0, x.x1, x.label),
      X(x.X), x_seg(x.x_seg) {}

XSpline::
~XSpline() {}

XSpline*
XSpline::
clone() const
{
    return new XSpline(*this);
}

// Vector3
// XSpline::
// eval(double t) const
// {
//     double x = map_t2x(t);
//     return xeval(x);
// }

// Vector3
// XSpline::
// dpdt(double t) const
// {
//     double x = map_t2x(t);
//     double ydash = dydx(x);

//     // Scale derivatives appropriately
//     double dpdx = x1 - x0;
//     double dpdy = (x1 - x0)*ydash;

//     return Vector3(dpdx, dpdy, 0.0);

// }


Vector3
XSpline::
xeval(double x) const
{
    // Search the segments until we reach the appropriate range of x.
    size_t n = X.size();
    size_t i;
    for ( i = 0; i < n; ++i ) {
	if ( x <= x_seg[i] ) break;
    }
    if ( i >= n ) i = n - 1;  // last segment
    // At this point, x_seg[i-1] < x <= x_seg[i] (we hope)
    // Have assumed that the x breakpoints are well behaved.
    // Ready to evaluate.
    return X[i].xeval(x);
}

double
XSpline::
dydx(double x) const
{
    // Search the segments until we reach the appropriate range of x.
    size_t n = X.size();
    size_t i;
    for ( i = 0; i < n; ++i ) {
	if ( x <= x_seg[i] ) break;
    }
    if ( i >= n ) i = n - 1;  // last segment
    // At this point, x_seg[i-1] < x <= x_seg[i] (we hope)
    // Have assumed that the t breakpoints are well behaved.
    // Ready to evaluate.
    return X[i].dydx(x);
}

double
XSpline::
d2ydx2(double x) const
{
    // Search the segments until we reach the appropriate range of x.
    size_t n = X.size();
    size_t i;
    for ( i = 0; i < n; ++i ) {
	if ( x <= x_seg[i] ) break;
    }
    if ( i >= n ) i = n - 1;  // last segment
    // At this point, x_seg[i-1] < x <= x_seg[i] (we hope)
    // Have assumed that the t breakpoints are well behaved.
    // Ready to evaluate.
    return X[i].d2ydx2(x);
}

// A global piece of data so that we can
// create a function to pass to the secant method.
static XSpline *gspline;
static double dydx_desired;

double error_in_dydx(double x)
{
    return dydx_desired - gspline->dydx(x);
}

Vector3
XSpline::
locate_dydx(double dydx, int &result_flag) const
{
    const double tolerance = 1.0e-6;
    // Initialise spline gspline.
    delete gspline;
    gspline = clone();

    // Initialise dydx_desired.
    dydx_desired = dydx;

    // Use the secant method to find answer.
    double x = secant_solve(error_in_dydx, x0, x1, result_flag, tolerance);
    // Finished with gspline.
    delete gspline;
    return xeval(x);
}

string
XSpline::
lua_str(const string name) const
{
    ostringstream ost;
    ostringstream tmp_name;
    ost << "tmpXPoly = XPolyList(" << X.size() << ")\n";
    for( size_t i = 0; i < X.size(); ++i ) {
	tmp_name.str("");
	tmp_name << name << "_" << i;
	ost << X[i].lua_str(tmp_name.str());
	ost << "tmpXPoly[" << i << "] = " << tmp_name.str() << endl;
    }
    
    ost << name << " = XSpline(" << x0 << ", " << x1 << ", tmpXPoly)\n";
    return ost.str();

}

// --------------------------------------------------------
// XBezier

XBezier::
XBezier( double x0, double x1, const vector<Vector3> &B,
	 const string label )
    : XPath(x0, x1, label), bez(B, label)
{
    // The re-parameterization in terms of x
    // only makes sense if the t --> x mapping is
    // linear.  This will occur when the x values
    // for the control points are equidistant.
    // Check this is true.  If not, consider just
    // using a normal Bezier.
    const double zero_tol = 1.0e-9;

    if( bez.B.size() >= 3 ) {
	double dx = bez.B[1].x - bez.B[0].x;
	for( size_t i = 2; i < bez.B.size()-1; ++i ) {
	    double test_dx = B[i].x - B[i-1].x;
	    if( fabs(dx - test_dx) > zero_tol ) {
		throw invalid_argument("The control points given to XBezier are NOT equidistant in x.");
	    }
	}
    }
}

XBezier::
XBezier( const XBezier &x )
    : XPath(x.x0, x.x1, x.label), bez(x.bez) {}

XBezier::
~XBezier() {}

XBezier*
XBezier::
clone() const
{
    return new XBezier(*this);
}

Vector3
XBezier::
eval( double t ) const
{
    return bez.eval(t);
}

Vector3
XBezier::
dpdt( double t ) const
{
    return bez.dpdt(t);
}

string
XBezier::
str() const
{
    ostringstream ost;
    ost << "XBezier([";
    for ( size_t i = 0; i < bez.B.size(); ++i ) {
	ost << bez.B[i];
	if (i < bez.B.size()-1 ) ost << ", "; else ost << "], ";
    }
    ost << "\"" << label << "\", " << x0 << ", " << x1 << ")";
    return ost.str();
}

Vector3
XBezier::
xeval(double x) const
{
    return bez.eval(map_x2t(x));
}

double
XBezier::
dydx(double x) const
{
    return (1.0/(x1 - x0))*bez.dpdt(map_x2t(x)).y;
}

double
XBezier::
d2ydx2(double x) const
{
    return (1.0/(x1 - x0))*bez.d2pdt2(map_x2t(x)).y;
}

string
XBezier::
lua_str(const string name) const
{
    ostringstream ost;
    ost << "tmpd = Vector3List(" << bez.B.size() << ")\n";
    for( size_t i = 0; i < bez.B.size(); ++i ) {
	ost << "tmpd[" << i << "] = " << bez.B[i] << endl;
    }
    ost << name << " = XBezier(" << x0 << ", " << x1 << ", tmpd)\n";
    return ost.str();
}

// --------------------------------------------------------
// XPolyline

XPolyline::
XPolyline(const vector<XPath*> &segments,
	  const string label)
    : XPath(0.0, 1.0, label)
{
    for(size_t i = 0; i < segments.size(); ++i ) {
	seg.push_back(segments[i]->clone());
	x_seg.push_back(segments[i]->eval(1.0).x);
    }
    x0 = seg[0]->eval(0.0).x;
    x1 = seg.back()->eval(1.0).x;
}

XPolyline::
XPolyline(const XPolyline &x)
    : XPath(x.x0, x.x1, label), x_seg(x.x_seg)
{
    for( size_t i = 0; i < x.seg.size(); ++i ) {
	seg.push_back(x.seg[i]->clone());
    }
}

XPolyline::
~XPolyline()
{
    for(size_t i = 0; i < seg.size(); ++i ) {
	delete seg[i];
    }
}

XPolyline*
XPolyline::
clone() const
{
    return new XPolyline(*this);
}

Vector3
XPolyline::
xeval(double x) const
{
    
    // Search the segments until we reach the appropriate range of x.
    size_t n = seg.size();
    size_t i;
    for ( i = 0; i < n; ++i ) {
	if ( x <= x_seg[i] ) break;
    }
    if ( i >= n ) i = n - 1;  // last segment
    // At this point, x_seg[i-1] < x <= x_seg[i] (we hope)
    // Have assumed that the x breakpoints are well behaved.
    // Ready to evaluate.
    return seg[i]->xeval(x);
}

double
XPolyline::
dydx(double x) const
{
    // Search the segments until we reach the appropriate range of x.
    size_t n = seg.size();
    size_t i;
    for ( i = 0; i < n; ++i ) {
	if ( x <= x_seg[i] ) break;
    }
    if ( i >= n ) i = n - 1;  // last segment
    // At this point, x_seg[i-1] < x <= x_seg[i] (we hope)
    // Have assumed that the t breakpoints are well behaved.
    // Ready to evaluate.
    return seg[i]->dydx(x);
}

double
XPolyline::
d2ydx2(double x) const
{
    // Search the segments until we reach the appropriate range of x.
    size_t n = seg.size();
    size_t i;
    for ( i = 0; i < n; ++i ) {
	if ( x <= x_seg[i] ) break;
    }
    if ( i >= n ) i = n - 1;  // last segment
    // At this point, x_seg[i-1] < x <= x_seg[i] (we hope)
    // Have assumed that the t breakpoints are well behaved.
    // Ready to evaluate.
    return seg[i]->d2ydx2(x);
}

string
XPolyline::
lua_str(const string name) const
{
    ostringstream ost;
    ostringstream tmp_name;
    ost << "tmpXPath = XPathList(" << seg.size() << ")\n";
    for( size_t i = 0; i < seg.size(); ++i ) {
	tmp_name.str("");
	tmp_name << name << "_" << i;
	ost << seg[i]->lua_str(tmp_name.str());
	ost << "tmpXPath[" << i << "] = " << tmp_name.str() << endl;
    }
    
    ost << name << " = XPolyline(tmpXPath)\n";
    return ost.str();

}
