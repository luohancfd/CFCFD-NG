/** gpath.d
 * Geometry-building elements for our 3D world -- one-parameter elements.
 * Note that these are geometric paths, as distinct from the file-system paths.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 */

module gpath;

import std.conv;
import std.math;
import geom;


class Path {
public:
    double t0; // to subrange t, when evaluating a point on the Path
    double t1;
    abstract Path dup() const;
    abstract Vector3 opCall(double t) const;
    abstract override string toString() const;
} // end class Path


class Line : Path {
public:
    Vector3 p0; // end-point at t=0
    Vector3 p1; // end-point at t=1
    this(in Vector3 p0, in Vector3 p1, double t0=0.0, double t1=1.0)
    {
	this.p0 = p0; this.p1 = p1;
	this.t0 = t0; this.t1 = t1;
    }
    this(ref const(Line) other)
    {
	p0 = other.p0; p1 = other.p1;
	t0 = other.t0; t1 = other.t1;
    }
    override Line dup() const
    {
	return new Line(p0, p1, t0, t1);
    }
    override Vector3 opCall(double t) const 
    {
	double tdsh = t0 + (t1-t0)*t; // subrange t
	return (1.0-tdsh)*p0 + tdsh*p1;
    }
    override string toString() const
    {
	return "Line(p0=" ~ to!string(p0) ~ ", p1=" ~ to!string(p1) ~
	    ", t0=" ~ to!string(t0) ~ ", t1=" ~ to!string(t1) ~ ")";
    }
} // end class Line


unittest {
    auto a = Vector3([1.0, 2.2, 3.0]);
    auto b = Vector3(1.0);
    auto ab = new Line(a, b);
    auto c = ab(0.5);
    assert(approxEqualVectors(c, Vector3(1.0, 1.1, 1.5)), "Line");
    auto ab2 = ab.dup();
    auto d = ab2(0.5);
    assert(approxEqualVectors(c, d), "Line.dup");
}


class Arc : Path {
public:
    Vector3 a; // beginning point (t = 0)
    Vector3 b; // end point (t = 1)
    Vector3 c; // centre of curvature
    // Arc constructed from start-point a, end-point b and centre of curvature c.
    this(in Vector3 a, in Vector3 b, in Vector3 c, double t0=0.0, double t1=1.0)
    {
	this.a = a; this.b = b; this.c = c;
	this.t0 = t0; this.t1 = t1;
    }
    this(ref const(Arc) other)
    {
	a = other.a; b = other.b; c = other.c;
	t0 = other.t0; t1 = other.t1;
    }
    override Arc dup() const
    {
	return new Arc(a, b, c, t0, t1);
    }
    override Vector3 opCall(double t) const 
    {
	double tdsh = t0 + (t1-t0)*t; // subrange t
	double L;
	Vector3 p;
	evaluate_position_and_length(t, p, L);
	return p;
    }
    override string toString() const
    {
	return "Arc(a=" ~ to!string(a) ~ ", b=" ~ to!string(b) ~ ", c=" ~ to!string(c)
	    ~ ", t0=" ~ to!string(t0) ~ ", t1=" ~ to!string(t1) ~ ")";
    }
    
    void evaluate_position_and_length(in double t, out Vector3 loc, out double L) const
    {
	// Both the position of the point and the length of the full arc are evaluated
	// using mostly the same process of transforming to the plane local to the arc.
	Vector3 ca, cb, tangent1, tangent2, n, cb_local;
	double ca_mag, cb_mag, theta;

	L = 0.0;
	ca = a - c; ca_mag = abs(ca);
	cb = b - c; cb_mag = abs(cb);
	if ( fabs(ca_mag - cb_mag) > 1.0e-5 ) {
	    throw new Error(text("Arc.evaluate(): radii do not match ca=",ca," cb=",cb));
	}
	// First vector in plane.
	tangent1 = Vector3(ca); tangent1.normalize(); 
	// Compute unit normal to plane of all three points.
	n = cross(ca, cb);
	if ( abs(n) > 0.0 ) {
	    n.normalize();
	} else {
	    throw new Error(text("Arc.evaluate(): cannot find plane of three points."));
	}
	// Third (orthogonal) vector is in the original plane.
	tangent2 = cross(n, tangent1); 
	// Now transform to local coordinates so that we can do 
	// the calculation of the point along the arc in 
	// the local xy-plane, with ca along the x-axis.
	cb_local = cb;
	Vector3 zero = Vector3(0.0,0.0,0.0);
	cb_local.transform_to_local_frame(tangent1, tangent2, n, zero);
	if ( fabs(cb_local.z) > 1.0e-6 ) {
	    throw new Error(text("Arc.evaluate(): problem with transformation cb_local=", cb_local));
	}
	// Angle of the final point on the arc is in the range -pi < th <= +pi.
	theta = atan2(cb_local.y, cb_local.x);
	// The length of the circular arc.
	L = theta * cb_mag;
	// Move the second point around the arc in the local xy-plane.
	theta *= t;
	loc.refx = cos(theta) * cb_mag;
	loc.refy = sin(theta) * cb_mag;
	loc.refz = 0.0;
	// Transform back to global xyz coordinates
	// and remember to add the centre coordinates.
	loc.transform_to_global_frame(tangent1, tangent2, n, c);
    } // end evaluate_position_and_length()
} // end class Arc


unittest {
    auto a = Vector3([2.0, 2.0, 0.0]);
    auto b = Vector3([1.0, 2.0, 1.0]);
    auto c = Vector3([1.0, 2.0, 0.0]);
    auto abc = new Arc(a, b, c);
    auto d = abc(0.5);
    assert(approxEqualVectors(d, Vector3(1.7071068, 2.0, 0.7071068)), "Arc");
}
