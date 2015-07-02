/** gpath.d
 * Geometry-building elements for our 3D world -- one-parameter elements.
 * Note that these are geometric paths, as distinct from the file-system paths.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 *          2015-04-21 added Arc, Bezier, Polyline
 *          2015-07-02 simplify classes
 */

module gpath;

import std.conv;
import std.math;
import geom;
import bbla;
//import numeric : findRoot;

class Path {
public:
    abstract Path dup() const;
    abstract Vector3 opCall(double t) const;
    Vector3 dpdt(double t) const
    {
	// Obtain the derivative approximately, via a finite-difference.
	double dt = 0.001;
	Vector3 p0 = this.opCall(t);
	Vector3 derivative;
	if ( t+dt > 1.0 ) {
	    // t is close to the t=1.0 boundary, use a one-sided difference.
	    Vector3 pminus1 = this.opCall(t-dt);
	    derivative = (p0 - pminus1) / dt;
	} else if ( t-dt < 0.0 ) {
	    // s is close to the s=0 boundary, use a one-sided difference.
	    Vector3 pplus1 = this.opCall(t+dt);
	    derivative = (pplus1 - p0) / dt;
	} else {
	    // Not near a boundary, use central-difference.
	    Vector3 pminus1 = this.opCall(t-dt);
	    Vector3 pplus1 = this.opCall(t+dt);
	    derivative = (pplus1 - pminus1) / (2.0 * dt);
	}
	return derivative;
    }
    Vector3 d2pdt2(double t) const
    {
	// Obtain the derivative approximately, via a finite-difference.
	double dt = 0.001;
	Vector3 p0 = this.opCall(t);
	Vector3 derivative;
	if ( t+dt > 1.0 ) {
	    // t is close to the t=1.0 boundary, use a one-sided difference.
	    Vector3 pminus1 = this.opCall(t-dt);
	    Vector3 pminus2 = this.opCall(t-2*dt);
	    derivative = (p0 - 2*pminus1 + pminus2) / pow(dt,2);
	} else if ( t-dt < 0.0 ) {
	    // s is close to the s=0 boundary, use a one-sided difference.
	    Vector3 pplus1 = this.opCall(t+dt);
	    Vector3 pplus2 = this.opCall(t+2*dt);
	    derivative = (pplus2 - 2*pplus1 + p0) / pow(dt,2);
	} else {
	    // Not near a boundary, use central-difference.
	    Vector3 pminus1 = this.opCall(t-dt);
	    Vector3 pplus1 = this.opCall(t+dt);
	    derivative = (pplus1 - 2*p0 + pminus1) / pow(dt,2);
	}
	return derivative;
    }
    double partial_length(double ta, double tb) const
    {
	if( tb < ta ) {
	    double tmp = ta; ta = tb; tb = tmp;
	}
	double L = 0.0;
	int n = 100;
	double dt = (tb - ta) / n;
	Vector3 p0 = this.opCall(ta);
	Vector3 p1;
	foreach (i; 1 .. n+1) {
	    p1 = this.opCall(ta + dt * i);
	    L += abs(p1 - p0);
	    p0 = p1;
	}
	return L;
    }
    double length() const
    {
	return partial_length(0.0, 1.0);
    }
    Vector3 point_from_length(double length, out double t) const
    {
	double L = 0.0;
	int n = 1000;
	double dt = 1.0 / n;
	Vector3 p0 = this.opCall(0.0);
	Vector3 p1;
	foreach (i; 1 .. n+1) {
	    p1 = this.opCall(dt * i);
	    L += abs(p1 - p0);
	    p0 = p1;
	    if( L > length ) {
		t = dt * i;
		return p1;
	    }
	}
	t = dt * n;
	return p1;
    }
    abstract override string toString() const;
} // end class Path


class Line : Path {
public:
    Vector3 p0; // end-point at t=0
    Vector3 p1; // end-point at t=1
    this(in Vector3 p0, in Vector3 p1)
    {
	this.p0 = p0; this.p1 = p1;
    }
    this(ref const(Line) other)
    {
	p0 = other.p0; p1 = other.p1;
    }
    override Line dup() const
    {
	return new Line(p0, p1);
    }
    override Vector3 opCall(double t) const 
    {
	return (1.0-t)*p0 + t*p1;
    }
    override Vector3 dpdt(double t) const 
    {
	return p1 - p0;
    }
    override Vector3 d2pdt2(double t) const 
    {
	return Vector3(0.0,0.0,0.0);
    }
    override string toString() const
    {
	return "Line(p0=" ~ to!string(p0) ~ ", p1=" ~ to!string(p1) ~ ")";
    }
    override double partial_length(double ta, double tb) const
    {
	return fabs(tb - ta) * abs(p1 - p0);
    }
    override Vector3 point_from_length(double length, out double t) const
    {
	t = length/(abs(p1 - p0));
	return this.opCall(t);
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
    this(in Vector3 a, in Vector3 b, in Vector3 c)
    {
	this.a = a; this.b = b; this.c = c;
    }
    this(ref const(Arc) other)
    {
	a = other.a; b = other.b; c = other.c;
    }
    override Arc dup() const
    {
	return new Arc(a, b, c);
    }
    override Vector3 opCall(double t) const 
    {
	double L;
	Vector3 p;
	evaluate_position_and_length(t, p, L);
	return p;
    }
    override string toString() const
    {
	return "Arc(a=" ~ to!string(a) ~ ", b=" ~ to!string(b) ~ ", c=" ~ to!string(c) ~ ")";
    }
    override double length() const
    {
	// Returns the geometric length.
	double L;
	Vector3 p;
	evaluate_position_and_length(1.0, p, L);
	return L;
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

class Arc3 : Arc {
    // Arc constructed from start-point a, end-point b and another intermediate point m.
    // Internally it is stored as an Arc object.
    Vector3 m;
    this(in Vector3 a, in Vector3 m, in Vector3 b)
    {
	Vector3 n = cross(m - a, m - b); // normal to plane of arc
	if (abs(n) <= 1.0e-11) {
	    throw new Error(text("Arc3: Points appear colinear.",
				 " a=", to!string(a),
				 " m=", to!string(m),
				 " b=", to!string(b)));
	}
	// The centre of the circle lies along the bisector of am and 
	// the bisector of mb.
	Vector3 mid_am = 0.5 * (a + m);
	Vector3 bisect_am = cross(n, a - m);
	Vector3 mid_mb = 0.5 * (b + m);
	Vector3 bisect_mb = cross(n, b - m);
	// Solve least-squares problem to get s_am, s_mb.
	auto amatrix = new Matrix([[bisect_am.x, bisect_mb.x],
				   [bisect_am.y, bisect_mb.y],
				   [bisect_am.z, bisect_mb.z]]);
	Vector3 diff_mid = mid_mb - mid_am;
	auto rhs = new Matrix([diff_mid.x, diff_mid.y, diff_mid.z], "column");
	auto s_values = lsqsolve(amatrix, rhs);
	double s_am = s_values[0,0];
	double s_mb = s_values[1,0];
        Vector3 c = mid_am + s_am * bisect_am;
	Vector3 c_check = mid_mb + s_mb * bisect_mb;
	if (abs(c_check - this.c) > 1.0e-9) {
	    throw new Error(text("Arc3: Points inconsistent centre estimates.",
				 " c=", to!string(this.c),
				 " c_check=", to!string(c_check)));
	}
	super(a, b, c);
	this.m = m;
    }
    this(ref const(Arc3) other)
    {
	this(other.a, other.m, other.b);
    }
    override Arc3 dup() const
    {
	return new Arc3(a, m, b);
    }
} // end class Arc3

unittest {
    auto a = Vector3([2.0, 2.0, 0.0]);
    auto b = Vector3([1.0, 2.0, 1.0]);
    auto c = Vector3([1.0, 2.0, 0.0]);
    auto abc = new Arc(a, b, c);
    auto d = abc(0.5);
    assert(approxEqualVectors(d, Vector3(1.7071068, 2.0, 0.7071068)), "Arc");
    auto adb = new Arc3(a, d, b);
    assert(approxEqualVectors(d, adb(0.5)), "Arc3");
}


class Bezier : Path {
public:
    Vector3[] B; // collection of control points
    this(in Vector3[] B)
    {
	if (B.length == 0) {
	    throw new Error(text("Bezier() No control points present."));
	}
	this.B = B.dup();
    }
    this(ref const(Bezier) other)
    {
	this(other.B);
    }
    override Bezier dup() const
    {
	return new Bezier(B);
    }
    override Vector3 opCall(double t) const 
    {
	// Evaluate B(t)
	if ( B.length == 1 ) return B[0];
	size_t n_order = B.length - 1;
	// Apply de Casteljau's algorithm. 
	Vector3[] Q = B.dup(); // work array will be overwritten
	foreach (k; 0 .. n_order) {
	    foreach (i; 0 .. n_order-k) {
		Q[i] = (1.0 - t) * Q[i] + t * Q[i+1];
	    }
	}
	return Q[0];
    } // end opCall()
    override string toString() const
    {
	return "Bezier(B=" ~ to!string(B) ~ ")";
    }
} // end class Bezier


unittest {
    auto a = Vector3([2.0, 2.0, 0.0]);
    auto b = Vector3([1.0, 2.0, 1.0]);
    auto c = Vector3([1.0, 2.0, 0.0]);
    auto abc = new Arc(a, b, c);
    auto d = abc(0.5);
    auto adb = new Bezier([a, d, b]);
    auto e = adb(0.5);
    assert(approxEqualVectors(e, Vector3(1.60355, 2, 0.603553)), "Bezier");
}


class Polyline : Path {
public:
    Path[] segments; // collection of Path segments
    double[] t_values; // collection of segment break-points (in parameter t)
    this(in Path[] segments)
    {
	if (segments.length == 0) {
	    throw new Error(text("Polyline() No segments present."));
	}
	foreach (myseg; segments) this.segments ~= myseg.dup(); 
	t_values.length = this.segments.length;
	reset_breakpoints();
    }
    this(ref const(Polyline) other)
    {
	this(other.segments);
    }
    override Polyline dup() const
    {
	return new Polyline(segments);
    }
    override Vector3 opCall(double t) const 
    {
	// Evaluate B(t) without considering arc_length parameterization flag
	// or subrange.
	auto n = segments.length;
	if ( n == 1 ) return segments[0](t);
	size_t i;
	for ( i = 0; i < n; ++i ) {
	    if ( t <= t_values[i] ) break;
	}
	if ( i >= n ) i = n - 1;  // last segment
	// At this point, t_values[i-1] < t <= t_values[i] (we hope)
	// Have assumed that the t breakpoints are well behaved.
	double t_local;
	if ( i == 0 ) {
	    t_local = t / t_values[i];
	} else {
	    t_local = (t - t_values[i-1]) / (t_values[i] - t_values[i-1]);
	}
	return segments[i](t_local);
    } // end opCall()
    override string toString() const
    {
	return "Polyline(segments=" ~ to!string(segments) ~ ")";
    }

private:
    void reset_breakpoints()
    {
	// Set up the parameter breakpoints based on cumulative length.
	t_values[0] = segments[0].length();
	foreach (i; 1 .. segments.length) t_values[i] = t_values[i-1] + segments[i].length(); 
	double L_total = t_values[$-1];
	foreach (i; 0 .. segments.length) t_values[i] /= L_total; 
    } // end reset_breakpoints()
} // end class Polyline


unittest {
    auto a = Vector3([2.0, 2.0, 0.0]);
    auto b = Vector3([1.0, 2.0, 1.0]);
    auto c = Vector3([1.0, 2.0, 0.0]);
    auto abc = new Arc(a, b, c);
    auto polyline = new Polyline([abc, new Line(b, c)]);
    auto f = polyline(0.5);
    assert(approxEqualVectors(f, Vector3(1.28154, 2, 0.95955)), "Polyline");
}


class ArcLengthParameterizedPath : Path {
public:
    Path underlying_path;
    double arc_length;
    this(const Path other)
    {
	underlying_path = other.dup();
	set_arc_length_vector(100);
    }
    this(ref const(ArcLengthParameterizedPath) other)
    {
	this(other.underlying_path);
    }
    override ArcLengthParameterizedPath dup() const
    {
	return new ArcLengthParameterizedPath(this.underlying_path);
    }
    override Vector3 opCall(double t) const 
    {
	double tdsh = t_from_arc_length(t);
	return underlying_path(tdsh);
    }
    override string toString() const
    {
	return "ArcLengthParameterizedPath(underlying_path=" ~ to!string(underlying_path) ~ ")";
    }

protected:
    double[] arc_length_vector;
    void set_arc_length_vector(int N)
    {
	// Compute the arc_lengths for a number of sample points 
	// so that these can later be used to do a reverse interpolation
	// on the evaluation parameter.
	arc_length = underlying_path.length();
	arc_length_vector.length = 0;
	if ( N == 0 ) return;
	double dt = 1.0 / N;
	double L = 0.0;
	arc_length_vector ~= L;
	Vector3 p0 = underlying_path(0.0);
	Vector3 p1;
	foreach (i; 1 .. N+1) {
	    p1 = underlying_path(dt * i);
	    L += abs(p1 - p0);
	    arc_length_vector ~= L;
	    p0 = p1;
	}
    } // end set_arc_length_vector()
    double t_from_arc_length(double t) const
    {
	// The incoming parameter value, t, is proportional to arc_length.
	// Do a reverse look-up from the arc_length to the original t parameter
	// of the underlying Path.
	double L_target = t * arc_length_vector[$-1];
	// Starting from the right-hand end,
	// let's try to find a point to the left of L_target.
	// If the value is out of range, this should just result in
	// us extrapolating one of the end segments -- that's OK.
	int i = to!int(arc_length_vector.length) - 1;
	double dt = 1.0 / arc_length_vector.length;
	while ( L_target < arc_length_vector[i] && i > 0 ) i--;
	double frac = (L_target - arc_length_vector[i]) /
	    (arc_length_vector[i+1] - arc_length_vector[i]);
	return (1.0 - frac) * dt*i + frac * dt*(i+1);
    } // end t_from_arc_length()
} // end class ArcLengthParameterizedPath


class SubRangedPath : Path {
public:
    Path underlying_path;
    double t0;
    double t1;
    this(const Path other, double newt0, double newt1)
    {
	underlying_path = other.dup();
	t0 = newt0;
	t1 = newt1;
    }
    this(ref const(SubRangedPath) other)
    {
	this(other.underlying_path, other.t0, other.t1);
    }
    override SubRangedPath dup() const
    {
	return new SubRangedPath(this.underlying_path, this.t0, this.t1);
    }
    override Vector3 opCall(double t) const 
    {
	double tdsh = t0 + (t1 - t0) * t;
	return underlying_path(tdsh);
    }
    override string toString() const
    {
	return "SubRangedPath(underlying_path=" ~ to!string(underlying_path) 
	    ~ ", t0=" ~ to!string(t0) ~ " t1=" ~ to!string(t1) ~ ")";
    }
} // end class SubRangedPath


class ReversedPath : SubRangedPath {
    // Just a particular case of SubRangedPath
    this(const Path other)
    {
	super(other, 1.0, 0.0);
    }
} // end class ReversedPath
