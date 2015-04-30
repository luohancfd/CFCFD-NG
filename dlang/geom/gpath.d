/** gpath.d
 * Geometry-building elements for our 3D world -- one-parameter elements.
 * Note that these are geometric paths, as distinct from the file-system paths.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 *          2015-04-21 added Arc, Bezier, Polyline
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
    double length() const
    {
	// Returns the geometric length of the Path.
	double L = 0.0;
	int n = 20;
	double dt = 1.0 / n;
	Vector3 p0 = this.opCall(0.0);
	Vector3 p1;
	foreach (i; 1 .. n) {
	    p1 = this.opCall(dt * i);
	    L += abs(p1 - p0);
	    p0 = p1;
	}
	return L;
    } // end length()
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
    override double length() const
    {
	// Returns the geometric length.
	return fabs(t1 - t0) * abs(p1 - p0);
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
	evaluate_position_and_length(tdsh, p, L);
	return p;
    }
    override string toString() const
    {
	return "Arc(a=" ~ to!string(a) ~ ", b=" ~ to!string(b) ~ ", c=" ~ to!string(c)
	    ~ ", t0=" ~ to!string(t0) ~ ", t1=" ~ to!string(t1) ~ ")";
    }
    override double length() const
    {
	// Returns the geometric length.
	double L;
	Vector3 p;
	evaluate_position_and_length(1.0, p, L);
	return fabs(t1 - t0) * L;
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


class Bezier : Path {
public:
    Vector3[] B; // collection of control points
    bool arc_length_param_flag;
    this(in Vector3[] B, double t0=0.0, double t1=1.0, bool arc_length_p=false)
    {
	this.B = B.dup();
	this.t0 = t0; this.t1 = t1;
	arc_length_param_flag = arc_length_p;
	if ( arc_length_param_flag ) set_arc_length_vector(100);
    }
    this(ref const(Bezier) other)
    {
	B = other.B.dup();
	t0 = other.t0; t1 = other.t1;
	arc_length_param_flag = other.arc_length_param_flag;
	if ( arc_length_param_flag ) arc_length = other.arc_length.dup();
    }
    override Bezier dup() const
    {
	return new Bezier(B, t0, t1);
    }

    Vector3 raw_eval(double t) const 
    {
	// Evaluate B(t) without considering arc_length parameterization flag
	// or subrange.
	if ( B.length == 0 ) {
	    throw new Error(text("Bezier.raw_eval() No control points present."));
	}
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
    } // end raw_eval()

    override Vector3 opCall(double t) const 
    {
	// Evaluate B(t) considering arc_length parameterization flag 
	// and possible subrange of t.
	if ( B.length == 0 ) {
	    throw new Error(text("Bezier.opCall() No control points present."));
	}
	if ( B.length == 1 ) return B[0];
    t = t0 + (t1 - t0)*t; // subrange t
	if ( arc_length_param_flag ) {
	    t = t_from_arc_length(t);
	}
	return raw_eval(t);
    } // end OpCall()

    override string toString() const
    {
	return "Bezier(B=" ~ to!string(B) 
	    ~ ", t0=" ~ to!string(t0) 
	    ~ ", t1=" ~ to!string(t1) 
	    ~ ", arc_length_p=" ~ to!string(arc_length_param_flag)
	    ~ ")";
    }

private:
    double[] arc_length;
    void set_arc_length_vector(int N)
    {
	// Compute the arc_lengths for a number of sample points 
	// so that these can later be used to do a reverse interpolation
	// on the evaluation parameter.
	arc_length.length = 0;
	if ( N == 0 ) return;
	double dt = 1.0 / N;
	double L = 0.0;
	arc_length ~= L;
	Vector3 p0 = raw_eval(0.0);
	Vector3 p1;
	foreach (i; 1 .. N+1) {
	    p1 = raw_eval(dt * i);
	    L += abs(p1 - p0);
	    arc_length ~= L;
	    p0 = p1;
	}
    } // end set_arc_length_vector()

    double t_from_arc_length(double t) const
    {
	// The incoming parameter value, t, is proportional to arc_length.
	// Do a reverse look-up from the arc_length to the original t parameter
	// of the Bezier curve.
	double L_target = t * arc_length[$-1];
	// Starting from the right-hand end,
	// let's try to find a point to the left of L_target.
	// If the value is out of range, this should just result in
	// us extrapolating one of the end segments -- that's OK.
	int i = to!int(arc_length.length) - 1;
	double dt = 1.0 / arc_length.length;
	while ( L_target < arc_length[i] && i > 0 ) i--;
	double frac = (L_target - arc_length[i]) / (arc_length[i+1] - arc_length[i]);
	return (1.0 - frac) * dt*i + frac * dt*(i+1);
    } // end t_from_arc_length()

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
    bool arc_length_param_flag;
    this(in Path[] segments, double t0=0.0, double t1=1.0, bool arc_length_p=false)
    {
	foreach (myseg; segments) this.segments ~= myseg.dup(); 
	t_values.length = this.segments.length;
	this.t0 = t0; this.t1 = t1;
	arc_length_param_flag = arc_length_p;
	if ( arc_length_param_flag ) set_arc_length_vector(100);
	reset_breakpoints();
    }
    this(ref const(Polyline) other)
    {
	foreach (myseg; other.segments) segments ~= myseg.dup(); 
	t_values = other.t_values.dup();
	t0 = other.t0; t1 = other.t1;
	arc_length_param_flag = other.arc_length_param_flag;
	if ( arc_length_param_flag ) arc_length = other.arc_length.dup();
	reset_breakpoints();
    }
    override Polyline dup() const
    {
	return new Polyline(segments, t0, t1);
    }

    Vector3 raw_eval(double t) const 
    {
	// Evaluate B(t) without considering arc_length parameterization flag
	// or subrange.
	auto n = segments.length;
	if ( n == 0 ) {
	    throw new Error(text("Polyline.raw_eval() No segments present."));
	}
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
    } // end raw_eval()

    override Vector3 opCall(double t) const 
    {
	// First, map to the subrange of the full polyline.
	t = t0 + t * (t1 - t0);
	if ( arc_length_param_flag ) t = t_from_arc_length(t);
	return raw_eval(t);
    } // end OpCall()

    override string toString() const
    {
	return "Polyline(segments=" ~ to!string(segments) 
	    ~ ", t0=" ~ to!string(t0) 
	    ~ ", t1=" ~ to!string(t1) 
	    ~ ", arc_length_p=" ~ to!string(arc_length_param_flag)
	    ~ ")";
    }

private:
    double[] arc_length;
    void reset_breakpoints()
    {
	// Set up the parameter breakpoints based on cumulative length.
	t_values[0] = segments[0].length();
	foreach (i; 1 .. segments.length) t_values[i] = t_values[i-1] + segments[i].length(); 
	double L_total = t_values[$-1];
	foreach (i; 0 .. segments.length) t_values[i] /= L_total; 
    } // end reset_breakpoints()

    void set_arc_length_vector(int N)
    {
	// Compute the arc_lengths for a number of sample points 
	// so that these can later be used to do a reverse interpolation
	// on the evaluation parameter.
	arc_length.length = 0;
	if ( N == 0 ) return;
	double dt = 1.0 / N;
	double L = 0.0;
	arc_length ~= L;
	Vector3 p0 = raw_eval(0.0);
	Vector3 p1;
	foreach (i; 1 .. N+1) {
	    p1 = raw_eval(dt * i);
	    L += abs(p1 - p0);
	    arc_length ~= L;
	    p0 = p1;
	}
    } // end set_arc_length_vector()

    double t_from_arc_length(double t) const
    {
	// The incoming parameter value, t, is proportional to arc_length.
	// Do a reverse look-up from the arc_length to the original t parameter
	// of the Bezier curve.
	double L_target = t * arc_length[$-1];
	// Starting from the right-hand end,
	// let's try to find a point to the left of L_target.
	// If the value is out of range, this should just result in
	// us extrapolating one of the end segments -- that's OK.
	int i = to!int(arc_length.length) - 1;
	double dt = 1.0 / arc_length.length;
	while ( L_target < arc_length[i] && i > 0 ) i--;
	double frac = (L_target - arc_length[i]) / (arc_length[i+1] - arc_length[i]);
	return (1.0 - frac) * dt*i + frac * dt*(i+1);
    } // end t_from_arc_length()

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
