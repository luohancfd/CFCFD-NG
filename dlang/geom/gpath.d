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
import bbla;
//import numeric : findRoot;

interface Parameterizable {
    // any public functions that use or change 't'
    // intended to have raw_ equivalents that can be overwritten if desired
    Vector3 opCall(double t) const;
    Vector3 dpdt(double t) const;
    Vector3 d2pdt2(double t) const;
    //double locate(ref const Vector3 p, double tolerance, int max_iterations, out int result_flag) const;
    double length() const;
    double partial_length(double ta, double tb) const;
    Vector3 point_from_length(double length, out double t) const;
    //Parameterizable reverse();
}

//interface Path : Parameterizable {
//    Path dup() const;
//    //Path translate(ref const Vector3 v);
//    //Path translate(double vx, double vy, double vz);
//    //Path mirror_image(ref const Vector3 point, ref const Vector3 normal);
//    //Path rotate_about_zaxis(double dtheta);
//}

class Path : Parameterizable{
public:
    double t0; // to subrange t, when evaluating a point on the Path
    double t1;
    bool arc_length_param_flag;
    double arc_length;
    this(double t0=0.0, double t1=1.0, bool arc_length_param=false){
	this.t0 = t0;
	this.t1 = t1;
	this.arc_length_param_flag = arc_length_param;
	if (arc_length_param){
	    arc_length = raw_length();
	}
    }
    this(const(Path) other)
    {
	this(other.t0, other.t1, other.arc_length_param_flag);
    }
    abstract Path dup() const;
    
    final Vector3 opCall(double t) const
    {
	t = get_raw_t(t);
	return raw_eval(t);
    }
    final Vector3 dpdt(double t) const
    {
	t = get_raw_t(t);
	return raw_dpdt(t)*d_raw_t_dt(t);
    }
    final Vector3 d2pdt2(double t) const
    {
	t = get_raw_t(t);
	return raw_d2pdt2(t)*pow(d_raw_t_dt(t), 2) + raw_dpdt(t)*d2_raw_t_dt2(t);
    }
    
    final double partial_length(double ta, double tb) const
    {
	double raw_ta = get_raw_t(ta);
	double raw_tb = get_raw_t(tb);
	return raw_partial_length(ta, tb);
    }
    final double length() const
    {
	return partial_length(0, 1);
    }
    final Vector3 point_from_length(double length, out double t) const
    {
	double raw_length = length + raw_partial_length(0, t0);
	double raw_t;
	Vector3 point = raw_point_from_length(raw_length, raw_t);
	t = get_t_from_raw(raw_t);
	return point;
    }
    //final Path reverse() 
    //{
    //    double t0_old = t0;
    //    double t1_old = t1;
    //    t0 = 1.0 - t1_old;
    //    t1 = 1.0 - t0_old;
    //    raw_reverse();
    //    return this;
    //}
    
    abstract Vector3 raw_eval(double t) const;
    Vector3 raw_dpdt(double t) const
    {
	// Obtain the derivative approximately, via a finite-difference.
	double dt = 0.001;
	Vector3 p0 = raw_eval(t);
	Vector3 derivative;
	if ( t+dt > 1.0 ) {
	    // t is close to the t=1.0 boundary, use a one-sided difference.
	    Vector3 pminus1 = raw_eval(t-dt);
	    derivative = (p0 - pminus1) / dt;
	} else if ( t-dt < 0.0 ) {
	    // s is close to the s=0 boundary, use a one-sided difference.
	    Vector3 pplus1 = raw_eval(t+dt);
	    derivative = (pplus1 - p0) / dt;
	} else {
	    // Not near a boundary, use central-difference.
	    Vector3 pminus1 = raw_eval(t-dt);
	    Vector3 pplus1 = raw_eval(t+dt);
	    derivative = (pplus1 - pminus1) / (2.0 * dt);
	}
	return derivative;
    }
    Vector3 raw_d2pdt2(double t) const
    {
	// Obtain the derivative approximately, via a finite-difference.
	double dt = 0.001;
	Vector3 p0 = raw_eval(t);
	Vector3 derivative;
	if ( t+dt > 1.0 ) {
	    // t is close to the t=1.0 boundary, use a one-sided difference.
	    Vector3 pminus1 = raw_eval(t-dt);
	    Vector3 pminus2 = raw_eval(t-2*dt);
	    derivative = (p0 - 2*pminus1 + pminus2) / pow(dt,2);
	} else if ( t-dt < 0.0 ) {
	    // s is close to the s=0 boundary, use a one-sided difference.
	    Vector3 pplus1 = raw_eval(t+dt);
	    Vector3 pplus2 = raw_eval(t+2*dt);
	    derivative = (pplus2 - 2*pplus1 + p0) / pow(dt,2);
	} else {
	    // Not near a boundary, use central-difference.
	    Vector3 pminus1 = raw_eval(t-dt);
	    Vector3 pplus1 = raw_eval(t+dt);
	    derivative = (pplus1 - 2*p0 + pminus1) / pow(dt,2);
	}
	return derivative;
    }
    double raw_partial_length(double ta, double tb) const
    {
	if( tb < ta ) {
	    double tmp = ta; ta = tb; tb = tmp;
	}
	double L = 0.0;
	int n = 100;
	double dt = (tb - ta) / n;
	Vector3 p0 = this.raw_eval(ta);
	Vector3 p1;
	foreach (i; 1 .. n+1) {
	    p1 = this.raw_eval(ta + dt * i);
	    L += abs(p1 - p0);
	    p0 = p1;
	}
	return L;
    }
    double raw_length() const
    {
	return raw_partial_length(0, 1);
    }
    Vector3 raw_point_from_length(double length, out double t) const
    {
	double L = 0.0;
	int n = 1000;
	double dt = 1.0 / n;
	Vector3 p0 = this.raw_eval(0.0);
	Vector3 p1;
	foreach (i; 1 .. n+1) {
	    p1 = this.raw_eval(dt * i);
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
    //abstract Path raw_reverse();
    
    // A global piece of data so that we can
    // create a minimum function for the optimiser
    // used to locate a point.

    // Note this is a minimisation problem because we are
    // trying to approach zero.  It is hard to write this
    // as a zero-finding problem because we do not cross the
    // axis, that is, the scalar length between the guessed
    // value and the desired value is *always* a positive number.

    //static Vector3 p_sought;
    //static Path *path = 0;
    //double error_in_point(double t)
    //{
    //    return vabs(p_sought - path->eval(t));
    //}

    // Note this function can only find *one* location on
    // the path where the path crosses p.  This may not
    // be suitable if your path crosses back on itself
    // AND you happen to be looking for the particular
    // points where the path crosses.
    // TODO need optimisation library
    //double locate(ref const Vector3 p, double tolerance, int max_iterations, out int result_flag) const;
    //{
    //    // Do this to mimic "closure"
    //    p_sought = p;
    //    path = this;

    //    double t = golden_section_search(error_in_point, 0.0, 1.0, result_flag, tolerance, max_iterations);

    //    return t;
    //}
    //abstract Path translate(ref const Vector3 v);
    //abstract Path translate(double vx, double vy, double vz);
    //abstract Path mirror_image(ref const Vector3 point, ref const Vector3 normal);
    //abstract Path rotate_about_zaxis(double dtheta);
    abstract override string toString() const;
protected:
    double[] arc_length_vector;
    void set_arc_length_vector(int N)
    {
	// Compute the arc_lengths for a number of sample points 
	// so that these can later be used to do a reverse interpolation
	// on the evaluation parameter.
	arc_length = raw_length();
	arc_length_vector.length = 0;
	if ( N == 0 ) return;
	double dt = 1.0 / N;
	double L = 0.0;
	arc_length_vector ~= L;
	Vector3 p0 = raw_eval(0.0);
	Vector3 p1;
	foreach (i; 1 .. N+1) {
	    p1 = raw_eval(dt * i);
	    L += abs(p1 - p0);
	    arc_length_vector ~= L;
	    p0 = p1;
	}
    } // end set_arc_length_vector()
    double t_from_arc_length_f(double t) const
    {
	// The incoming parameter value, t, is proportional to arc_length.
	// Do a reverse look-up from the arc_length to the original t parameter
	// of the Bezier curve.
	double L_target = t * arc_length_vector[$-1];
	// Starting from the right-hand end,
	// let's try to find a point to the left of L_target.
	// If the value is out of range, this should just result in
	// us extrapolating one of the end segments -- that's OK.
	int i = to!int(arc_length_vector.length) - 1;
	double dt = 1.0 / arc_length_vector.length;
	while ( L_target < arc_length_vector[i] && i > 0 ) i--;
	double frac = (L_target - arc_length_vector[i]) / (arc_length_vector[i+1] - arc_length_vector[i]);
	return (1.0 - frac) * dt*i + frac * dt*(i+1);
    } // end t_from_arc_length_f()
private:
    double get_raw_t(double t) const
    {
	t = t_from_subrange(t);
	if ( arc_length_param_flag ) {
	    t = t_from_arc_length_f(t);
	}
	return t;
    }
    double get_t_from_raw(double t) const
    {
	if ( arc_length_param_flag ) {
	    t = raw_partial_length(0,t)/arc_length;
	}
	t = t_to_subrange(t);
	return t;
    }
    double d_raw_t_dt(double t) const
    {
	if ( arc_length_param_flag ) {
	    return d_t_from_subrange_dt(t)*d_t_from_arc_f_dt(t);
	}
	return d_t_from_subrange_dt(t);
    }
    double d2_raw_t_dt2(double t) const
    {
	// d2_t_from_subrange_dt2 == 0
	if ( arc_length_param_flag ) {
	    return pow(d_t_from_subrange_dt(t),2)*d2_t_from_arc_f_dt2(t);
	}
	return 0.0;
    }
    double t_from_subrange(double t) const
    {
	return t0 + (t1 - t0)*t;
    }
    double d_t_from_subrange_dt(double t) const
    {
	return t1 - t0;
    }
    double t_to_subrange(double t) const
    {
	return (t - t0)/(t1 - t0);
    }
    double d_arc_f_dt(double t) const
    {
	return abs(raw_dpdt(t))/arc_length;
    }
    double d2_arc_f_dt2(double t) const
    {
	//chain rule on d_arc_f_dt
	return dot(unit(raw_dpdt(t)),raw_d2pdt2(t))/arc_length;
    }
    double d_t_from_arc_f_dt(double t) const
    {
	// input is parameter t not arclength fraction
	// derivative of inverse fn
	return 1.0/d_arc_f_dt(t);
    }
    double d2_t_from_arc_f_dt2(double t) const
    {
	// input is parameter t not arclength fraction
	// derivative of inverse fn
	return -d2_arc_f_dt2(t)/pow(d_arc_f_dt(t),3);
    }
} // end class Path


class Line : Path {
public:
    Vector3 p0; // end-point at t=0
    Vector3 p1; // end-point at t=1
    this(in Vector3 p0, in Vector3 p1, double t0=0.0, double t1=1.0, bool arc_length_param=false)
    {
	super(t0, t1, arc_length_param);
	this.p0 = p0; this.p1 = p1;
    }
    this(ref const(Line) other)
    {
	super(other);
	p0 = other.p0; p1 = other.p1;
    }
    override Line dup() const
    {
	return new Line(p0, p1, t0, t1);
    }
    override Vector3 raw_eval(double t) const 
    {
	return (1.0-t)*p0 + t*p1;
    }
    override Vector3 raw_dpdt(double t) const 
    {
	return p1 - p0;
    }
    override Vector3 raw_d2pdt2(double t) const 
    {
	return Vector3(0.0,0.0,0.0);
    }
    override string toString() const
    {
	return "Line(p0=" ~ to!string(p0) ~ ", p1=" ~ to!string(p1) ~
	    ", t0=" ~ to!string(t0) ~ ", t1=" ~ to!string(t1) ~ ")";
    }
    override double raw_partial_length(double ta, double tb) const
    {
	return fabs(tb - ta) * abs(p1 - p0);
    }
    override Vector3 raw_point_from_length(double length, out double t) const
    {
	t = length/(abs(p1 - p0));
	return raw_eval(t);
    }
    //override Line raw_reverse()
    //{
    //    Vector3 tmp;
    //    tmp = p0; p0 = p1; p1 = tmp;
    //    return this;
    //}
protected:
    override void set_arc_length_vector(int N){}
    override double t_from_arc_length_f(double t) const {return t;}
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
    override Vector3 raw_eval(double t) const 
    {
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
    override double raw_length() const
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
    this(in Vector3 a, in Vector3 m, in Vector3 b, double t0=0.0, double t1=1.0)
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
	super(a, b, c, t0, t1);
	this.m = m;
    }
    this(ref const(Arc3) other)
    {
	this(other.a, other.m, other.b, other.t0, other.t1);
    }
    override Arc3 dup() const
    {
	return new Arc3(a, m, b, t0, t1);
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
	this(other.B, other.t0, other.t1, other.arc_length_param_flag);
    }
    override Bezier dup() const
    {
	return new Bezier(B, t0, t1, arc_length_param_flag);
    }
    override Vector3 raw_eval(double t) const 
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
    override string toString() const
    {
	return "Bezier(B=" ~ to!string(B) 
	    ~ ", t0=" ~ to!string(t0) 
	    ~ ", t1=" ~ to!string(t1) 
	    ~ ", arc_length_p=" ~ to!string(arc_length_param_flag)
	    ~ ")";
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
	this(other.segments, other.t0, other.t1, other.arc_length_param_flag);
    }
    override Polyline dup() const
    {
	return new Polyline(segments, t0, t1, arc_length_param_flag);
    }

    override Vector3 raw_eval(double t) const 
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
    override string toString() const
    {
	return "Polyline(segments=" ~ to!string(segments) 
	    ~ ", t0=" ~ to!string(t0) 
	    ~ ", t1=" ~ to!string(t1) 
	    ~ ", arc_length_p=" ~ to!string(arc_length_param_flag)
	    ~ ")";
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
