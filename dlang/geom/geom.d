/**
 * geom.d  Geometric (vector) primitives for our 3D world.
 *
 * Author: Peter J.
 * Version: 2014-06-16 first cut.
 */
module geom;

import std.conv;
import std.stdio;
import std.math;

class Vector3 {
    double[] _p;

    this() { 
	_p = new double[3];
	_p[0] = _p[1] = _p[2] = 0.0;
    }

    this(in double[] p) {
	_p = new double[3];
	switch ( p.length ) {
	case 0: _p[0] = _p[1] = _p[2] = 0.0; break;
	case 1: _p[0] = p[0]; _p[1] = _p[2] = 0.0; break;
	case 2: _p[0] = p[0]; _p[1] = p[1]; _p[2] = 0.0; break;
	default: _p[0] = p[0]; _p[1] = p[1]; _p[2] = p[2]; break;
	}
    }

    this(in double p0, in double p1=0.0, in double p2=0.0) {
	_p = new double[3];
	_p[0] = p0;
	_p[1] = p1;
	_p[2] = p2;
    }

    this(in Vector3 other) {
	_p = other._p.dup;
    }

    // Note that the following three properties hand out references
    // to the elements, so we can't control changes to the elements.
    @property ref double x() { return _p[0]; }
    @property ref double y() { return _p[1]; }
    @property ref double z() { return _p[2]; }

    @property Vector3 dup() {
	return new Vector3(this);
    }

    override string toString() {
	return "Vector3(" ~ to!string(_p) ~ ")";
    }

    // Some operators, at least those that make sense.
    // Note that they make new values and do not alter the original object.
    // Note 2, the const does not refer to the returned (new) object.
    const Vector3 opUnary(string op)() if (op == "+") {
	Vector3 result = new Vector3();
	result._p[] = this._p[];
	return result;
    }

    const Vector3 opUnary(string op)() if (op == "-") {
	Vector3 result = new Vector3();
	result._p[] = - this._p[];
	return result;
    }

    const Vector3 opBinary(string op)(in Vector3 rhs) if (op == "+") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] + rhs._p[];
	return result;
    }

    const Vector3 opBinary(string op)(in Vector3 rhs) if (op == "-") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] - rhs._p[];
	return result;
    }

    const Vector3 opBinary(string op)(in double rhs) if (op == "*") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] * rhs;
	return result;
    }

    const Vector3 opBinaryRight(string op)(in double lhs) if (op == "*") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] * lhs;
	return result;
    }

    const Vector3 opBinary(string op)(in double rhs) if (op == "/") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] / rhs;
	return result;
    }

    // Combined assignment operators do change the original object.
    ref Vector3 opOpAssign(string op)(in Vector3 rhs) if (op == "+") {
	this._p[] += rhs._p[];
	return this;
    }

    ref Vector3 opOpAssign(string op)(in Vector3 rhs) if (op == "-") {
	this._p[] -= rhs._p[];
	return this;
    }

    ref Vector3 opOpAssign(string op)(in double rhs) if (op == "*") {
	this._p[] *= rhs;
	return this;
    }

    ref Vector3 opOpAssign(string op)(in double rhs) if (op == "/") {
	this._p[] /= rhs;
	return this;
    }

    // Other vector-specific operations.

    /**
     * Scales the vector to unit magnitude.
     */
    ref Vector3 normalize() {
        double magnitude = abs(this);
	if ( magnitude > 0.0 ) {
	    this /= magnitude;
	} else {
	    // Clean up, in case dot() underflows.
	    this._p[0] = this._p[1] = this._p[2] = 0.0;
	}
	return this;
    }
} // end class Vector3


/**
 * Returns the scalar dot product of two vectors.
 */   
double dot(in Vector3 v1, in Vector3 v2) {
    double result = 0.0;
    // Maybe we should be careful with underflow and overflow...
    foreach(i; 0 .. 3) result += v1._p[i] * v2._p[i];
    return result;
}

/**
 * Returns magnitude of the vector.
 */
double abs(in Vector3 v) {
    return sqrt(v.dot(v));
}

/**
 * Returns a unit vector in the same direction as v.
 */
Vector3 unit(in Vector3 v) {
    return new Vector3(v).normalize();
}

/**
 * Vector cross product.
 */
Vector3 cross(in Vector3 v1, in Vector3 v2) {
    Vector3 v3 = new Vector3();
    v3._p[0] = v1._p[1] * v2._p[2] - v2._p[1] * v1._p[2];
    v3._p[1] = v2._p[0] * v1._p[2] - v1._p[0] * v2._p[2];
    v3._p[2] = v1._p[0] * v2._p[1] - v2._p[0] * v1._p[1];
    return v3;
}

// Transform functions used to reorient vector values in the CFD codes.

/**
 * Rotate v from the global xyz coordinate system into the local frame
 * defined by the orthogonal unit vectors n,t1,t2.
 *
 * We assume, without checking, that these vectors do nicely define 
 * such a local system.
 */
void to_local_frame(ref Vector3 v, in Vector3 n, in Vector3 t1, in Vector3 t2)
{
    double v_x = dot(v, n); // normal component
    double v_y = dot(v, t1); // tangential component 1
    double v_z = dot(v, t2); // tangential component 2
    v._p[0] = v_x;
    v._p[1] = v_y;
    v._p[2] = v_z;
}

/**
 * Rotate v back into the global coordinate system.
 */
void to_xyz_frame(ref Vector3 v, in Vector3 n, in Vector3 t1, in Vector3 t2)
{
    double v_x = v._p[0]*n._p[0] + v._p[1]*t1._p[0] + v._p[2]*t2._p[0]; // global-x
    double v_y = v._p[0]*n._p[1] + v._p[1]*t1._p[1] + v._p[2]*t2._p[1]; // global-y
    double v_z = v._p[0]*n._p[2] + v._p[1]*t1._p[2] + v._p[2]*t2._p[2]; // global-z
    v._p[0] = v_x;
    v._p[1] = v_y;
    v._p[2] = v_z;
}

/**
 * Returns true is all of the components are approximately equal.
 */
bool approxEqualVectors(in Vector3 v1, in Vector3 v2) {
    return (approxEqual(v1._p[0], v2._p[0]) && 
	    approxEqual(v1._p[1], v2._p[1]) &&
	    approxEqual(v1._p[2], v2._p[2]));
}

unittest {
    // Check that we have separate data with the correct values.
    Vector3 a = new Vector3([1.0, 2.2, 3.0]);
    Vector3 b = new Vector3(1.0);
    assert(a.x == 1.0, "a.x fail");
    assert(a.y == 2.2, "a.y fail");
    assert(a.z == 3.0, "a.z fail");
    assert(a.x == b.x, "a.x == b.x fail");
    assert(b.y == 0.0, "b.y fail");
    assert(b.z == 0.0, "b.z fail");
    // Check operators
    b = -a;
    assert(b.x == -a.x && b.y == -a.y && b.z == -a.z, "unary negation fail");

    b = new Vector3(1.0);
    Vector3 c = a + b;
    assert(c.y == a.y+b.y, "vector addition fail");
    c = a - b;
    assert(c.y == a.y-b.y, "vector subtraction fail");
    Vector3 d = a.dup;
    a.y = 99.0;
    assert(a.y == 99.0 && d.y == 2.2, "dup followed by vector addition fail");

    Vector3 e = a * 2.0;
    Vector3 f = 3.0 * d;
    assert(e.z == 6.0 && f.z == 9.0, "scalar multiply fail");
    Vector3 g = d / 3.0;
    assert(g.z == 1.0, "scalar division fail");

    g += f;
    assert(g.z == 10.0, "plus-assign fail");
    g /= 2.0;
    assert(g.z == 5.0, "divide-assign fail");

    Vector3 u = unit(g);
    assert(approxEqual(abs(u), 1.0), "unit, dot, abs fail");

    Vector3 x = new Vector3(1.0, 0.0, 0.0);
    Vector3 y = new Vector3(0.0, 1.0, 0.0);
    Vector3 z = cross(x,y);
    Vector3 zref = new Vector3(0.0,0.0,1.0);
    assert(approxEqualVectors(z, zref), "cross product fail");

    Vector3 n = unit(new Vector3(1.0,1.0,0.0));
    Vector3 t1 = unit(new Vector3(-1.0,1.0,0.0));
    Vector3 t2 = cross(n, t1);
    Vector3 h = new Vector3(1.0,0.0,1.0);
    Vector3 h_ref = new Vector3(h);
    to_local_frame(h, n, t1, t2);
    assert(approxEqualVectors(h, new Vector3(sqrt(1.0/2.0), -sqrt(1.0/2.0), 1.0)),
	   "to_local_frame fail");
    to_xyz_frame(h, n, t1, t2);
    assert(approxEqualVectors(h, h_ref), "to_xyz_frame fail");
}


//------------------------------------------------------------------------
// Geometry functions that make use of Vector3 objects...

/**
 * Project the point q onto the plane, along the vector qr.
 */
int project_onto_plane(ref Vector3 q, in Vector3 qr,
		       in Vector3 a, in Vector3 b, in Vector3 c)
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
} // end project_onto_plane()

unittest {
    Vector3 a = new Vector3(1.0, 0.0, 0.0); // plane through a,b,c
    Vector3 b = new Vector3(1.0, 1.0, 0.0);
    Vector3 c = new Vector3(0.5, 0.0, 0.0);
    Vector3 qr = new Vector3(3.0, 3.0, -3.0); // direction
    Vector3 q = new Vector3(0.0, 0.0, 1.0); // start point
    int flag =  project_onto_plane(q, qr, a, b, c);
    assert(approxEqualVectors(q, new Vector3(1.0,1.0,0.0)), "project_onto_plane fail");
}
