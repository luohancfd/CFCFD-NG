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
    Vector3 opUnary(string op)() if (op == "+") {
	Vector3 result = new Vector3();
	result._p[] = this._p[];
	return result;
    }

    Vector3 opUnary(string op)() if (op == "-") {
	Vector3 result = new Vector3();
	result._p[] = - this._p[];
	return result;
    }

    Vector3 opBinary(string op)(in Vector3 rhs) if (op == "+") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] + rhs._p[];
	return result;
    }

    Vector3 opBinary(string op)(in Vector3 rhs) if (op == "-") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] - rhs._p[];
	return result;
    }

    Vector3 opBinary(string op)(in double rhs) if (op == "*") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] * rhs;
	return result;
    }

    Vector3 opBinaryRight(string op)(in double lhs) if (op == "*") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] * lhs;
	return result;
    }

    Vector3 opBinary(string op)(in double rhs) if (op == "/") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] / rhs;
	return result;
    }

    // Combined assignment operators
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
    
double dot(in Vector3 v1, in Vector3 v2) {
    double result = 0.0;
    // Maybe we should be careful with underflow and overflow...
    foreach(i; 0 .. 3) result += v1._p[i] * v2._p[i];
    return result;
}

double abs(in Vector3 v) {
    return sqrt(v.dot(v));
}

Vector3 unit(in Vector3 v) {
    return new Vector3(v).normalize();
}

Vector3 cross(in Vector3 v1, in Vector3 v2) {
    Vector3 v3 = new Vector3();
    v3._p[0] = v1._p[1] * v2._p[2] - v2._p[1] * v1._p[2];
    v3._p[1] = v2._p[0] * v1._p[2] - v1._p[0] * v2._p[2];
    v3._p[2] = v1._p[0] * v2._p[1] - v2._p[0] * v1._p[1];
    return v3;
}

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
}
