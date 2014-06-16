/**
 * geom.d  Geometric primitives for our 3D world.
 *
 * Author: Peter J.
 * Version: 2014-06-16 first cut.
 */

import std.conv;
import std.stdio;

class Vector3 {
    double[] _p;

    this() { 
	_p = new double[3];
	_p[0] = _p[1] = _p[2] = 0.0;
    }

    this(double[] p) {
	_p = new double[3];
	switch ( p.length ) {
	case 0: _p[0] = _p[1] = _p[2] = 0.0; break;
	case 1: _p[0] = p[0]; _p[1] = _p[2] = 0.0; break;
	case 2: _p[0] = p[0]; _p[1] = p[1]; _p[2] = 0.0; break;
	default: _p[0] = p[0]; _p[1] = p[1]; _p[2] = p[2]; break;
	}
    }

    this(double p0, double p1=0.0, double p2=0.0) {
	_p = new double[3];
	_p[0] = p0;
	_p[1] = p1;
	_p[2] = p2;
    }

    this(Vector3 other) {
	_p = other._p.dup;
    }

    @property ref double x() { return _p[0]; }
    @property ref double y() { return _p[1]; }
    @property ref double z() { return _p[2]; }

    @property Vector3 dup() {
	return new Vector3(this);
    }

    override string toString() {
	return "Vector3(" ~ to!string(_p) ~ ")";
    }

    Vector3 opBinary(string op)(Vector3 rhs) if (op == "+") {
	Vector3 result = new Vector3();
	result._p[] = this._p[] + rhs._p[];
	return result;
    }
} // end class Vector3

unittest {
    // Check that we have separate data with the correct values.
    Vector3 a = new Vector3([1.0, 2.2, 3.0]);
    Vector3 b = new Vector3(1.0);
    assert( a.x == 1.0, "a.x failed" );
    assert( a.y == 2.2, "a.y failed" );
    assert( a.z == 3.0, "a.z failed" );
    assert( a.x == b.x, "a.x == b.x failed" );
    assert( b.y == 0.0, "b.y failed" );
    assert( b.z == 0.0, "b.z failed" );

    Vector3 c = a + b;
    assert( c.y == a.y+b.y );
    Vector3 d = a.dup;
    a.y = 99.0;
    assert( a.y == 99.0 && d.y == 2.2 );   
}
