/** surface.d
 * Geometry-building elements for our 3D world -- two-parameter surfaces.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 */

module surface;

import std.conv;
import geom;
import gpath;


class ParametricSurface {
public:
    double r0; // to subrange r, when evaluating a point on the surface
    double r1;
    double s0;
    double s1;

    Vector3 eval(double r, double s) const
    {
	return Vector3(0.0, 0.0, 0.0);
    }
    ParametricSurface dup() const
    {
	return new ParametricSurface();
    }
    override string toString() const
    {
	return "ParametricSurface()";
    }
} // end class ParametricSurface


class CoonsPatch : ParametricSurface {
public:
    //         north
    //     p01-------p11 s=1
    //      |         |
    // west |         | east
    //      |         |
    //     p00-------p10 s=0
    //         south
    //     r=0       r=1
    Path north, east, south, west; // bounding paths
    Vector3 p00, p10, p11, p01;    // corners

    this(in Vector3 p00, in Vector3 p10, in Vector3 p11, in Vector3 p01,
	 double r0=0.0, double r1=1.0, double s0=0.0, double s1=1.0)
    {
	this.p00 = p00; this.p10 = p10;
	this.p11 = p11; this.p01 = p01;
	this.r0 = r0; this.r1 = r1;
	this.s0 = s0; this.s1 = s1;
	north = new Line(p01, p11);
	east = new Line(p10, p11);
	south = new Line(p00, p10);
	west = new Line(p00, p01);
    }

    this(in Path north, in Path east, in Path south, in Path west,
	 double r0=0.0, double r1=1.0, double s0=0.0, double s1=1.0)
    {
	this.north = north.dup();
	this.east = east.dup();
	this.south = south.dup();
	this.west = west.dup();
	this.r0 = r0; this.r1 = r1;
	this.s0 = s0; this.s1 = s1;
	p00 = south.eval(0.0);
	p10 = south.eval(1.0);
	p01 = north.eval(0.0);
	p11 = north.eval(1.0);
	// TODO check alternate evaluation of corners for consistency.
    }

    this(ref const(CoonsPatch) other)
    {
	this.north = other.north.dup();
	this.east = other.east.dup();
	this.south = other.south.dup();
	this.west = other.west.dup();
	this.r0 = other.r0; this.r1 = other.r1;
	this.s0 = other.s0; this.s1 = other.s1;
	p00 = other.p00;
	p10 = other.p10;
	p01 = other.p01;
	p11 = other.p11;
    }

    override CoonsPatch dup() const
    {
	return new CoonsPatch(this.north, this.east, this.south, this.west,
			      r0, r1, s0, s1);
    }

    override Vector3 eval(double r, double s) const 
    {
	r = r0 + (r1-r0)*r; // subrange
	s = s0 + (s1-s0)*s;
	Vector3 south_r = south.eval(r); 
	Vector3 north_r = north.eval(r);
	Vector3 west_s = west.eval(s); 
	Vector3 east_s = east.eval(s);
	Vector3 p = (1.0-s)*south_r + s*north_r + (1.0-r)*west_s + r*east_s - 
	    ( (1.0-r)*(1.0-s)*p00 + (1.0-r)*s*p01 + r*(1.0-s)*p10 + r*s*p11 );
	return p;
    }

    override string toString() const
    {
	return "CoonsPatch(north=" ~ to!string(north) ~
	    ", east=" ~ to!string(east) ~
	    ", south=" ~ to!string(south) ~
	    ", west=" ~ to!string(west) ~
	    ", r0=" ~ to!string(r0) ~ ", r1=" ~ to!string(r1) ~
	    ", s0=" ~ to!string(s0) ~ ", s1=" ~ to!string(s1) ~
	    ")";
    }
} // end class CoonsPatch


unittest {
    auto p00 = Vector3([0.0, 0.1, 3.0]);
    auto p10 = Vector3(1.0, 0.1, 3.0);
    auto p11 = Vector3(1.0, 1.1, 3.0);
    auto p01 = Vector3(0.0, 1.1, 3.0);
    auto my_patch = new CoonsPatch(p00, p10, p11, p01);
    auto c = my_patch.eval(0.5, 0.5);
    assert(approxEqualVectors(c, Vector3(0.5, 0.6, 3.0)), "CoonsPatch.eval mid");
    c = my_patch.eval(0.1, 0.1);
    assert(approxEqualVectors(c, Vector3(0.1, 0.2, 3.0)), "CoonsPatch.eval lower-left");
}
