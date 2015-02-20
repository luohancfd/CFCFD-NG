/** gpath.d
 * Geometry-building elements for our 3D world -- one-parameter elements.
 * Note that these are geometric paths, as distinct from the file-system paths.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 */

module gpath;

import std.conv;
import geom;


class Path {
public:
    double t0; // to subrange t, when evaluating a point on the Path
    double t1;
    Path dup() const { return new Path(); }
    Vector3 opCall(double t) const
    {
	return Vector3(0.0, 0.0, 0.0);
    }
    override string toString() const
    {
	return "Path()";
    }
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

//-----------------------------------------------------------------------------
// Connection to Lua following Rowan's lead.

import luad.all;

struct LineParams {
    int[] plist;
}

// A little database of Path objects
static Path[] paths;

int makeLine(int p0, int p1)
// Makes a new Path object and returns its index.
{
    paths ~= new Line(points[p0], points[p1]);
    return paths.length - 1;
}

void evalLine(ref LuaTable table, int i, double s)
// Places the position of the parametric point s into the table t.
{
    if (i >= 0 && i < paths.length) {
	auto p = paths[i](s);
	table["x"] = p.x; table["y"] = p.y; table["z"] = p.z;
    } else {
	table["x"] = nil; table["y"] = nil; table["z"] = nil;
    }

    return;
}

void registerPath(LuaState lua)
// Register the Line service functions with the Lua interpreter.
{
    lua["Line"] = &makeLine;
    lua["evalLine"] = &evalLine;
}
