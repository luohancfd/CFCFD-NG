/**
 * luageom.d  LuaD connection for the geometric elements.
 *
 * Author: Peter J and Rowan G.
 * Version: 2014-Feb-21 First pieces of the experiment 
 *          extracted from geom.d and gpath.d.
 */

module luageom;

import std.conv;
import std.stdio;
import std.math;
import geom;
import gpath;
import luad.all;

struct Vector3Params {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
}

// A little database of Vector3 objects
static Vector3[] points;

void Vector3FillTable(ref LuaTable table, int i) {
    if (i >= 0 && i < points.length) {
	table["x"] = points[i].x;
	table["y"] = points[i].y;
	table["z"] = points[i].z;
    } else {
	table["x"] = nil; table["y"] = nil; table["z"] = nil;
    }
    return;
}
double[] Vector3AsArray(int i)
{
    return points[i]._p;
}
double Vector3GetX(int i)
{
    if (i >= 0 && i < points.length)
	return points[i].x;
    else
	return 0.0;
}
double Vector3GetY(int i)
{
    if (i >= 0 && i < points.length)
	return points[i].y;
    else
	return 0.0;
}
double Vector3GetZ(int i)
{
    if (i >= 0 && i < points.length)
	return points[i].z;
    else
	return 0.0;
}

int makeVector3(LuaTable t)
// Returns index of newly added point.
{
    auto data = t.toStruct!Vector3Params();
    points ~= Vector3(data.x, data.y, data.z);
    return cast(int) points.length - 1;
}
int makeVector3FromArray(double[] xyz)
{
    points ~= Vector3(xyz);
    return cast(int) points.length - 1;
}

Vector3 checkVector3Table(ref LuaTable t)
{
    // If we are going use it as a table
    // fill in the missing spots
    if ( t["x"].isNil ) t["x"] = 0.0;
    if ( t["y"].isNil ) t["y"] = 0.0;
    if ( t["z"].isNil ) t["z"] = 0.0;
    auto data = t.toStruct!Vector3Params();
    return Vector3(data.x, data.y, data.z);
}

LuaTable makeVector3FromTable(LuaTable t)
{
    checkVector3Table(t);
    return t;
}

void registerVector3(LuaState lua)
// Register the Vector3 service functions with the Lua interpreter.
{
    lua["Vector"] = &makeVector3;
    lua["Vector3Value"] = &Vector3FillTable;
    lua["Vector3Array"] = &Vector3AsArray;
    lua["VectorA"] = &makeVector3FromArray;
    lua["getX"] = &Vector3GetX;
    lua["getY"] = &Vector3GetY;
    lua["getZ"] = &Vector3GetZ;
    lua["Vector3"] = &makeVector3FromTable;
}

//---------------------------------------------------------------------

struct LineParams {
    int[] plist;
}

// A little database of Path objects
static Path[] paths;

int makeLine(int p0, int p1)
// Makes a new Path object and returns its index.
{
    paths ~= new Line(points[p0], points[p1]);
    return cast(int) paths.length - 1;
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

// RJG additions to experiment
Line checkLineTable(ref LuaTable t)
{
    if ( t["t0"].isNil ) t["t0"] = 0.0;
    if ( t["t1"].isNil ) t["t1"] = 1.0;
    auto p0Tab = t.get!LuaTable("p0");
    auto p0 = checkVector3Table(p0Tab);
    auto p1Tab = t.get!LuaTable("p1");
    auto p1 = checkVector3Table(p1Tab);
    return new Line(p0, p1, t.get!double("t0"), t.get!double("t1"));
}

static LuaState my_lua_ref;
LuaTable makeLine2(LuaTable p0Tab, LuaTable p1Tab)
{
    auto lTab = my_lua_ref.newTable();
    lTab["p0"] = p0Tab;
    lTab["p1"] = p1Tab;
    checkLineTable(lTab);
    return lTab;
}

LuaTable evalLine2(ref LuaTable lTab, double s)
{
    auto line = checkLineTable(lTab);
    auto p = line(s);
    auto t = my_lua_ref.newTable();
    t["x"] = p.x; t["y"] = p.y; t["z"] = p.z;
    return t;
}

void registerPath(LuaState lua)
// Register the Line service functions with the Lua interpreter.
{
    my_lua_ref = lua;
    lua["Line"] = &makeLine;
    lua["evalLine"] = &evalLine;
    lua["Line2"] = &makeLine2;
    lua["evalLine2"] = &evalLine2;
}
