/** volume.d
 * Geometry-building elements for our 3D world -- three-parameter volumes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-04-07 first code
 */

module volume;

import std.math;
import std.stdio;
import std.conv;
import geom;
import gpath;
import surface;

// Nomenclature for the parametric distances, bounding surfaces, paths and corners.
//
// t=1 at top surface
//         north
//    p011-------p111 s=1
//      |         |
// west |   Top   | east
//      |         |
//    p001-------p101 s=0
//         south
//     r=0       r=1
//
//
// t=0 at Bottom surface
//         north
//    p010-------p110 s=1
//      |         |
// west |  Bottom | east
//      |         |
//    p000-------p100 s=0
//         south
//     r=0       r=1
//
// Faces:
// North = face[0]; East = face[1]; South = face[2]; West = face[3];
// Top = face[4]; Bottom = face[5]
//
// Corners:
// Bottom surface: p000 == p[0]; p100 == p[1]; p110 == p[2]; p010 == p[3]
// Top surface   : p001 == p[4]; p101 == p[5]; p111 == p[6]; p011 == p[7]
//
// Edges:
// edge[0] p[0] --> p[1] around Bottom surface
//     [1] p[1] --> p[2]
//     [2] p[3] --> p[2]
//     [3] p[0] --> p[3]
//
//     [4] p[4] --> p[5] around Top surface
//     [5] p[5] --> p[6]
//     [6] p[7] --> p[6]
//     [7] p[4] --> p[7]
//
//     [8] p[0] --> p[4] connecting Bottom to Top
//     [9] p[1] --> p[5]
//    [10] p[2] --> p[6]
//    [11] p[3] --> p[7]
//
// We'll try to use this notation consistently in the classes below.

class ParametricVolume {
public:
    double r0; // to subrange r, when evaluating a point on the surface
    double r1;
    double s0;
    double s1;
    double t0;
    double t1;

    Vector3 opCall(double r, double s, double t) const
    {
	return Vector3(0.0, 0.0, 0.0);
    }
    ParametricVolume dup() const
    {
	return new ParametricVolume();
    }
    override string toString() const
    {
	return "ParametricVolume()";
    }
} // end class ParametricVolume

class TFIVolume : ParametricVolume {
public:
    ParametricSurface[6] faces;
    Vector3[8] p;
    // Path[12] edges;

    // Generic, 6-face constructor
    this(in ParametricSurface[] faceArray,
	 double r0=0.0, double r1=1.0,
	 double s0=0.0, double s1=1.0,
	 double t0=0.0, double t1=1.0)
    {
	foreach(i; 0 .. 6) faces[i] = faceArray[i].dup();
	this.r0 = r0; this.r1 = r1;
	this.s0 = s0; this.s1 = s1;
	this.t0 = t0; this.t1 = t1;
	// Check that the corners of the faces coincide.
	Vector3 p000 = faces[Face.bottom](0.0,0.0);
    }

    // Wire-Frame constructor TODO

    // Simple-Box constructor

    // Surface-extrusion constructor TODO

}

class MeshVolume : ParametricVolume {
}
