/**
 * fvvertex.d
 * Finite-volume cell-vertex class for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvvertex;

import std.string;
import std.conv;
import fvcore;
import geom;
import gasmodel;

class FVVertex {
public:
    size_t id;  // allows us to work out where, in the block, the vertex is
    // Geometry
    Vector3[] pos;  // x,y,z-Coordinates for time-levels, m
    Vector3[] vel;  // vertex velocity for time-levels, m/s
    double areaxy;    // x,y-plane area of secondary cells (for spatial derivatives)
    double volume;  // volume of 3D secondary cells (for spatial derivatives)
    // Derivatives of primary-cell variables.
    double[][] grad_vel; // velocity derivatives stored as a second-order tensor
                         // [[du/dx du/dy du/dz]
                         //  [dv/dx dv/dy dv/dz]
                         //  [dw/dx dw/dy dw/dz]]
    Vector3 grad_T;      // Temperature derivatives (static temperature only)
    Vector3 grad_tke;    // turbulence kinetic energy
    Vector3 grad_omega;  // pseudo vorticity for k-omega turbulence
    Vector3[] grad_f;    // mass fraction derivatives
    Vector3 grad_pe;     // electron pressure derivatives

    this(in GasModel gm, size_t id_init=0)
    {
	id = id_init;
	pos.length = n_time_levels;
	vel.length = n_time_levels;
	grad_vel.length = 3;
	foreach(ref e; grad_vel) e.length = 3;
	grad_f.length = gm.n_species;
    }

    this(in FVVertex other)
    {
	id = other.id;
	pos = other.pos.dup;
	vel = other.vel.dup;
	areaxy = other.areaxy;
	volume = other.volume;
	grad_vel.length = 3;
	foreach(i; 0 .. 3) grad_vel[i] = other.grad_vel[i].dup; 
	grad_T = other.grad_T;
	grad_tke = other.grad_tke;
	grad_omega = other.grad_omega;
	foreach(i; 0 .. grad_f.length) grad_f[i] = other.grad_f[i];
	grad_pe = other.grad_pe;
    }

    void copy_values_from(in FVVertex other)
    {
	if ( !(this is other) ) {
	    id = other.id;
	    pos[] = other.pos[];
	    vel[] = other.vel[];
	    areaxy = other.areaxy;
	    volume = other.volume;
	    foreach(i; 0 .. vel.length) grad_vel[i][] = other.grad_vel[i][]; 
	    grad_T = other.grad_T;
	    grad_tke = other.grad_tke;
	    grad_omega = other.grad_omega;
	    grad_f[] = other.grad_f[];
	    grad_pe = other.grad_pe;
	}
    }

    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
	pos[to_level] = pos[from_level];
	vel[to_level] = vel[from_level];
    }

    override string toString()
    {
	char[] repr;
	repr ~= "FVVertex(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ", pos=" ~ to!string(pos);
	repr ~= ", vel=" ~ to!string(vel);
	repr ~= ", areaxy=" ~ to!string(areaxy);
	repr ~= ", volume=" ~ to!string(volume);
	repr ~= ", grad_vel=" ~ to!string(grad_vel);
	repr ~= ", grad_T=" ~ to!string(grad_T);
	repr ~= ", grad_tke=" ~ to!string(grad_tke);
	repr ~= ", grad_omega=" ~ to!string(grad_omega);
	repr ~= ", grad_f=" ~ to!string(grad_f);
	repr ~= ", grad_pe=" ~ to!string(grad_pe);
	repr ~= ")";
	return to!string(repr);
    }

/+ [TODO]
    int copy_grid_level_to_level(size_t from_level, size_t to_level);
+/
} // end class FVVertex
