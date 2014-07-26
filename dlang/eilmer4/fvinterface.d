/**
 * fvinterface.d
 * Finite-volume cell-interface class, for use in the CFD codes.
 * Fluxes of conserved quantities are transported (between cells) across cell interfaces.

 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
 */

module fvinterface;

import std.conv;
import geom;
import gasmodel;
import fvcore;
import flowstate;
import conservedquantities;

class FVInterface {
public:
    size_t id;  // allows us to work out where, in the block, the interface is
    // Geometry
    Vector3 pos;           // position of the (approx) midpoint
    Vector3 vel;           // face (grid) velocity, m/s
    double Ybar;           // Y-coordinate of the mid-point
    double length;         // Interface length in the x,y-plane
    double[] area;         // Area m**2 for each time-level.
                           // Area per radian in axisymmetric geometry
    Vector3 n;             // Direction cosines for unit normal
    Vector3 t1;            // tangent vector 1 (aka p)
    Vector3 t2;            // tangent vector 2 (aka q)
    // Flow
    FlowState fs;          // Flow properties
    ConservedQuantities F; // Flux conserved quantity per unit area
    // [TODO] Point-implicit variables

    this(in GasModel gm, size_t id_init=0)
    {
	id = id_init;
	area.length = n_time_levels;
	fs = new FlowState(gm, 100.0e3, [300.0,], Vector3(0.0,0.0,0.0));
	F = new ConservedQuantities(gm);
    }

    this(in FVInterface other, in GasModel gm)
    {
	id = other.id;
	pos = other.pos;
	vel = other.vel;
	Ybar = other.Ybar;
	length = other.length;
	area = other.area.dup;
	n = other.n;
	t1 = other.t1;
	t2 = other.t2;
	fs = new FlowState(other.fs, gm);
	F = new ConservedQuantities(other.F);
    }

    void copy_values_from(in FVInterface other, uint type_of_copy)
    {
	switch ( type_of_copy ) {
	case CopyDataOption.flow:
	    fs.copy_values_from(other.fs);
	    F.copy_values_from(other.F);
	    break;
	case CopyDataOption.grid:
	    pos = other.pos;
	    vel = other.vel;
	    Ybar = other.Ybar;
	    length = other.length;
	    area[] = other.area[];
	    n = other.n;
	    t1 = other.t1;
	    t2 = other.t2;
	    break;
	case CopyDataOption.all: 
	default:
	    id = other.id;
	    pos = other.pos;
	    vel = other.vel;
	    Ybar = other.Ybar;
	    length = other.length;
	    area[] = other.area[];
	    n = other.n;
	    t1 = other.t1;
	    t2 = other.t2;
	    fs.copy_values_from(other.fs);
	    F.copy_values_from(other.F);
	} // end switch
    }

    void copy_grid_level_to_level(uint from_level, uint to_level)
    {
	area[to_level] = area[from_level];
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FVInterface(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ", pos=" ~ to!string(pos);
	repr ~= ", vel=" ~ to!string(vel);
	repr ~= ", Ybar=" ~ to!string(Ybar);
	repr ~= ", length=" ~ to!string(length);
	repr ~= ", area=" ~ to!string(area);
	repr ~= ", n=" ~ to!string(n);
	repr ~= ", t1=" ~ to!string(t1);
	repr ~= ", t2=" ~ to!string(2);
	repr ~= ", fs=" ~ to!string(fs);
	repr ~= ", F=" ~ to!string(F);
	repr ~= ")";
	return to!string(repr);
    }
} // end of class FV_Interface
