/**
 * sgrid.d -- structured-grid functions
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-02-22 First code ported from e3_grid.py
 */

module sgrid;

import std.string;
import std.array;
import std.conv;
import std.stdio;
import std.format;
import std.math;
import gzip;

import geom;
import gpath;
import surface;
// import volume;
import univariatefunctions;

//-----------------------------------------------------------------

enum GridFileFormat {
    text,
    gziptext,
    vtk,
    vtkxml
}

class StructuredGrid {
public:
    int niv, njv, nkv;
    Vector3[][][] grid;
    string label;

    this(const ParametricSurface surf, int niv, int njv,
	 const(UnivariateFunction)[] clusterf,
	 string label="empty-label")
    {
	this.niv = niv; this.njv = njv; this.nkv = 1;
	// Any unspecified clustering functions defaults to the linear identity.
	while (clusterf.length < 4) clusterf ~= new LinearFunction(0.0, 1.0);
	resize_array();
	make_grid_from_surface(surf, clusterf);
    }

    this(string fileName, GridFileFormat fmt)
    {
	final switch (fmt) {
	case GridFileFormat.text:
	    read_grid_from_text_file(fileName);
	    break;
	case GridFileFormat.gziptext:
	    throw new Error("Reading from GZip text format not implemented.");
	case GridFileFormat.vtk:
	    read_grid_from_text_file(fileName, true);
	    break;
	case GridFileFormat.vtkxml:
	    throw new Error("Reading from VTK XML format not implemented.");
	}
    }

    void resize_array()
    {
	grid.length = niv;
	foreach (i; 0 .. niv) {
	    grid[i].length = njv;
	    foreach (j; 0 .. njv) {
		grid[i][j].length = nkv;
	    }
	}
    } // end resize_array

    void make_grid_from_surface(const ParametricSurface surf,
				const(UnivariateFunction)[] clusterf)
    {
	// First, set up clustered parameter values along each edge.
        double[] rNorth = clusterf[0].distribute_parameter_values(niv);
        double[] sEast = clusterf[1].distribute_parameter_values(njv);
        double[] rSouth = clusterf[2].distribute_parameter_values(niv);
        double[] sWest = clusterf[3].distribute_parameter_values(njv);
	// Now, work through the mesh, one point at a time,
        // blending the stretched parameter values
        // and creating the actual vertex coordinates in Cartesian space.
        int k = 0;
        foreach (j; 0 .. njv) {
            double s = to!double(j) / (njv - 1);
	    foreach (i; 0 .. niv) {
                double r = to!double(i) / (niv - 1);
                double sdash = (1.0-r) * sWest[j] + r * sEast[j]; 
                double rdash = (1.0-s) * rSouth[i] + s * rNorth[i];
                Vector3 p = surf(rdash, sdash);
                grid[i][j][k].refx = p.x;
                grid[i][j][k].refy = p.y;
                grid[i][j][k].refz = p.z;
	    }
	}
    } // end make_grid_from_surface()

    void write_to_text_file(string fileName, bool vtkHeader=true)
    {
	auto f = File(fileName, "w");
	if (vtkHeader) {
	    f.writeln("# vtk DataFile Version 2.0");
	    f.writeln(label);
	    f.writeln("ASCII");
	    f.writeln("");
	    f.writeln("DATASET STRUCTURED_GRID");
	    f.writefln("DIMENSIONS %d %d %d", niv, njv, nkv);
	    f.writefln("POINTS %d float", (niv * njv * nkv));
	} else {
	    f.writefln("%d %d %d  # ni nj nk", niv, njv, nkv);
	}
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    f.writefln("%20.16e %20.16e %20.16e", 
			       grid[i][j][k].x, grid[i][j][k].y, grid[i][j][k].z);
		}
	    }
	}
    } // end write_to_text_file()

    void write_grid_to_gzip_file(string fileName)
    // This function essentially defines the Eilmer4 native format.
    {
	auto f = new GzipOut(fileName);
	auto writer = appender!string();
	formattedWrite(writer, "%d %d %d  # ni nj nk\n", niv, njv, nkv);
	f.compress(writer.data);
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    writer = appender!string();
		    formattedWrite(writer, "%20.16e %20.16e %20.16e\n", 
				   grid[i][j][k].x, grid[i][j][k].y, grid[i][j][k].z);
		    f.compress(writer.data);
		}
	    }
	}
	f.finish();
    } // end write_grid_to_gzip_file()

    void read_grid_from_text_file(string fileName, bool vtkHeader=true)
    {
	string[] tokens;
	auto f = File(fileName, "r");
	if (vtkHeader) {
	    read_VTK_header_line("vtk", f);
	    label = f.readln().strip();
	    read_VTK_header_line("ASCII", f);
	    read_VTK_header_line("STRUCTURED_GRID", f);
	    tokens = read_VTK_header_line("DIMENSIONS", f);
	} else {
	    tokens = f.readln().strip().split();
	}
	niv = to!int(tokens[0]);
	njv = to!int(tokens[1]);
	nkv = to!int(tokens[2]);
	resize_array();
	foreach (k; 0 .. nkv) {
	    foreach (j; 0 .. njv) {
		foreach (i; 0 .. niv) {
		    tokens = f.readln().strip.split();
		    try {
			grid[i][j][k].refx = to!double(tokens[0]);
			grid[i][j][k].refy = to!double(tokens[1]);
			grid[i][j][k].refz = to!double(tokens[2]);
		    } catch (Exception e) {
			throw new Error(text("Failed to read grid file at "
					     "i=", i, " j=", j, " k=", k,
					     "tokens=", tokens, "exception=", e));
		    }
		} // foreach i
	    } // foreach j
	} // foreach k
    } // end read_grid_from_text_file()

} // end class SGrid

//-----------------------------------------------------------------
// Helper functions

// Set very small quantities to zero, exactly.
//
// This is intended primarily to avoid the bad behaviour of VTK
// when it is reading Float32 values that are *too* small.
// We have also come across unreasonably-small float values 
// in the context of reading GridPro files.
double uflowz(double q, double tiny=1.0e-30)
{
    return (fabs(q) > tiny) ? q: 0.0;
}

// Locate the line containing target and return the tokens on that line.
// Legacy-format VTK lines use spaces as delimiters.
string[] read_VTK_header_line(string target, File f)
{
    bool found = false; 
    string[] tokens;
    while (!found) {
        auto line = f.readln();
        if (line.length == 0) break; // presume end of file
	line = strip(line);
        if (indexOf(line, target, CaseSensitive.no) > -1) {
            tokens = split(line);
            found = true; break;
	}
    }
    if (!found) { 
        throw new Error(text("Did not find ", target, " while reading VTK grid file"));
    }
    return tokens;
} // end locate_VTK_header_line()

//-----------------------------------------------------------------

unittest {
    auto p00 = Vector3(0.0, 0.1);
    auto p10 = Vector3(1.0, 0.1);
    auto p11 = Vector3(1.0, 1.1);
    auto p01 = Vector3(0.0, 1.1);
    auto my_patch = new CoonsPatch(p00, p10, p11, p01);
    auto cf = [new LinearFunction(), new LinearFunction(), 
	       new LinearFunction(), new LinearFunction()];
    auto my_grid = new StructuredGrid(my_patch, 11, 21, cf);
    assert(approxEqualVectors(my_grid.grid[5][5][0], Vector3(0.5, 0.35, 0.0)),
			      "StructuredGrid sample point");
}
