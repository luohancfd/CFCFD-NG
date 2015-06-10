/** flowsolution.d
 * Eilmer4 compressible-flow simulation code, postprocessing functions.
 *
 * The role of the post-processing functions is just to pick up data
 * from a previously-run simulation and either write plotting files
 * or extract interesting pieces of data.  To do this, we really don't
 * need or want all of the code machinery in the classes built for the 
 * simulation code so we rebuild a bit of the core data handling here.
 * This allows us to be a bit flexible in what variables are held.
 * We might want to add flow variables, such as Mach number or Pitot
 * pressure, whish are not normally held in the main simulation data
 * structures.  It also frees us of the internal numbering of cells in
 * the simulation code that must allow for ghost-cells.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-06-09
 */

module flowsolution;

import std.stdio;
import std.conv;
import std.format;
import std.string;
import std.algorithm;
import std.array;
import std.math;
import gzip;
import fileutil;
import geom;
import sgrid;
import gas;
import globalconfig;


class FlowSolution {
    // The collection of flow blocks and grid blocks that define the flow
    // and the domain at one particular instant in time.
public:
    double sim_time;
    size_t nBlocks;
    SBlockFlow[] flowBlocks;
    StructuredGrid[] gridBlocks;
    
    this(string jobName, string dir, int tindx, size_t nBlocks)
    {
	foreach (ib; 0 .. nBlocks) {
	    string fileName;
	    if (GlobalConfig.moving_grid) {
		fileName = make_file_name!"grid"(jobName, to!int(ib), tindx);
	    } else {
		fileName = make_file_name!"grid"(jobName, to!int(ib), 0);
	    }
	    fileName = dir ~ "/" ~ fileName;
	    gridBlocks ~= new StructuredGrid(fileName, GridFileFormat.gziptext);
	    fileName = make_file_name!"flow"(jobName, to!int(ib), tindx);
	    fileName = dir ~ "/" ~ fileName;
	    flowBlocks ~= new SBlockFlow(fileName);
	} // end foreach ib
	this.nBlocks = nBlocks;
	sim_time = flowBlocks[0].sim_time;
    } // end constructor

    size_t[] find_nearest_cell_centre(double x, double y, double z=0.0)
    {
	size_t[] nearest = [0, 0, 0, 0]; // blk_id, i, j, k
	double dx = x - flowBlocks[0]["pos.x", 0, 0, 0];
	double dy = y - flowBlocks[0]["pos.y", 0, 0, 0];
	double dz = z - flowBlocks[0]["pos.z", 0, 0, 0];
	double minDist = sqrt(dx*dx + dy*dy + dz*dz);
	foreach (ib; 0 .. nBlocks) {
	    auto flow = flowBlocks[ib];
	    foreach (k; 0 .. flow.nkc) {
		foreach (j; 0 .. flow.njc) {
		    foreach (i; 0 .. flow.nic) {
			dx = x - flowBlocks[ib]["pos.x", i, j, k];
			dy = y - flowBlocks[ib]["pos.y", i, j, k];
			dz = z - flowBlocks[ib]["pos.z", i, j, k];
			double dist = sqrt(dx*dx + dy*dy + dz*dz);
			if (dist < minDist) {
			    minDist = dist;
			    nearest = [ib, i, j, k];
			}
		    } // foreach i
		} // foreach j
	    } // foreach k
	} // foreach ib
	return nearest;
    } // end find_nearest_cell_centre()
} // end class FlowSolution

class SBlockFlow {
    // Much like the Python library for postprocessing in Eilmer3,
    // we are going to handle the data as a big chunk of numbers,
    // with the label for each variable coming the top of the file.
public:
    size_t nic;
    size_t njc;
    size_t nkc;
    string[] variableNames;
    size_t[string] variableIndex;
    double sim_time;

    this(string filename)
    {
	// Read in the flow data for a single structured block.
	string[] tokens;
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, " %g", &sim_time);
	line = byLine.front; byLine.popFront();
	variableNames = line.strip().split();
	foreach (ref var; variableNames) { var = removechars(var, "\""); }
	foreach (i; 0 .. variableNames.length) { variableIndex[variableNames[i]] = i; }
	line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &nic, &njc, &nkc);
	_data.length = nic;
	// Resize the storage for our block of data.
	foreach (i; 0 .. nic) {
	    _data[i].length = njc;
	    foreach (j; 0 .. njc) {
		_data[i][j].length = nkc;
		foreach (k; 0 .. nkc) {
		    _data[i][j][k].length = variableNames.length;
		} // foreach k
	    } // foreach j
	} // foreach i
	// Scan the remainder of the file, extracting our data.
	foreach (k; 0 .. nkc) {
	    foreach (j; 0 .. njc) {
		foreach (i; 0 ..nic) {
		    line = byLine.front; byLine.popFront();
		    tokens = line.strip().split();
		    assert(tokens.length == variableNames.length, "wrong number of items");
		    foreach (ivar; 0 .. variableNames.length) {
			_data[i][j][k][ivar] = to!double(tokens[ivar]);
		    }
		} // foreach i
	    } // foreach j
	} // foreach k
    } // end constructor from file

    ref double opIndex(string varName, size_t i, size_t j, size_t k=0)
    {
	return _data[i][j][k][variableIndex[varName]];
    }

    string variable_names_as_string()
    {
	auto writer = appender!string();
	formattedWrite(writer, "#");
	foreach(name; variableNames) {
	    formattedWrite(writer, " \"%s\"", name);
	}
	return writer.data;
    }

    string values_as_string(size_t i, size_t j, size_t k)
    {
	auto writer = appender!string();
	formattedWrite(writer, "%e", _data[i][j][k][0]);
	foreach (ivar; 1 .. variableNames.length) {
	    formattedWrite(writer, " %e", _data[i][j][k][ivar]);
	}
	return writer.data;
    }

private:
    double[][][][] _data;
} // end class SBlockFlowData
