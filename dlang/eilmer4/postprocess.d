/** postprocess.d
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

module postprocess;

import std.stdio;
import std.conv;
import std.format;
import std.string;
import gzip;
import fileutil;
import geom;
import sgrid;
import gas;
import globalconfig;
import readconfig;

SBlockFlow[] flowBlocks;
StructuredGrid[] gridBlocks;

void post_process(string tindxPlot, bool vtkxmlFlag)
{
    writeln("Begin post_process()...");
    read_config_file();
    auto jobName = GlobalConfig.base_file_name;
    if (tindxPlot == "all") {
	// Read the times file for all tindx values.
    }
    if (vtkxmlFlag) {
	begin_Visit_file(jobName, GlobalConfig.nBlocks);
	begin_PVD_file(jobName);
    }
    int[] tindx_list = [to!int(tindxPlot)]; // [TODO]
    foreach (tindx; tindx_list) {
	foreach (ib; 0 .. GlobalConfig.nBlocks) {
	    string fileName;
	    if (GlobalConfig.moving_grid) {
		fileName = make_file_name!"grid"(jobName, ib, tindx);
	    } else {
		fileName = make_file_name!"grid"(jobName, ib, 0);
	    }
	    writeln("grid fileName=", fileName);
	    gridBlocks ~= new StructuredGrid(fileName, GridFileFormat.gziptext);
	    fileName = make_file_name!"flow"(jobName, ib, tindx);
	    writeln("flow fileName=", fileName);
	    flowBlocks ~= new SBlockFlow(fileName);
	}
	// write_VTK_XML_files(jobName, tindx, times_dict[tindx]); // [TODO]
    } // foreach tindx
    if (vtkxmlFlag) {
	finish_PVD_file(jobName);
    }
    writeln("Done post_process().");
} // end post_process()

string unquote(string s) {
    return removechars(s, "\"");
}

string quote(string s) {
    return "\"" ~ s ~ "\"";
}

class SBlockFlow {
    // Much like the Python library for postprocessing in Eilmer3,
    // we are going to handle the data as a big chunk of numbers,
    // with the label for each variable coming the top of the file.
public:
    size_t nic;
    size_t njc;
    size_t nkc;
    string[] variableNames;
    int[string] variableIndex;
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
	foreach (ref var; variableNames) var = unquote(var);
	writeln("variableNames=", variableNames);
	line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &nic, &njc, &nkc);
	writeln("nic=", nic, " njc=", njc, " nkc=", nkc);
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
    } // end this(fileName)

    ref double opIndex(string varName, size_t i, size_t j, size_t k=1)
    {
	return _data[i][j][k][variableIndex[varName]];
    }

private:
    double[][][][] _data;
} // end class SBlockFlowData


void begin_Visit_file(string rootName, int nblock)
{
    // Will be handy to have a Visit file, also.
    // For each time index, this justs lists the names of the files for individual blocks.
    string plotPath = "plot";
    ensure_directory_is_present(plotPath);
    string fileName = plotPath ~ "/" ~ rootName ~ ".visit";
    auto visitFile = File(fileName, "w");
    visitFile.writef("!NBLOCKS %d\n", nblock);
    visitFile.close();
    return;
}

void begin_PVD_file(string rootName)
{
    // Will be handy to have a Paraview collection file, also.
    // For each time index, this justs lists the name of the top-level .pvtu file.
    string plotPath = "plot";
    ensure_directory_is_present(plotPath);
    string fileName = plotPath ~ "/" ~ rootName ~ ".pvd";
    auto pvdFile = File(fileName, "w");
    pvdFile.write("<?xml version=\"1.0\"?>\n");
    pvdFile.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    pvdFile.write("<Collection>\n");
    pvdFile.close();
    return;
}

void finish_PVD_file(string rootName)
{
    string plotPath = "plot";
    ensure_directory_is_present(plotPath);
    string fileName = plotPath ~ "/" ~ rootName ~ ".pvd";
    auto pvdFile = File(fileName, "a");
    pvdFile.write("</Collection>\n");
    pvdFile.write("</VTKFile>\n");
    pvdFile.close();
    return;
}
