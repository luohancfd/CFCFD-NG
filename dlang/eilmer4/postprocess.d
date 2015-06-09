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
import std.algorithm;
import gzip;
import fileutil;
import geom;
import sgrid;
import gas;
import globalconfig;
import readconfig;

// Some module-level data that is shared across post-processing functions.
string jobName;
SBlockFlow[] flowBlocks;
StructuredGrid[] gridBlocks;
double[int] times_dict;
int[] tindx_list;
string plotPath = "plot";

void post_process(string tindxPlot, bool vtkxmlFlag, bool binary_format)
{
    writeln("Begin post_process()...");
    read_config_file();
    jobName = GlobalConfig.base_file_name;
    readTimesFile();
    writeln("tindx_list=", tindx_list);
    int[] tindx_list_to_plot;
    switch (tindxPlot) {
    case "all":
	tindx_list_to_plot = tindx_list.dup;
	break;
    case "9999":
    case "last":
	tindx_list_to_plot ~= tindx_list[$-1];
        break;
    default:
	// We assume that the command-line argument was an integer.
        tindx_list_to_plot ~= to!int(tindxPlot);
    } // end switch
    if (vtkxmlFlag) {
	begin_Visit_file(GlobalConfig.nBlocks);
	begin_PVD_file();
    }
    foreach (tindx; tindx_list_to_plot) {
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
	} // end foreach ib
	write_VTK_XML_files(tindx);
    } // foreach tindx
    if (vtkxmlFlag) {
	finish_PVD_file();
    }
    writeln("Done post_process().");
} // end post_process()

void readTimesFile()
{
    // Read the times file for all tindx values.
    auto timesFile = File(jobName ~ ".times");
    auto line = timesFile.readln().strip();
    while (line.length > 0) {
	if (line[0] != '#') {
	    // Process a non-comment line.
	    auto tokens = line.split();
	    times_dict[to!int(tokens[0])] = to!double(tokens[1]);
	}
	line = timesFile.readln().strip();
    }
    timesFile.close();
    tindx_list = times_dict.keys;
    sort(tindx_list);
} // end readTimesFile()

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
	foreach (ref var; variableNames) { var = unquote(var); }
	foreach (i; 0 .. variableNames.length) { variableIndex[variableNames[i]] = i; }
	writeln("variableNames=", variableNames);
	writeln("variableIndex=", variableIndex);
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


void begin_Visit_file(int nblock)
{
    // Will be handy to have a Visit file, also.
    // For each time index, this justs lists the names of the files for individual blocks.
    ensure_directory_is_present(plotPath);
    string fileName = plotPath ~ "/" ~ jobName ~ ".visit";
    auto visitFile = File(fileName, "w");
    visitFile.writef("!NBLOCKS %d\n", nblock);
    visitFile.close();
    return;
}

void begin_PVD_file()
{
    // Will be handy to have a Paraview collection file, also.
    // For each time index, this justs lists the name of the top-level .pvtu file.
    ensure_directory_is_present(plotPath);
    string fileName = plotPath ~ "/" ~ jobName ~ ".pvd";
    auto pvdFile = File(fileName, "w");
    pvdFile.write("<?xml version=\"1.0\"?>\n");
    pvdFile.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    pvdFile.write("<Collection>\n");
    pvdFile.close();
    return;
}

void finish_PVD_file()
{
    ensure_directory_is_present(plotPath);
    string fileName = plotPath ~ "/" ~ jobName ~ ".pvd";
    auto pvdFile = File(fileName, "a");
    pvdFile.write("</Collection>\n");
    pvdFile.write("</VTKFile>\n");
    pvdFile.close();
    return;
}

void write_VTK_XML_files(int tindx, bool binary_format=false)
{
    ensure_directory_is_present(plotPath);
    string fileName = plotPath~"/"~jobName~format(".t%04d", tindx)~".pvtu";
    writeln("fileName=", fileName);
    auto pvtuFile = File(fileName, "w");
    pvtuFile.write("<VTKFile type=\"PUnstructuredGrid\">\n");
    pvtuFile.write("<PUnstructuredGrid GhostLevel=\"0\">");
    pvtuFile.write("<PPoints>\n");
    pvtuFile.write(" <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    pvtuFile.write("</PPoints>\n");
    pvtuFile.write("<PCellData>\n");
    foreach (var; flowBlocks[0].variableNames) {
        pvtuFile.writef(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"/>\n", var);
    }
    pvtuFile.write(" <PDataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    if (canFind(flowBlocks[0].variableNames,"c.x")) {
	pvtuFile.write(" <PDataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    }
    if (canFind(flowBlocks[0].variableNames,"B.x")) {
	pvtuFile.write(" <PDataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    }
    pvtuFile.write("</PCellData>\n");
    foreach (jb; 0 .. flowBlocks.length) {
        fileName = jobName ~ format(".b%04d.t%04d.vtu", jb, tindx);
        // We write the short version of the fileName into the pvtu file.
        pvtuFile.writef("<Piece Source=\"%s\"/>\n", fileName);
        // but use the long version to actually open it.
        fileName = plotPath ~ "/" ~ fileName;
        write_VTK_XML_unstructured_file(fileName, jb, binary_format);
    }
    pvtuFile.write("</PUnstructuredGrid>\n");
    pvtuFile.write("</VTKFile>\n");
    pvtuFile.close();
    // Will be handy to have a Visit file, also.
    // This justs lists the names of the files for individual blocks.
    fileName = plotPath ~ "/" ~ jobName ~ ".visit";
    // Note that we append to the visit file for each tindx.
    auto visitFile = File(fileName, "a");
    foreach (jb; 0 .. flowBlocks.length) {
        fileName = jobName ~ format(".b%04d.t%04d.vtu", jb, tindx);
        visitFile.writef("%s\n", fileName);
    }
    visitFile.close();
    // Will be handy to have a Paraview PVD file, also.
    // This justs lists the top-level .pvtu files.
    fileName = plotPath ~ "/" ~ jobName ~ ".pvd";
    // Note that we append to the .pvd file for each tindx.
    auto pvdFile = File(fileName, "a");
    fileName = jobName ~ format(".t%04d.pvtu", tindx);
    pvdFile.writef("<DataSet timestep=\"%e\" group=\"\" part=\"0\" file=\"%s\"/>\n",
		   times_dict[tindx], fileName);
    pvdFile.close();
    return;
}

void write_VTK_XML_unstructured_file(string fileName, size_t jb, bool binary_format)
// Write the cell-centred flow data from a single block (index jb)
// as an unstructured grid of finite-volume cells.
{
    auto fp = File(fileName, "wb"); // We may be writing some binary data.
    auto flow = flowBlocks[jb];
    auto grid = gridBlocks[jb];
    if (binary_format) {
        //binary_data_string = "";
        //binary_data_offset = 0;
    }
    size_t niv = grid.niv; size_t njv = grid.njv; size_t nkv = grid.nkv;
    size_t nic = flow.nic; size_t njc = flow.njc; size_t nkc = flow.nkc;
    bool two_D = (nkv == 1);
    size_t NumberOfPoints = niv * njv * nkv;
    size_t NumberOfCells = nic * njc * nkc;
    fp.write("<VTKFile type=\"UnstructuredGrid\" byte_order=\"BigEndian\">\n");
    fp.write("<UnstructuredGrid>");
    fp.writef("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
	      NumberOfPoints, NumberOfCells);
    //
    fp.write("<Points>\n");
    fp.write(" <DataArray type=\"Float32\" NumberOfComponents=\"3\"");
    if (binary_format) {
        //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        //binary_data = r""
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    size_t vtx_number = 0;
    size_t[][][] vtx_id;
    vtx_id.length = niv;
    foreach (i; 0 .. niv) {
	vtx_id[i].length = njv;
	foreach (j; 0 .. njv) {
	    vtx_id[i][j].length = nkv;
	}
    }
    foreach (k; 0 .. nkv) {
        foreach (j; 0 .. njv) {
            foreach (i; 0 .. niv) {
                vtx_id[i][j][k] = vtx_number;
                float x = uflowz(grid.grid[i][j][k].x);
		float y = uflowz(grid.grid[i][j][k].y);
		float z = uflowz(grid.grid[i][j][k].z);
                if (binary_format) {
                    // binary_data += struct.pack('> f f f', x, y, z)
                } else {
                    fp.writef(" %e %e %e\n", x,y,z);
		}
                vtx_number += 1;
	    }
	}
    }
    fp.write(" </DataArray>\n");
    if (binary_format) {
        //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        //binary_data_string += binary_data_count
        //binary_data_string += binary_data
        //binary_data_offset += len(binary_data_count) + len(binary_data)
    }
    fp.write("</Points>\n");
    //
    fp.write("<Cells>\n");
    fp.write(" <DataArray type=\"Int32\" Name=\"connectivity\"");
    if (binary_format) {
        //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        //binary_data = r""
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (k; 0 .. nkc) {
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
                if (two_D) {
                    auto ids = [vtx_id[i][j][k], vtx_id[i+1][j][k],
				vtx_id[i+1][j+1][k], vtx_id[i][j+1][k]];
                    if (binary_format) {
                        //binary_data += struct.pack('> i i i i', ids[0], ids[1], ids[2], ids[3])
                    } else {
                        fp.writef(" %d %d %d %d\n", ids[0], ids[1], ids[2], ids[3]);
		    }
                } else {
                    auto ids = [vtx_id[i][j][k], vtx_id[i+1][j][k], 
				vtx_id[i+1][j+1][k], vtx_id[i][j+1][k],
				vtx_id[i][j][k+1], vtx_id[i+1][j][k+1], 
				vtx_id[i+1][j+1][k+1], vtx_id[i][j+1][k+1]];
                    if (binary_format) {
                        //binary_data += struct.pack('> i i i i i i i i', ids[0], ids[1], ids[2],
                        //                           ids[3], ids[4], ids[5], ids[6], ids[7])
                    } else {
                        fp.writef(" %d %d %d %d %d %d %d %d\n", ids[0], ids[1], ids[2],
				  ids[3], ids[4], ids[5], ids[6], ids[7]);
		    }
		}
	    } // end foreach i
	} // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        //binary_data_string += binary_data_count
        //binary_data_string += binary_data
        //binary_data_offset += len(binary_data_count) + len(binary_data)
    }
    //
    fp.write(" <DataArray type=\"Int32\" Name=\"offsets\"");
    if (binary_format) {
        //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        //binary_data = r""
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    // Since all of the point-lists are concatenated, these offsets into the connectivity
    // array specify the end of each cell.
    foreach (k; 0 .. nkc) {
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
		size_t conn_offset;
                if (two_D) {
                    conn_offset = 4*(1+i+j*nic);
                } else {
                    conn_offset = 8*(1+i+j*nic+k*(nic*njc));
		}
		if (binary_format) {
                    //binary_data += struct.pack('> i', conn_offset)
                } else {
                    fp.writef(" %d\n", conn_offset);
		}
	    } // end foreach i
	} // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        //binary_data_string += binary_data_count
        //binary_data_string += binary_data
        //binary_data_offset += len(binary_data_count) + len(binary_data)
    }
    //
    fp.write(" <DataArray type=\"UInt8\" Name=\"types\"");
    if (binary_format) {
        //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        //binary_data = r""
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    int type_value;
    if (two_D) {
        type_value = 9; // VTK_QUAD
    } else {
        type_value = 12; // VTK_HEXAHEDRON
    }
    foreach (k; 0 .. nkc) {
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
                if (binary_format) {
                    //binary_data += struct.pack('> B', type_value)
                } else {
                    fp.writef(" %d\n", type_value);
		}
	    } // end foreach i
	} // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        //binary_data_string += binary_data_count
        //binary_data_string += binary_data
        //binary_data_offset += len(binary_data_count) + len(binary_data)
    }
    fp.write("</Cells>\n");
    //
    fp.write("<CellData>\n");
    // Write variables from the dictionary.
    foreach (var; flow.variableNames) {
        fp.writef(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"", var);
        if (binary_format) {
            //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
            //binary_data = r""
        } else {
            fp.write(" format=\"ascii\">\n");
	}
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    if (binary_format) {
                        //binary_data += struct.pack('> f', uflowz(flow.data[var][i,j,k]))
                    } else {
                        fp.writef(" %e\n", uflowz(flow[var,i,j,k]));
		    }
		} // end foreach i
	    } // end foreach j
	} // end foreach k
	fp.write(" </DataArray>\n");
        if (binary_format) {
            //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
            //binary_data_string += binary_data_count
            //binary_data_string += binary_data
            //binary_data_offset += len(binary_data_count) + len(binary_data)
	}
    } // end foreach var
    //
    // Write the special variables:
    // i.e. variables constructed from those in the dictionary.
    fp.write(" <DataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
    if (binary_format) {
        //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
        //binary_data = r""
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (k; 0 .. nkc) {
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
                float x = uflowz(flow["vel.x",i,j,k]);
                float y = uflowz(flow["vel.y",i,j,k]);
                float z = uflowz(flow["vel.z",i,j,k]);
                if (binary_format) {
                    //binary_data += struct.pack('> f f f', x, y, z)
                } else {
                    fp.writef(" %e %e %e\n", x, y, z);
		}
	    } // end foreach i
	} // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
        //binary_data_string += binary_data_count
        //binary_data_string += binary_data
        //binary_data_offset += len(binary_data_count) + len(binary_data)
    }
    //
    if (canFind(flow.variableNames, "c.x")) {
	fp.write(" <DataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
        if (binary_format) {
            //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
            //binary_data = r""
        } else {
            fp.write(" format=\"ascii\">\n");
	}
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
		    float x = uflowz(flow["c.x",i,j,k]);
		    float y = uflowz(flow["c.y",i,j,k]);
		    float z = uflowz(flow["c.z",i,j,k]);
		    if (binary_format) {
			//binary_data += struct.pack('> f f f', x, y, z)
		    } else {
			fp.writef(" %e %e %e\n", x, y, z);
		    }
		} // end foreach i
	    } // end foreach j
	} // end foreach k
	fp.write(" </DataArray>\n");
	if (binary_format) {
            //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
            //binary_data_string += binary_data_count
            //binary_data_string += binary_data
            //binary_data_offset += len(binary_data_count) + len(binary_data)
	}
    } // if canFind c.x
    //
    if (canFind(flow.variableNames, "B.x")) {
	fp.write(" <DataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
        if (binary_format) {
            //fp.write(" format=\"appended\" offset=\"%d\">" % binary_data_offset)
            //binary_data = r""
        } else {
            fp.write(" format=\"ascii\">\n");
	}
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
		    float x = uflowz(flow["B.x",i,j,k]);
		    float y = uflowz(flow["B.y",i,j,k]);
		    float z = uflowz(flow["B.z",i,j,k]);
		    if (binary_format) {
			//binary_data += struct.pack('> f f f', x, y, z)
		    } else {
			fp.writef(" %e %e %e\n", x, y, z);
		    }
		} // end foreach i
	    } // end foreach j
	} // end foreach k
	fp.write(" </DataArray>\n");
	if (binary_format) {
            //binary_data_count = struct.pack('> I', len(binary_data)) # 4-byte count of bytes
            //binary_data_string += binary_data_count
            //binary_data_string += binary_data
            //binary_data_offset += len(binary_data_count) + len(binary_data)
	}
    } // if canFind B.x
    //
    fp.write("</CellData>\n");
    fp.write("</Piece>\n");
    fp.write("</UnstructuredGrid>\n");
    if (binary_format) {
        fp.write("<AppendedData encoding=\"raw\">\n");
        //fp.write('_'+binary_data_string);
        fp.write("</AppendedData>\n");
    }
    fp.write("</VTKFile>\n");
    fp.close();
    return;
} // end write_VTK_XML_unstructured_file()
