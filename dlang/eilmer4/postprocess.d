/** postprocess.d
 * Eilmer4 compressible-flow simulation code, postprocessing functions.
 *
 * The role of the post-processing functions is just to pick up data
 * from a previously-run simulation and either write plotting files
 * or extract interesting pieces of data.
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
import std.bitmanip;
import std.stdint;
import gzip;
import fileutil;
import geom;
import sgrid;
import gas;
import globalconfig;
import readconfig;
import flowsolution;


void post_process(string plotDir, bool listInfoFlag, string tindxPlot,
		  string addVarsStr, string luaRefSoln,
		  bool vtkxmlFlag, bool binary_format,
		  string outputFileName, string sliceListStr, string probeStr,
		  string normsStr, string regionStr)
{
    read_config_file();
    string jobName = GlobalConfig.base_file_name;
    //
    string[] addVarsList;
    addVarsStr = addVarsStr.strip();
    addVarsStr = removechars(addVarsStr, "\"");
    if (addVarsStr.length > 0) {
	addVarsList = addVarsStr.split(",");
    }
    //
    auto times_dict = readTimesFile(jobName);
    auto tindx_list = times_dict.keys;
    sort(tindx_list);
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
    //
    if (listInfoFlag) {
	writeln("Some information about this simulation.");
	writeln("  nBlocks= ", GlobalConfig.nBlocks);
	writeln("  last tindx= ", tindx_list[$-1]);
	writeln("  Variables:");
	// Dip into the top of a solution file that is likely to be present
	// to get the variable names, as saved by the simulation.
	string fileName = make_file_name!"flow"(jobName, to!int(0), 0);
	auto byLine = new GzipByLine(fileName);
	auto line = byLine.front; byLine.popFront();
	double sim_time;
	formattedRead(line, " %g", &sim_time);
	line = byLine.front; byLine.popFront();
	auto variableNames = line.strip().split();
	foreach (ref var; variableNames) { var = removechars(var, "\""); }
	foreach (i; 0 .. variableNames.length) {
	    writeln(format("%4d %s", i, variableNames[i]));
	}
    } // end if listInfoFlag
    //
    if (vtkxmlFlag) {
	writeln("writing VTK-XML files to directory \"", plotDir, "\"");
	begin_Visit_file(jobName, plotDir, GlobalConfig.nBlocks);
	begin_PVD_file(jobName, plotDir);
	foreach (tindx; tindx_list_to_plot) {
	    writeln("  tindx= ", tindx);
	    auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nBlocks);
	    soln.add_aux_variables(addVarsList);
	    if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
	    write_VTK_XML_files(jobName, plotDir, soln, tindx, binary_format);
	} // foreach tindx
	finish_PVD_file(jobName, plotDir);
    } // end if vtkxml
    //
    if (probeStr.length > 0) {
	writeln("Probing flow solution at specified points.");
	// The output may go to a user-specified file, or stdout.
	File outFile;
	if (outputFileName.length > 0) {
	    outFile = File(outputFileName, "w");
	    writeln("Output will be sent to File: ", outputFileName);
	} else {
	    outFile = stdout;
	}
	probeStr = probeStr.strip();
	probeStr = removechars(probeStr, "\"");
	double[] xp, yp, zp;
	foreach(triple; probeStr.split(";")) {
	    auto items = triple.split(",");
	    xp ~= to!double(items[0]);
	    yp ~= to!double(items[1]);
	    zp ~= to!double(items[2]);
	}
	foreach (tindx; tindx_list_to_plot) {
	    writeln("  tindx= ", tindx);
	    auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nBlocks);
	    soln.add_aux_variables(addVarsList);
	    if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
	    outFile.writeln(soln.flowBlocks[0].variable_names_as_string());
	    foreach (ip; 0 .. xp.length) {
		auto nearest = soln.find_nearest_cell_centre(xp[ip], yp[ip], zp[ip]);
		size_t ib = nearest[0]; size_t i = nearest[1];
		size_t j = nearest[2]; size_t k = nearest[3];
		outFile.writeln(soln.flowBlocks[ib].values_as_string(i,j,k));
	    }
	} // end foreach tindx
    } // end if probeStr
    //
    if (sliceListStr.length > 0) {
	writeln("Extracting slices of the flow solution.");
	// The output may go to a user-specified file, or stdout.
	File outFile;
	if (outputFileName.length > 0) {
	    outFile = File(outputFileName, "w");
	    writeln("Output will be sent to File: ", outputFileName);
	} else {
	    outFile = stdout;
	}
	foreach (tindx; tindx_list_to_plot) {
	    writeln("  tindx= ", tindx);
	    auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nBlocks);
	    soln.add_aux_variables(addVarsList);
	    if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
	    //
	    outFile.writeln(soln.flowBlocks[0].variable_names_as_string());
	    foreach (sliceStr; sliceListStr.split(";")) {
		auto rangeStrings = sliceStr.split(",");
		auto blk_range = decode_range_indices(rangeStrings[0], 0, soln.nBlocks);
		foreach (ib; blk_range[0] .. blk_range[1]) {
		    auto blk = soln.flowBlocks[ib];
		    // We need to do the decode in the context of each block because
		    // the upper limits to the indices are specific to the block.
		    auto i_range = decode_range_indices(rangeStrings[1], 0, blk.nic);
		    auto j_range = decode_range_indices(rangeStrings[2], 0, blk.njc);
		    auto k_range = decode_range_indices(rangeStrings[3], 0, blk.nkc);
		    foreach (k; k_range[0] .. k_range[1]) {
			foreach (j; j_range[0] .. j_range[1]) {
			    foreach (i; i_range[0] .. i_range[1]) {
				outFile.writeln(blk.values_as_string(i,j,k));
			    }
			}
		    }
		} // end foreach ib
	    } // end foreach sliceStr
	} // end foreach tindx
    } // end if sliceListStr
    //
    if (normsStr.length > 0) {
	writeln("Norms for variables.");
	normsStr = normsStr.strip();
	normsStr = removechars(normsStr, "\"");
	foreach (tindx; tindx_list_to_plot) {
	    writeln("  tindx= ", tindx);
	    auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nBlocks);
	    soln.add_aux_variables(addVarsList);
	    if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
	    //
	    foreach (varName; normsStr.split(",")) {
		auto norms = soln.compute_volume_weighted_norms(varName, regionStr);
		write("    variable= ", varName);
		write(" L1= ", norms[0], " L2= ", norms[1]);
		write(" Linf= ", norms[2], " x= ", norms[3], " y= ", norms[4], " z= ", norms[5]);
		write("\n");
	    }
	} // end foreach tindx
    } // end if normsStr
    //
} // end post_process()

//-----------------------------------------------------------------------

size_t[] decode_range_indices(string rangeStr, size_t first, size_t endplus1)
// Decode strings such as "0:$", ":", "0:3",
// On input, first and endplus1 represent the largest, available range.
// Return the pair of numbers that can be used in a foreach loop range.
{
    if (rangeStr == ":") {
	return [first, endplus1];
    }
    if (canFind(rangeStr, ":")) {
	// We have a range specification to pull apart.
	auto items = rangeStr.split(":");
	first = to!size_t(items[0]);
	if (items.length > 1 && items[1] != "$") {
	    // Presume that we have a second integer.
	    size_t new_endplus1 = to!size_t(items[1]);
	    if (new_endplus1 < endplus1) endplus1 = new_endplus1; 
	}
    } else {
	// Presume that we have a single integer.
	first = to!size_t(rangeStr);
	if (first < endplus1) endplus1 = first+1;
    }
    return [first, endplus1];
} // end decode_range_indices()

double[int] readTimesFile(string jobName)
{
    double[int] times_dict;
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
    return times_dict;
} // end readTimesFile()

//-----------------------------------------------------------------------

void begin_Visit_file(string jobName, string plotDir, int nBlocks)
{
    // Will be handy to have a Visit file, also.
    // For each time index, this justs lists the names of the files for individual blocks.
    ensure_directory_is_present(plotDir);
    string fileName = plotDir ~ "/" ~ jobName ~ ".visit";
    auto visitFile = File(fileName, "w");
    visitFile.writef("!NBLOCKS %d\n", nBlocks);
    visitFile.close();
    return;
} // end begin_Visit_file()

void begin_PVD_file(string jobName, string plotDir)
{
    // Will be handy to have a Paraview collection file, also.
    // For each time index, this justs lists the name of the top-level .pvtu file.
    ensure_directory_is_present(plotDir);
    string fileName = plotDir ~ "/" ~ jobName ~ ".pvd";
    auto pvdFile = File(fileName, "w");
    pvdFile.write("<?xml version=\"1.0\"?>\n");
    pvdFile.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    pvdFile.write("<Collection>\n");
    pvdFile.close();
    return;
} // end begin_PVD_file()

void finish_PVD_file(string jobName, string plotDir)
{
    ensure_directory_is_present(plotDir);
    string fileName = plotDir ~ "/" ~ jobName ~ ".pvd";
    auto pvdFile = File(fileName, "a");
    pvdFile.write("</Collection>\n");
    pvdFile.write("</VTKFile>\n");
    pvdFile.close();
    return;
} // end finish_PVD_file()

void write_VTK_XML_files(string jobName, string plotDir,
			 FlowSolution soln, int tindx,
			 bool binary_format=false)
{
    ensure_directory_is_present(plotDir);
    string fileName = plotDir~"/"~jobName~format(".t%04d", tindx)~".pvtu";
    auto pvtuFile = File(fileName, "w");
    pvtuFile.write("<VTKFile type=\"PUnstructuredGrid\">\n");
    pvtuFile.write("<PUnstructuredGrid GhostLevel=\"0\">");
    pvtuFile.write("<PPoints>\n");
    pvtuFile.write(" <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    pvtuFile.write("</PPoints>\n");
    pvtuFile.write("<PCellData>\n");
    foreach (var; soln.flowBlocks[0].variableNames) {
        pvtuFile.writef(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"/>\n", var);
    }
    pvtuFile.write(" <PDataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    if (canFind(soln.flowBlocks[0].variableNames,"c.x")) {
	pvtuFile.write(" <PDataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    }
    if (canFind(soln.flowBlocks[0].variableNames,"B.x")) {
	pvtuFile.write(" <PDataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    }
    pvtuFile.write("</PCellData>\n");
    foreach (jb; 0 .. soln.nBlocks) {
        fileName = jobName ~ format(".b%04d.t%04d.vtu", jb, tindx);
        // We write the short version of the fileName into the pvtu file.
        pvtuFile.writef("<Piece Source=\"%s\"/>\n", fileName);
        // but use the long version to actually open it.
        fileName = plotDir ~ "/" ~ fileName;
        write_VTK_XML_unstructured_file(soln, jb, fileName, binary_format);
    }
    pvtuFile.write("</PUnstructuredGrid>\n");
    pvtuFile.write("</VTKFile>\n");
    pvtuFile.close();
    // Will be handy to have a Visit file, also.
    // This justs lists the names of the files for individual blocks.
    fileName = plotDir ~ "/" ~ jobName ~ ".visit";
    // Note that we append to the visit file for each tindx.
    auto visitFile = File(fileName, "a");
    foreach (jb; 0 .. soln.nBlocks) {
        fileName = jobName ~ format(".b%04d.t%04d.vtu", jb, tindx);
        visitFile.writef("%s\n", fileName);
    }
    visitFile.close();
    // Will be handy to have a Paraview PVD file, also.
    // This justs lists the top-level .pvtu files.
    fileName = plotDir ~ "/" ~ jobName ~ ".pvd";
    // Note that we append to the .pvd file for each tindx.
    auto pvdFile = File(fileName, "a");
    fileName = jobName ~ format(".t%04d.pvtu", tindx);
    pvdFile.writef("<DataSet timestep=\"%e\" group=\"\" part=\"0\" file=\"%s\"/>\n",
		   soln.sim_time, fileName);
    pvdFile.close();
    return;
} // end write_VTK_XML_files()

void write_VTK_XML_unstructured_file(FlowSolution soln, size_t jb, 
				     string fileName, bool binary_format)
// Write the cell-centred flow data from a single block (index jb)
// as an unstructured grid of finite-volume cells.
{
    auto fp = File(fileName, "wb"); // We may be writing some binary data.
    auto flow = soln.flowBlocks[jb];
    auto grid = soln.gridBlocks[jb];
    ubyte[] binary_data_string;
    ubyte[] binary_data;
    int binary_data_offset = 0;
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
	fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
	binary_data.length=0;
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
		    binary_data ~= nativeToBigEndian(x);
		    binary_data ~= nativeToBigEndian(y);
		    binary_data ~= nativeToBigEndian(z);
                } else {
                    fp.writef(" %e %e %e\n", x,y,z);
		}
                vtx_number += 1;
	    }
	}
    }
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Points>\n");
    //
    fp.write("<Cells>\n");
    fp.write(" <DataArray type=\"Int32\" Name=\"connectivity\"");
    if (binary_format) {
	fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
	binary_data.length = 0;
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
			foreach (id; ids) { binary_data ~= nativeToBigEndian(to!int32_t(id)); }
                    } else {
                        fp.writef(" %d %d %d %d\n", ids[0], ids[1], ids[2], ids[3]);
		    }
                } else {
                    auto ids = [vtx_id[i][j][k], vtx_id[i+1][j][k], 
				vtx_id[i+1][j+1][k], vtx_id[i][j+1][k],
				vtx_id[i][j][k+1], vtx_id[i+1][j][k+1], 
				vtx_id[i+1][j+1][k+1], vtx_id[i][j+1][k+1]];
                    if (binary_format) {
			foreach (id; ids) { binary_data ~= nativeToBigEndian(to!int32_t(id)); }
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
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"Int32\" Name=\"offsets\"");
    if (binary_format) {
	fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
	binary_data.length = 0;
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
                    binary_data ~= nativeToBigEndian(to!int32_t(conn_offset));
                } else {
                    fp.writef(" %d\n", conn_offset);
		}
	    } // end foreach i
	} // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"UInt8\" Name=\"types\"");
    if (binary_format) {
        fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
        binary_data.length = 0;
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
                    binary_data ~= nativeToBigEndian(to!uint8_t(type_value));
                } else {
                    fp.writef(" %d\n", type_value);
		}
	    } // end foreach i
	} // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Cells>\n");
    //
    fp.write("<CellData>\n");
    // Write variables from the dictionary.
    foreach (var; flow.variableNames) {
        fp.writef(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"", var);
        if (binary_format) {
	    fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
	    binary_data.length = 0;
        } else {
            fp.write(" format=\"ascii\">\n");
	}
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    if (binary_format) {
                        binary_data ~= nativeToBigEndian(to!float(uflowz(flow[var,i,j,k])));
                    } else {
                        fp.writef(" %e\n", uflowz(flow[var,i,j,k]));
		    }
		} // end foreach i
	    } // end foreach j
	} // end foreach k
	fp.write(" </DataArray>\n");
        if (binary_format) {
	    uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
	    binary_data_string ~= nativeToBigEndian(binary_data_count);
	    binary_data_string ~= binary_data;
	    binary_data_offset += 4 + binary_data.length;
	}
    } // end foreach var
    //
    // Write the special variables:
    // i.e. variables constructed from those in the dictionary.
    fp.write(" <DataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
    if (binary_format) {
        fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
        binary_data.length = 0;
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
		    binary_data ~= nativeToBigEndian(x);
		    binary_data ~= nativeToBigEndian(y);
		    binary_data ~= nativeToBigEndian(z);
                } else {
                    fp.writef(" %e %e %e\n", x, y, z);
		}
	    } // end foreach i
	} // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    if (canFind(flow.variableNames, "c.x")) {
	fp.write(" <DataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
        if (binary_format) {
	    fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
	    binary_data.length = 0;
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
			binary_data ~= nativeToBigEndian(x);
			binary_data ~= nativeToBigEndian(y);
			binary_data ~= nativeToBigEndian(z);
		    } else {
			fp.writef(" %e %e %e\n", x, y, z);
		    }
		} // end foreach i
	    } // end foreach j
	} // end foreach k
	fp.write(" </DataArray>\n");
	if (binary_format) {
	    uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
	    binary_data_string ~= nativeToBigEndian(binary_data_count);
	    binary_data_string ~= binary_data;
	    binary_data_offset += 4 + binary_data.length;
	}
    } // if canFind c.x
    //
    if (canFind(flow.variableNames, "B.x")) {
	fp.write(" <DataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
        if (binary_format) {
	    fp.write(format(" format=\"appended\" offset=\"%d\">", binary_data_offset));
	    binary_data.length = 0;
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
			binary_data ~= nativeToBigEndian(x);
			binary_data ~= nativeToBigEndian(y);
			binary_data ~= nativeToBigEndian(z);
		    } else {
			fp.writef(" %e %e %e\n", x, y, z);
		    }
		} // end foreach i
	    } // end foreach j
	} // end foreach k
	fp.write(" </DataArray>\n");
	if (binary_format) {
	    uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
	    binary_data_string ~= nativeToBigEndian(binary_data_count);
	    binary_data_string ~= binary_data;
	    binary_data_offset += 4 + binary_data.length;
	}
    } // if canFind B.x
    //
    fp.write("</CellData>\n");
    fp.write("</Piece>\n");
    fp.write("</UnstructuredGrid>\n");
    if (binary_format) {
        fp.write("<AppendedData encoding=\"raw\">\n");
        fp.write('_');
        fp.rawWrite(binary_data_string);
        fp.write("</AppendedData>\n");
    }
    fp.write("</VTKFile>\n");
    fp.close();
    return;
} // end write_VTK_XML_unstructured_file()
