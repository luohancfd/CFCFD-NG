// sblock.d
// Class for structured blocks of cells, for use within Eilmer4.
// This is the "classic" block within the mbcns/Eilmer series 
// of flow simulation codes.

// Peter J. 2014-07-20 first cut.

module sblock;

import std.conv;
import std.file;
import std.json;
import std.stdio;
import std.format;
import std.string;
import std.array;
import gzip;
import geom;
import gasmodel;
import globalconfig;
import fvvertex;
import fvinterface;
import fvcell;
import block;

class SBlock: Block {
public:
    size_t nicell;
    size_t njcell;
    size_t nkcell;
    size_t imin, imax;
    size_t jmin, jmax;
    size_t kmin, kmax;
    size_t[] hicell, hjcell, hkcell; // locations of sample cells for history record
    size_t[] micell, mjcell, mkcell; // locations of monitor cells

private:
    // Total number of cells in each direction for this block.
    // these will be used in the array allocation routines.
    size_t _nidim;
    size_t _njdim;
    size_t _nkdim;
    // Most of the data is stored in the following arrays.
    // ctr = cell center values
    // ifi = east-facing face properties and fluxes (unit normal in the i-index direction)
    // ifj = north-facing face properties and fluxes (normal in the j-index direction)
    // ifk = top-facing 
    // vtx = cell vertex values (used for the viscous terms, mostly)
    // sifi, sifj and sifk are secondary-cell faces (whose corner nodes are the
    //                     the primary-cell centres.
    FVCell[] _ctr;
    FVInterface[] _ifi;
    FVInterface[] _ifj;
    FVInterface[] _ifk;
    FVVertex[] _vtx;
    FVInterface[] _sifi;
    FVInterface[] _sifj;
    FVInterface[] _sifk;

public:
    this(int id, in char[] file_name)
    // Need to create blocks in the context of the GlobalConfig.
    {
	this.id = id;
	if ( file_name.length == 0 ) {
	    throw new Error("We need a file from which to read the block parameters.");
	}
	auto text = cast(string) read(file_name);
	auto items = parseJSON(text);
	nicell = to!size_t(items["nicell"].integer);
	njcell = to!size_t(items["njcell"].integer);
	nkcell = to!size_t(items["nkcell"].integer);
	_nidim = nicell + 2 * nghost;
	_njdim = njcell + 2 * nghost;
	// Indices, in each grid direction for the active cells.
	// These limits are inclusive. The mincell and max cell
	// are both within the active set of cells.
	imin = nghost; imax = imin + nicell - 1;
	jmin = nghost; jmax = jmin + njcell - 1;
	if ( GlobalConfig.dimensions == 2 ) {
	    // In 2D simulations, the k range is from 0 to 0 for the
	    // storage arrays of cells and relevant faces.
	    if ( nkcell != 1 ) {
		writeln("Warning: inconsistent dimensions nkcell set to 1 for 2D");
		nkcell = 1;
	    }
	    _nkdim = 1;
	    kmin = 0; kmax = 0;
	} else {
	    // In 3D simulations the k index is just like the i and j indices.
	    _nkdim = nkcell + 2 * nghost;
	    kmin = nghost; kmax = kmin + nkcell - 1;
	}
    } // end constructor

    override string toString() const
    {
	char[] repr;
	repr ~= "SBlock(";
	repr ~= "id=" ~ to!string(id);
	repr ~= ", nicell=" ~ to!string(nicell);
	repr ~= ", njcell=" ~ to!string(njcell);
	repr ~= ", nkcell=" ~ to!string(nkcell);
	repr ~= ")";
	return to!string(repr);
    }

    size_t to_global_index(size_t i, size_t j, size_t k) const
    {
	if ( k < 0 || k >= _nkdim || j < 0 || j >= _njdim || i < 0 || i >= _nidim ) {
	    throw new Error(text("SBlock:to_global_index: index out of bounds for block[", id,
				 "] i=", i, " j=", j, " k=", k, " nidim=", _nidim, 
				 " njdim=", _njdim, " nkdim=", _nkdim));
	}
	return k * (_njdim * _nidim) + j * _nidim + i; 
    }

    size_t[] to_ijk_indices(size_t gid) const
    {
	size_t k = gid / (_njdim * _nidim);
	size_t j = (gid - k * (_njdim * _nidim)) / _nidim;
	size_t i = gid - k * (_njdim * _nidim) - j * _nidim;
	return [i, j, k];
    }

    ref FVCell get_cell(size_t i, size_t j, size_t k=0) { return _ctr[to_global_index(i,j,k)]; }
    ref FVInterface get_ifi(size_t i, size_t j, size_t k=0) { return _ifi[to_global_index(i,j,k)]; }
    ref FVInterface get_ifj(size_t i, size_t j, size_t k=0) { return _ifj[to_global_index(i,j,k)]; }
    ref FVInterface get_ifk(size_t i, size_t j, size_t k=0) { return _ifk[to_global_index(i,j,k)]; }
    ref FVVertex get_vtx(size_t i, size_t j, size_t k=0) { return _vtx[to_global_index(i,j,k)]; }
    ref FVInterface get_sifi(size_t i, size_t j, size_t k=0) { return _sifi[to_global_index(i,j,k)]; }
    ref FVInterface get_sifj(size_t i, size_t j, size_t k=0) { return _sifj[to_global_index(i,j,k)]; }
    ref FVInterface get_sifk(size_t i, size_t j, size_t k=0) { return _sifk[to_global_index(i,j,k)]; }

    void assemble_arrays()
    // We shouldn't be calling this until the essential bits of the GlobalConfig
    // have been set up.
    {
	if ( GlobalConfig.verbosity_level >= 2 ) 
	    writefln("assemble_arrays(): Begin for block %d", id);
	auto gm = GlobalConfig.gmodel;
	// Check for obvious errors.
	if ( _nidim <= 0 || _njdim <= 0 || _nkdim <= 0 ) {
	    throw new Error(text("resize_arrays(): invalid dimensions nidim=",
				 _nidim, " njdim=", _njdim, " nkdim=", _nkdim));
	}
	size_t ntot = _nidim * _njdim * _nkdim;
	try {
	    // Create the cell and interface objects for the entire block.
	    foreach (gid; 0 .. ntot) {
		_ctr ~= new FVCell(gm); _ctr[gid].id = gid;
		auto ijk = to_ijk_indices(gid);
		if ( ijk[0] >= imin && ijk[0] <= imax && 
		     ijk[1] >= jmin && ijk[1] <= jmax && 
		     ijk[2] >= kmin && ijk[2] <= kmax ) {
		    active_cells ~= _ctr[gid];
		}
		_ifi ~= new FVInterface(gm); _ifi[gid].id = gid;
		_ifj ~= new FVInterface(gm); _ifj[gid].id = gid;
		if ( GlobalConfig.dimensions == 3 ) {
		    _ifk ~= new FVInterface(gm); _ifk[gid].id = gid;
		}
		_vtx ~= new FVVertex(gm); _vtx[gid].id = gid;
		_sifi ~= new FVInterface(gm); _sifi[gid].id = gid;
		_sifj ~= new FVInterface(gm); _sifj[gid].id = gid;
		if ( GlobalConfig.dimensions == 3 ) {
		    _sifk ~= new FVInterface(gm); _sifk[gid].id = gid;
		}
	    } // gid loop
	} catch (Error err) {
	    writeln("Crapped out while assembling block arrays.");
	    writefln("nicell=%d njcell=%d nkcell=%d", nicell, njcell, nkcell);
	    writefln("nidim=%d njdim=%d nkdim=%d", _nidim, _njdim, _nkdim);
	    writeln("Probably ran out of memory.");
	    writeln("Be a little less ambitious and try a smaller grid next time.");
	    writefln("System message: %s", err.msg);
	    throw new Error("Block.assemble_arrays() failed.");
	}
	if ( GlobalConfig.verbosity_level >= 2 ) {
	    writefln("Done assembling arrays for %d cells.", ntot);
	}
    } // end of assemble_arrays()

    void bind_faces_and_vertices_to_cells()
    // There is a fixed order of faces and vertices for each cell.
    // Refer to fvcore.d
    {
	size_t kstart, kend;
	if ( GlobalConfig.dimensions == 3 ) {
	    kstart = kmin - 1;
	    kend = kmax + 1;
	} else {
	    kstart = 0;
	    kend = 0;
	}
	// With these ranges, we also do the first layer of ghost cells.
	for ( size_t k = kstart; k <= kend; ++k ) {
	    for ( size_t j = jmin-1; j <= jmax+1; ++j ) {
		for ( size_t i = imin-1; i <= imax+1; ++i ) {
		    FVCell cell = get_cell(i,j,k);
		    cell.iface ~= get_ifj(i,j+1,k); // north
		    cell.iface ~= get_ifi(i+1,j,k); // east
		    cell.iface ~= get_ifj(i,j,k); // south
		    cell.iface ~= get_ifi(i,j,k); // west
		    cell.vtx ~= get_vtx(i,j,k);
		    cell.vtx ~= get_vtx(i+1,j,k);
		    cell.vtx ~= get_vtx(i+1,j+1,k);
		    cell.vtx ~= get_vtx(i,j+1,k);
		    if ( GlobalConfig.dimensions == 3 ) {
			cell.iface ~= get_ifk(i,j,k+1); // top
			cell.iface ~= get_ifk(i,j,k); // bottom
			cell.vtx ~= get_vtx(i,j,k+1);
			cell.vtx ~= get_vtx(i+1,j,k+1);
			cell.vtx ~= get_vtx(i+1,j+1,k+1);
			cell.vtx ~= get_vtx(i,j+1,k+1);
		    } // end if
		} // for i
	    } // for j
	} // for k
    } // end bind_faces_and_vertices_to_cells()

    // to be ported from block.cxx
    void identify_reaction_zones(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void identify_turbulent_zones(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void clear_fluxes_of_conserved_quantities()
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    int count_invalid_cells(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void init_residuals()
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void compute_residuals(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void determine_time_step_size()
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void detect_shock_points()
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    // to be ported from block_geometry.cxx
    void compute_primary_cell_geometric_data(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void compute_distance_to_nearest_wall_for_all_cells(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void compute_secondary_cell_geometric_data(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void calc_volumes_2D(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void secondary_areas_2D(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void calc_faces_2D(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void calc_ghost_cell_geom_2D(int gtl)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    // to be ported from block_io.cxx
    void read_grid(string filename, size_t gtl=0)
    // Read the grid vertices from a gzip file.
    {
	size_t nivtx, njvtx, nkvtx;
	double x, y, z;
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("read_grid(): Start block ", id);
	}
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &nivtx, &njvtx, &nkvtx);
	if ( GlobalConfig.dimensions == 3 ) {
	    if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx-1 != nkcell ) {
		throw new Error(text("For block[", id, "] we have a mismatch in 3D grid size.",
                                     " Have read nivtx=", nivtx, " njvtx=", njvtx,
				     " nkvtx=", nkvtx));
	    }
	    for ( size_t k = kmin; k <= kmax+1; ++k ) {
		for ( size_t j = jmin; j <= jmax+1; ++j ) {
		    for ( size_t i = imin; i <= imax+1; ++i ) {
			line = byLine.front; byLine.popFront();
			// Note that the line starts with whitespace.
			formattedRead(line, " %g %g %g", &x, &y, &z);
			auto vtx = get_vtx(i,j,k);
			vtx.pos[gtl].refx = x;
			vtx.pos[gtl].refy = y;
			vtx.pos[gtl].refz = z;
		    } // for i
		} // for j
	    } // for k
	} else { // 2D case
	    if ( nivtx-1 != nicell || njvtx-1 != njcell || nkvtx != 1 ) {
		throw new Error(text("For block[", id, "] we have a mismatch in 2D grid size.",
				     " Have read nivtx=", nivtx, " njvtx=", njvtx,
				     " nkvtx=", nkvtx));
	    }
	    for ( size_t j = jmin; j <= jmax+1; ++j ) {
		for ( size_t i = imin; i <= imax+1; ++i ) {
		    line = byLine.front; byLine.popFront();
		    // Note that the line starts with whitespace.
		    formattedRead(line, " %g %g", &x, &y);
		    auto vtx = get_vtx(i,j);
		    vtx.pos[gtl].refx = x;
		    vtx.pos[gtl].refy = y;
		    vtx.pos[gtl].refz = 0.0;
		} // for i
	    } // for j
	}
    } // end read_grid()

    void write_grid(string filename, double sim_time, size_t gtl=0)
    {
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("write_grid(): Start block ", id);
	}
	size_t kmaxrange;
	auto outfile = new GzipOut(filename);
	auto writer = appender!string();
	if ( GlobalConfig.dimensions == 3 ) {
	    formattedWrite(writer, "%d %d %d  # ni nj nk\n", nicell+1, njcell+1, nkcell+1);
	    kmaxrange = kmax + 1;
	} else { // 2D case
	    formattedWrite(writer, "%d %d %d  # ni nj nk\n", nicell+1, njcell+1, nkcell);
	    kmaxrange = kmax;
	}
	outfile.compress(writer.data);
	for ( size_t k = kmin; k <= kmaxrange; ++k ) {
	    for ( size_t j = jmin; j <= jmax+1; ++j ) {
		for ( size_t i = imin; i <= imax+1; ++i ) {
		    auto vtx = get_vtx(i,j,k);
		    writer = appender!string();
		    formattedWrite(writer, "%20.12e %20.12e %20.12e\n", vtx.pos[gtl].x,
				   vtx.pos[gtl].y, vtx.pos[gtl].z);
		    outfile.compress(writer.data);
		} // for i
	    } // for j
	} // for k
	outfile.finish();
    } // end write_grid()

    double read_solution(string filename)
    // Note that the position data is read into grid-time-level 0
    // by scan_values_from_string(). 
    {
	size_t ni, nj, nk;
	double sim_time;
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("read_solution(): Start block ", id);
	}
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, " %g", &sim_time);
	line = byLine.front; byLine.popFront();
	// ignore second line; it should be just the names of the variables
	// [TODO] We should test the incoming strings against the current variable names.
	line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &ni, &nj, &nk);
	if ( ni != nicell || nj != njcell || 
	     nk != ((GlobalConfig.dimensions == 3) ? nkcell : 1) ) {
	    throw new Error(text("For block[", id, "] we have a mismatch in solution size.",
				 " Have read ni=", ni, " nj=", nj, " nk=", nk));
	}	
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t i = imin; i <= imax; ++i ) {
		    line = byLine.front; byLine.popFront();
		    get_cell(i,j,k).scan_values_from_string(line);
		} // for i
	    } // for j
	} // for k
	return sim_time;
    }

    // Returns sim_time from file.

    void write_solution(string filename, double sim_time)
    // Write the flow solution (i.e. the primary variables at the cell centers)
    // for a single block.
    // This is almost Tecplot POINT format.
    {
	if ( GlobalConfig.verbosity_level >= 1 && id == 0 ) {
	    writeln("write_solution(): Start block ", id);
	}
	auto outfile = new GzipOut(filename);
	auto writer = appender!string();
	formattedWrite(writer, "%20.12e\n", sim_time);
	outfile.compress(writer.data);
	writer = appender!string();
	foreach(varname; variable_list_for_cell()) {
	    formattedWrite(writer, " \"%s\"", varname);
	}
	formattedWrite(writer, "\n");
	outfile.compress(writer.data);
	writer = appender!string();
	formattedWrite(writer, "%d %d %d\n", nicell, njcell, nkcell);
	outfile.compress(writer.data);
	for ( size_t k = kmin; k <= kmax; ++k ) {
	    for ( size_t j = jmin; j <= jmax; ++j ) {
		for ( size_t i = imin; i <= imax; ++i ) {
		    outfile.compress(" " ~ get_cell(i,j,k).write_values_to_string() ~ "\n");
		} // for i
	    } // for j
	} // for k
	outfile.finish();
    }

    void write_history(string filename, double sim_time, bool write_header=false)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    // to be ported from invs.cxx
    void inviscid_flux()
    {
	throw new Error("[TODO] Not implemented yet.");
    }

} // end class SBlock


int find_nearest_cell(in Vector3 point, 
		      ref size_t jb_near, ref size_t i_near, ref size_t j_near, ref size_t k_near,
		      int gtl)
{
    int result_flag = 0;
    throw new Error("Not implemented yet.");
    return result_flag;
}

int locate_cell(in Vector3 point,
		ref size_t jb_found, ref size_t i_found, ref size_t j_found, ref size_t k_found,
		int gtl)
{
    int result_flag = 0;
    throw new Error("Not implemented yet.");
    return result_flag;
}



/** Indexing of the data in 2D.
 *
 * \verbatim
 * The following figure shows cell [i,j] and its associated
 * vertices and faces. 
 * (New arrangement, planned August 2006, implemented Nov 2006)
 *
 *
 *
 *     Vertex 3         North face           Vertex 2 
 *   vtx[i,j+1]         ifj[i,j+1]           vtx[i+1,j+1]
 *             +--------------x--------------+
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *   West      |         cell center         |  East 
 *   face      |          ctr[i,j]           |  face
 *   ifi[i,j]  x              o              x  ifi[i+1,j]
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             |                             |
 *             +--------------x--------------+
 *     Vertex 0           South face         Vertex 1
 *     vtx[i,j]           ifj[i,j]           vtx[i+1,j]
 *
 *
 * Thus...
 * ----
 * Active cells are indexed as ctr[i][i], where
 * imin <= i <= imax, jmin <= j <= jmax.
 *
 * Active east-facing interfaces are indexed as ifi[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax.
 *
 * Active north-facing interfaces are indexed as ifj[i][j], where
 * imin <= i <= imax, jmin <= j <= jmax+1.
 *
 * Active vertices are indexed as vtx[i][j], where
 * imin <= i <= imax+1, jmin <= j <= jmax+1.
 *
 * Space for ghost cells is available outside these ranges.
 *
 * Indexing for the 3D data -- see page 8 in 3D CFD workbook
 * \endverbatim
 */
