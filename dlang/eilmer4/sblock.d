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
import geom;
import gasmodel;
import globalconfig;
import fvvertex;
import fvinterface;
import fvcell;
import block;

class SBlock: Block {
public:
    int nicell;
    int njcell;
    int nkcell;
    int imincell, imaxcell;
    int jmincell, jmaxcell;
    int kmincell, kmaxcell;
    int[] hicell, hjcell, hkcell; // locations of sample cells for history record
    int[] micell, mjcell, mkcell; // locations of monitor cells

private:
    // Total number of cells in each direction for this block.
    // these will be used in the array allocation routines.
    int _nidim;
    int _njdim;
    int _nkdim;
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
	nicell = to!int(items["nicell"].integer);
	njcell = to!int(items["njcell"].integer);
	nkcell = to!int(items["nkcell"].integer);
	_nidim = nicell + 2 * nghost;
	_njdim = njcell + 2 * nghost;
	// Indices, in each grid direction for the active cells.
	// These limits are inclusive. The mincell and max cell
	// are both within the active set of cells.
	imincell = nghost; imaxcell = imincell + nicell - 1;
	jmincell = nghost; jmaxcell = jmincell + njcell - 1;
	if ( GlobalConfig.dimensions == 2 ) {
	    // In 2D simulations, the k range is from 0 to 0 for the
	    // storage arrays of cells and relevant faces.
	    if ( nkcell != 1 ) {
		writeln("Warning: inconsistent dimensions nkcell set to 1 for 2D");
		nkcell = 1;
	    }
	    _nkdim = 1;
	    kmincell = 0; kmaxcell = 0;
	} else {
	    // In 3D simulations the k index is just like the i and j indices.
	    _nkdim = nkcell + 2 * nghost;
	    kmincell = nghost; kmaxcell = kmincell + nkcell - 1;
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

    int to_global_index(int i, int j, int k) const
    {
	if ( k < 0 || k >= _nkdim || j < 0 || j >= _njdim || i < 0 || i >= _nidim ) {
	    throw new Error(text("SBlock:to_global_index: index out of bounds for block[", id,
				 "] i=", i, " j=", j, " k=", k, " nidim=", _nidim, 
				 " njdim=", _njdim, " nkdim=", _nkdim));
	}
	return k * (_njdim * _nidim) + j * _nidim + i; 
    }

    int[] to_ijk_indices(int gid) const
    {
	int k = gid / (_njdim * _nidim);
	int j = (gid - k * (_njdim * _nidim)) / _nidim;
	int i = gid - k * (_njdim * _nidim) - j * _nidim;
	return [i, j, k];
    }

    ref FVCell get_cell(int i, int j, int k=0) { return _ctr[to_global_index(i,j,k)]; }
    ref FVInterface get_ifi(int i, int j, int k=0) { return _ifi[to_global_index(i,j,k)]; }
    ref FVInterface get_ifj(int i, int j, int k=0) { return _ifj[to_global_index(i,j,k)]; }
    ref FVInterface get_ifk(int i, int j, int k=0) { return _ifk[to_global_index(i,j,k)]; }
    ref FVVertex get_vtx(int i, int j, int k=0) { return _vtx[to_global_index(i,j,k)]; }
    ref FVInterface get_sifi(int i, int j, int k=0) { return _sifi[to_global_index(i,j,k)]; }
    ref FVInterface get_sifj(int i, int j, int k=0) { return _sifj[to_global_index(i,j,k)]; }
    ref FVInterface get_sifk(int i, int j, int k=0) { return _sifk[to_global_index(i,j,k)]; }

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
	int ntot = _nidim * _njdim * _nkdim;
	try {
	    // Create the cell and interface objects for the entire block.
	    foreach (gid; 0 .. ntot) {
		_ctr ~= new FVCell(gm); _ctr[gid].id = gid;
		auto ijk = to_ijk_indices(gid);
		if ( ijk[0] >= imincell && ijk[0] <= imaxcell && 
		     ijk[1] >= jmincell && ijk[1] <= jmaxcell && 
		     ijk[2] >= kmincell && ijk[2] <= kmaxcell ) {
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
	    writefln("Done assembling for %d cells.", ntot);
	}
    } // end of assemble_arrays()

    void bind_faces_and_vertices_to_cells()
    // There is a fixed order of faces and vertices for each cell.
    // Refer to fvcore.d
    {
	size_t kstart, kend;
	if ( GlobalConfig.dimensions == 3 ) {
	    kstart = kmincell - 1;
	    kend = kmaxcell + 1;
	} else {
	    kstart = 0;
	    kend = 0;
	}
	// With these ranges, we also do the first layer of ghost cells.
	for ( int k = kstart; k <= kend; ++k ) {
	    for ( int j = jmincell-1; j <= jmaxcell+1; ++j ) {
		for ( int i = imincell-1; i <= imaxcell+1; ++i ) {
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
    void read_grid(string filename, bool zip_file=true, int gtl=0)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void write_grid(string filename, double sim_time, bool zip_file=true, int gtl=0)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    double read_solution(string filename, bool zip_file=true, int gtl=0)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    // Returns sim_time from file.

    void write_solution(string filename, double sim_time, bool zip_file=true, int gtl=0)
    {
	throw new Error("[TODO] Not implemented yet.");
    }

    void write_history(string filename, double sim_time, bool write_header=false, int gtl=0)
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
		      ref int jb_near, ref int i_near, ref int j_near, ref int k_near,
		      int gtl)
{
    int result_flag = 0;
    throw new Error("Not implemented yet.");
    return result_flag;
}

int locate_cell(in Vector3 point,
		ref int jb_found, ref int i_found, ref int j_found, ref int k_found,
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
