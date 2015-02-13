// bc_mapped_cell_exchange.d
//
// Acquire ghost cell flow data from cells that have been individually mapped.
// Peter J. 2014-07-26

import std.conv;

import bc;
import sblock;

class MappedCellExchangeBC: BoundaryCondition {
public:
    // For each ghost cell associated with the boundary.
    // the following data structure stores 4 integers,
    // specifying the mapped-cell block and its ijk-indices,
    int[][] mapped_cells;
    // These data are (i,j,k)-triples indexed by [other_block][other_face]
    // and are used by the distributed-memory mapped-cell copying functions
    // to marshall and send the requested data.
    int[][][] incoming_mapped_cells; 
    int[][][] outgoing_mapped_cells;

    this(int id, int boundary)
    {
	blk_id = id;
	which_boundary = boundary;
	type_code = BCCode.mapped_cell;
	ghost_cell_data_available = true;
	sets_conv_flux_directly = false;
	sets_visc_flux_directly = false;
	emissivity = 0.0;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FullFaceExchangeBC(";
	repr ~= "mapped_cells=" ~ to!string(mapped_cells);
	repr ~= ")";
	return to!string(repr);
    }

    override void apply_convective(double t)
    {
	throw new Error("Not implemented yet.");
    } // end apply_convective

    override void apply_viscous(double t)
    {
	throw new Error("Not implemented yet.");
    }  // end apply_viscous

} // end class MappedCellExchangeBC
