// bc/full_face_exchange.d
//
// Cell-to-cell exchange boundary condition.
// Peter J. 2014-07-26

import std.conv;
import std.stdio;

import globalconfig;
import globaldata;
import fvcore;
import fvcell;
import bc;
import sblock;

class FullFaceExchangeBC: BoundaryCondition {
public:
    int neighbour_block;
    int neighbour_face;
    int neighbour_orientation;

    this(int id, int boundary, 
	 int other_block, int other_face, int orient)
    {
	blk_id = id;
	which_boundary = boundary;
	type_code = BCCode.full_face_exchange;
	is_wall = false;
	ghost_cell_data_available = true;
	sets_conv_flux_directly = false;
	sets_visc_flux_directly = false;
	emissivity = 0.0;
	neighbour_block = other_block;
	neighbour_face = other_face;
	neighbour_orientation = orient;
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "FullFaceExchangeBC(";
	repr ~= "neighbour_block=" ~ to!string(neighbour_block);
	repr ~= ", neighbour_face=" ~ to!string(neighbour_face);
	repr ~= ", neighbour_orientation=" ~ to!string(neighbour_orientation);
	repr ~= ")";
	return to!string(repr);
    }

    override void apply_convective(double t)
    {
	auto blk = allBlocks[blk_id];
	blk.copy_into_ghost_cells(which_boundary, 
				  allBlocks[neighbour_block],
				  neighbour_face, neighbour_orientation,
				  CopyDataOption.minimal_flow, true);
    }

    override void apply_viscous(double t)
    {
	// do nothing
    }

    override void do_copy_into_boundary()
    {
	// TODO Check me!  This is a work-around.
	// We should be able to directly reference the BCs block as blk.
	auto blk = allBlocks[blk_id];
	blk.copy_into_ghost_cells(which_boundary, 
				  allBlocks[neighbour_block],
				  neighbour_face, neighbour_orientation,
				  CopyDataOption.all, true);
    }
} // end class FullFaceExchangeBC

