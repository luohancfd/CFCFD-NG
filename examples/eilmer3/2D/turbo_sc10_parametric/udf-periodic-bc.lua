-- udf-periodic-bc.lua
-- Lua script for the user-defined periodic BC
--
-- This particular example sets up peroidic boundary conditions
-- for the turbine-blade simulation.
-- When called, this boundary conditions looks up the flow data
-- in a cell that would overlay the ghost cell, 
-- shifted by 1 period in the y-direction.
-- We will assume that the boundary blocks are approximately aligned 
-- with the x,y-axes so that we simply add or subtract the y_period value.
--
-- PJ, 07-Mar-2008
--     03-Sep-2008 port to Elmer3

-- We will remember where we found the appropriate cells.
g1_src_blk = {}; g1_src_i = {}; g1_src_j = {}; g1_src_k = {}
g2_src_blk = {}; g2_src_i = {}; g2_src_j = {}; g2_src_k = {}
y_period = 1.0 -- as set by Hannes and Paul

function ghost_cell_position(xc, yc, xw, yw)
   -- c represents the cell-centre
   -- w represents the wall-interface position
   dx = xc - xw; dy = yc - yw
   return xc - 2*dx, yc - 2*dy
end

function ghost_cell(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   i = args.i; j = args.j; k = args.k
   x = args.x; y = args.y
   indx = k*(nni*nnj) + j*nni + i
   if g1_src_blk[indx] == nil then
      if args.which_boundary == NORTH then
	 -- Search for the cell corresponding to the ghost-cell,
	 -- offset by one period.
	 c = sample_flow(block_id, i, j, k)
	 xg, yg = ghost_cell_position(c.x, c.y, x, y)
	 yg = yg - y_period
	 g1_src_blk[indx], g1_src_i[indx], g1_src_j[indx], g1_src_k[indx] = 
	    locate_cell(xg, yg, 0.0)
	 -- Locate cell corresponding to second ghost cell similarly. 
	 j = j - 1
	 c = sample_flow(block_id, i, j, k)
	 xg, yg = ghost_cell_position(c.x, c.y, x, y)
	 yg = yg - y_period
	 g2_src_blk[indx], g2_src_i[indx], g2_src_j[indx], g2_src_k[indx] = 
	    locate_cell(xg, yg, 0.0)
      elseif args.which_boundary == EAST then
	 print("EAST boundary should not be periodic!")
      elseif args.which_boundary == SOUTH then
	 -- Search for the cell corresponding to the ghost-cell,
	 -- offset by one period.
	 c = sample_flow(block_id, i, j, k)
	 xg, yg = ghost_cell_position(c.x, c.y, x, y)
	 yg = yg + y_period
	 g1_src_blk[indx], g1_src_i[indx], g1_src_j[indx], g1_src_k[indx] = 
	    locate_cell(xg, yg, 0.0)
	 -- Locate cell corresponding to second ghost cell similarly. 
	 j = j + 1
	 c = sample_flow(block_id, i, j, k)
	 xg, yg = ghost_cell_position(c.x, c.y, x, y)
	 yg = yg + y_period
	 g2_src_blk[indx], g2_src_i[indx], g2_src_j[indx], g2_src_k[indx] = 
	    locate_cell(xg, yg, 0.0)
      elseif args.which_boundary == WEST then
	 print("WEST boundary should not be periodic!")
      end
      -- print("indx=", indx, "g1=", g1_src_blk[indx], g1_src_i[indx], g1_src_j[indx],
      --       "g2=", g2_src_blk[indx], g2_src_i[indx], g2_src_j[indx])
   end
   -- On subsequent calls, the array entries should be non-nil so
   -- we can immediately look up the flow data.
   cell1 = sample_flow(g1_src_blk[indx], g1_src_i[indx], g1_src_j[indx], g1_src_k[indx])
   cell2 = sample_flow(g2_src_blk[indx], g2_src_i[indx], g2_src_j[indx], g2_src_k[indx])
   return cell1, cell2
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
