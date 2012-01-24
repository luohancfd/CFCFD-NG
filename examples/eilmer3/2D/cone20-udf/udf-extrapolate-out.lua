-- udf-extrapolate-out.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example is defining the supersonic outflow
-- for the cone20 test case.


function ghost_cell(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   --
   -- Sample the flow field at the current cell 
   -- which is beside the boundary.
   cell = sample_flow(block_id, args.i, args.j, args.k)
   return cell, cell
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
