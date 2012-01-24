-- diaphragm-test.lua
-- Lua script for the user-defined functions 
-- called by the AdjacentPlusUDFBC.
--
-- This particular example partially closes-off 
-- the interblock connection for the cone20 test case.
-- For y > 0.95, use slip-wall conditions along both sides
-- of the interface, else just let the exchange data be used.
--
-- Maybe the arrangement could be used as the basis for 
-- a slowly-opening diaphragm in the X-tube simulations.
--
-- Adapted from udf-slip-wall.lua
-- PJ, 28-Apr-2009

function reflect_normal_velocity(ux, vy, cosX, cosY)
   -- Copied from cns_bc.h.
   un = ux * cosX + vy * cosY;     -- Normal velocity
   vt = -ux * cosY + vy * cosX;    -- Tangential velocity
   un = -un;                       -- Reflect normal component
   ux = un * cosX - vt * cosY;     -- Back to Cartesian coords
   vy = un * cosY + vt * cosX;
   return ux, vy
end

function ghost_cell(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   i = args.i; j = args.j; k = args.k
   cell1 = sample_flow(block_id, i, j, k)
   cell1.u, cell1.v = reflect_normal_velocity(cell1.u, cell1.v, args.csX, args.csY)
   if args.which_boundary == NORTH then
      j = j - 1
   elseif args.which_boundary == EAST then
      i = i - 1
   elseif args.which_boundary == SOUTH then
      j = j + 1
   elseif args.which_boundary == WEST then
      i = i + 1
   end
   cell2 = sample_flow(block_id, i, j, k)
   cell2.u, cell2.v = reflect_normal_velocity(cell2.u, cell2.v, args.csX, args.csY)
   if args.y > 0.95 then
      return cell1, cell2
   else
      return nil, nil
   end
end


function zero_normal_velocity(ux, vy, cosX, cosY)
   -- Just the interesting bits from reflect_normal_velocity().
   vt = -ux * cosY + vy * cosX;    -- Tangential velocity
   ux = -vt * cosY;                -- Back to Cartesian coords
   vy =  vt * cosX;                -- just tangential component
   return ux, vy
end

function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   cell = sample_flow(block_id, args.i, args.j, args.k)
   cell.u, cell.v = zero_normal_velocity(cell.u, cell.v, args.csX, args.csY)
   if args.y > 0.95 then
      return cell
   else
      return nil
   end
end
