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
--     03-Nov-2008 sc10_3D version
-- Max.F
--     14-Jul-2009 modification for simulation of Blade Radial Turbine Nasa 6.02-inch
-- Paul van der Laan
--     Okt-2009 testcase.py (mixing plane model)

-- We will remember where we found the appropriate cells.
g1_src_blk = {}; g1_src_i = {}; g1_src_j = {}; g1_src_k = {}; g1_src_theta = {};
g2_src_blk = {}; g2_src_i = {}; g2_src_j = {}; g2_src_k = {}; g2_src_theta = {};

theta_period = math.rad(360.0/11.0)

-- Here we have assumed that the midspan plane has mapped unstretched
-- to the midspan cylinder.

function ghost_cell_position(xc, yc, zc, xw, yw, zw)
   -- Returns the ghost cell centre position as a reflection of the original cell
   -- in the boundary surface.
   -- c represents the cell-centre
   -- w represents the wall-interface position
   -- We will work all movements in the radial (x,y)-plane.
   theta_w = math.atan2(yw,xw)
   rc = math.sqrt(xc*xc + yc*yc)
   theta_c = math.atan2(yc,xc)
   dtheta = theta_c - theta_w
   theta_ghost = theta_c - 2*dtheta
   return rc*math.cos(theta_ghost), rc*math.sin(theta_ghost), zc
end

function source_cell_position(xc, yc, zc, theta_change)
   -- Computes the cell position in the (x,y)-plane,
   -- theta_change radians around the axis.
   rc = math.sqrt(xc*xc + yc*yc)
   theta = math.atan2(yc,xc)
   new_theta = theta + theta_change
   return rc*math.cos(new_theta), rc*math.sin(new_theta), zc
end

function rotate_velocity_vector(cell, theta_old, theta_new)
   -- Rotates the (x,y)-plane velocity vector.
   -- First, transforms into radial and tangential components, then back.
   vx = cell.u; vy = cell.v

   vt = -vx * math.sin(theta_old) + vy * math.cos(theta_old)
   vr = vx * math.cos(theta_old) + vy * math.sin(theta_old)
   cell.u = -vt * math.sin(theta_new) + vr * math.cos(theta_new)
   cell.v = vt * math.cos(theta_new) + vr * math.sin(theta_new)
   return cell.u, cell.v
end

function ghost_cell(args)
   -- Function that returns the flow state for a ghost cell
   -- for use in the inviscid flux calculations.
   i = args.i; j = args.j; k = args.k
   x = args.x; y = args.y; z = args.z
   indx = j*nnk + k
   if g1_src_blk[indx] == nil then  --When a new block is searching for it's ghost cell, g1_src_blk[indx] = nil again
      if args.which_boundary == NORTH then
	 c = sample_flow(block_id, i, j, k)
  	 xg, yg, zg = ghost_cell_position(c.x, c.y, c.z, x, y, z)
	 xs, ys, zs = source_cell_position(xg, yg, zg, theta_period)
	 g1_src_blk[indx], g1_src_i[indx], g1_src_j[indx], g1_src_k[indx]
	    = locate_cell(xs, ys, zs)
	 g1_src_theta[indx] = math.atan2(ys,xs)
	 -- Locate cell corresponding to second ghost cell similarly. 
	 j = j - 1
	 c = sample_flow(block_id, i, j, k)
	 xg, yg, zg= ghost_cell_position(c.x, c.y, c.z, x, y, z)
	 xs, ys, zs = source_cell_position(xg, yg, zg, theta_period)
	 g2_src_blk[indx], g2_src_i[indx], g2_src_j[indx], g2_src_k[indx]
	    = locate_cell(xs, ys, zs)
         g2_src_theta[indx] = math.atan2(ys,xs)
      elseif args.which_boundary == EAST then
	 print("EAST boundary should not be periodic!")
      elseif args.which_boundary == SOUTH then
      --elseif args.which_boundary == BOTTOM then
	 -- Search for the cell corresponding to the ghost-cell,
	 -- offset by one period.
	 c = sample_flow(block_id, i, j, k)
	 xg, yg, zg = ghost_cell_position(c.x, c.y, c.z, x, y, z)
	 xs, ys, zs = source_cell_position(xg, yg, zg, -theta_period)
	 g1_src_blk[indx], g1_src_i[indx], g1_src_j[indx], g1_src_k[indx] 
	    = locate_cell(xs, ys, zs)
 	 g1_src_theta[indx] = math.atan2(ys,xs)
	 -- Locate cell corresponding to second ghost cell similarly. 
	 j = j + 1
	 c = sample_flow(block_id, i, j, k)
	 xg, yg, zg = ghost_cell_position(c.x, c.y, c.z, x, y, z)
	 xs, ys, zs = source_cell_position(xg, yg, zg, -theta_period)
	 g2_src_blk[indx], g2_src_i[indx], g2_src_j[indx], g2_src_k[indx]
	    = locate_cell(xs, ys, zs)
 	 g2_src_theta[indx] = math.atan2(ys,xs)
      elseif args.which_boundary == WEST then
	 print("WEST boundary should not be periodic!")
      end
   end
   -- On subsequent calls, the array entries should be non-nil so
   -- we can immediately look up the flow data.
   cell1 = sample_flow(g1_src_blk[indx], i, g1_src_j[indx], g1_src_k[indx])
   cell2 = sample_flow(g2_src_blk[indx], i, g2_src_j[indx], g2_src_k[indx])
   if args.which_boundary == NORTH then
      theta_old = g1_src_theta[indx]
      theta_new = theta_old - theta_period
      cell1.u, cell1.v = rotate_velocity_vector(cell1, theta_old, theta_new)
      theta_old = g2_src_theta[indx]
      theta_new = theta_old - theta_period
      cell2.u, cell2.v = rotate_velocity_vector(cell2, theta_old, theta_new)
   elseif args.which_boundary == SOUTH then
   --elseif args.which_boundary == BOTTOM then
      theta_old = g1_src_theta[indx]
      theta_new = theta_old + theta_period
      cell1.u, cell1.v = rotate_velocity_vector(cell1, theta_old, theta_new)
      theta_old = g2_src_theta[indx]
      theta_new = theta_old + theta_period
      cell2.u, cell2.v = rotate_velocity_vector(cell2, theta_old, theta_new)
   end
   return cell1, cell2
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
