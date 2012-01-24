-- udf-subsonic-bc.lua
-- Lua script for the user-defined subsonic inflow for sc10_3D profile 
-- called by the UserDefinedGhostCell BC.
-- PJ, 03-Nov-2008
-- Modified by M.P. van der Laan Nov-2009 for MP-model
-- simplified by PJ 15-Nov-2009

-- constants and definitions:
omega_ini = 2000.0
p_inflow = 80000.0
T_inflow = 250.0
w_inflow = 400.0

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   -- Set constant conditions across the whole boundary.
   -- We assume that the inflow plane is a (y,z)-plane and that
   -- the velocity vector lies in the cylindrical surface normal
   -- to the inflow plane.

   x = args.x; y = args.y; z = args.z
   r = math.sqrt(x*x + y*y)  -- r is radial coordinate in x,y-plane

   c_theta = omega_ini*r
   c_r = 0.0
   w_x = -y/r*c_theta + x/r*c_r 
   w_y =  x/r*c_theta + y/r*c_r 

   ghost = {}
   ghost.p = p_inflow   	-- Pressure [Pa]
   ghost.T = {}
   ghost.T[0] = T_inflow        -- Temperature [K] 
   ghost.u, ghost.v = w_x, w_y  -- Cartesian velocities
   ghost.w = w_inflow 		-- z-velocity [m/s]
   ghost.massf = {}             -- Mass fractions to be provided as a table
   ghost.massf[0] = 1.0         -- Mass fractions are indexed from 0 to nsp-1  
   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
