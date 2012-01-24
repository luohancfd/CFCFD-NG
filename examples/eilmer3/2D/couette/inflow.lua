-- inflow.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example is defining the linear velocity profile
-- for (dodgy, incompressible) Couette flow.

y_max = 0.010
u_max = 100.0
p_inf = 100.0e3
T_inf = 1000.0

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- but we don't happen to us any of them.
   --
   -- Set constant conditions across the whole boundary.
   ghost = {}
   ghost.p = p_inf -- pressure, Pa
   ghost.T = {}  -- temperatures, K (as a table)
   ghost.T[0] = T_inf * args.y / y_max + 300.0
   u = u_max * args.y / y_max
   ghost.u = u
   ghost.v = u
   ghost.w = 0.0
   ghost.massf = {} -- mass fractions to be provided as a table
   ghost.massf[0] = 1.0 -- mass fractions are indexed from 0 to nsp-1
   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- but we don't happen to us any of them.
   wall = {}
   u = u_max * args.y / y_max
   wall.u = u
   wall.v = u
   wall.w = 0.0
   wall.T_wall = T_inf * args.y / y_max + 300.0
   wall.massf = {}
   wall.massf[0] = 1.0
   return wall
end
