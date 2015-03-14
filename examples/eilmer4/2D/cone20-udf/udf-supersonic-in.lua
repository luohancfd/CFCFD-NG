-- udf-supersonic-in.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example is defining the constant supersonic inflow
-- for the cone20 test case.
--
-- RG & PJ 2015-03-14 : ported from PJs eilmer3 example

function ghostCells(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, timeStep, timeLevel, dt, x, y, z, csX, csY, csZ, i, j, k
   -- but we don't happen to us any of them.
   --
   -- Set constant conditions across the whole boundary.
   ghost = {}
   ghost.p = 95.84e3 -- pressure, Pa
   ghost.T = {1103.0}  -- temperatures, K (as a table)
   ghost.massf = {1.0} -- mass fractions to be provided as a table
   ghost.velx = 1000.0  -- x-velocity, m/s
   ghost.vely = 0.0     -- y-velocity, m/s
   ghost.velz = 0.0
   return ghost, ghost
end

