-- udf-supersonic-in.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This particular example is defining the constant supersonic inflow
-- for the cone20 test case.


function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- but we don't happen to us any of them.
   --
   -- Set constant conditions across the whole boundary.
   -- print("Hello from function ghost_cell.")
   ghost = {}
   ghost.p = 95.84e3 -- pressure, Pa
   ghost.T = {}  -- temperatures, K (as a table)
   ghost.T[0] = 1103.0
   ghost.u = 1000.0  -- x-velocity, m/s
   ghost.v = 0.0     -- y-velocity, m/s
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
   -- print("Hello from function interface.")
   face = {}
   face.u = 1000.0
   face.v = 0.0
   face.w = 0.0
   face.T = {[0]=1103.0,}
   face.massf = {[0]=1.0,}
   return face
end


function flux(args)
   -- Function that returns the fluxes of conserved quantities.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   --
   -- Set constant conditions across the whole boundary.
   -- print("Hello from function flux.")
   R = 287          -- gas constant J/(kg.K)
   g = 1.4          -- ratio of specific heats
   Cv = R / (g - 1) -- specific-heat, constant volume
   p = 95.84e3      -- pressure, Pa
   T = 1103.0       -- temperature, K
   rho = p/(R*T)    -- density, kg/m**3
   u = 1000.0       -- x-velocity, m/s
   v = 0.0          -- y-velocity, m/s
   w = 0.0
   massf = {}       -- mass fractions to be provided as a table
   massf[0] = 1.0   -- mass fractions are indexed from 0 to nsp-1
   -- Assemble flux vector
   F = {}
   F.mass = rho * (u*args.csX + v*args.csY) -- kg/s/m**2
   F.momentum_x = p * args.csX + u * F.mass
   F.momentum_y = p * args.csY + v * F.mass
   F.momentum_z = 0.0
   F.total_energy = F.mass * (Cv*T + 0.5*(u*u+v*v) + p/rho)
   F.species = {}
   F.species[0] = F.mass * massf[0]
   F.renergies = {}
   F.renergies[0] = F.mass * (Cv*T)
   return F
end
