-- udf-subsonic-sc10.lua
-- Lua script for the user-defined subsonic inflow for sc10 profile 
-- called by the UserDefinedGhostCell BC.

-- input parameters:
T0 = 300                                   -- total temp in [K]
p0 = 100.0e3                               -- total pressure [Pa]
alpha = math.rad(55)                       -- inflow angle [rad]

-- constants and definitions:
R = 287.0                                  -- gas constant [J/(kg.K)]
g = 1.4                                    -- ratio of specific heats [-]
Cp = g*R/(g-1)                             -- specific-heat, constant volume [J/(kg.K)]
h0 = Cp*T0                                 -- total enthalpy [J/kg]

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   -- Set constant conditions across the whole boundary.

   cell_flow = sample_flow(block_id, args.i, args.j, args.k) -- adjacent cell properties
   vel_sq = cell_flow.u^2+cell_flow.v^2      -- square of inflow velocity [m^2/s^2]
   vel = math.sqrt(vel_sq)                   -- inflow velocity [m/s]
   M_sq = vel_sq/(h0-0.5*vel_sq)/(g-1)       -- square of Mach number [-]
   ratio = 1+0.5*(g-1)*M_sq                  -- T0/T [-]

   ghost = {}
   ghost.p = p0/math.pow(ratio,(g/(g-1)))    -- pressure [Pa]
   ghost.T = {}
   ghost.T[0] = (h0-0.5*vel_sq)/Cp           -- temperature [K]
   ghost.u = vel*math.cos(alpha)             -- x-velocity [m/s]
   ghost.v = vel*math.sin(alpha)             -- y-velocity [m/s]
   ghost.w = 0.0
   ghost.massf = {}                          -- mass fractions to be provided as a table
   ghost.massf[0] = 1.0                      -- mass fractions are indexed from 0 to nsp-1
   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
