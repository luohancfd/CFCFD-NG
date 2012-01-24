-- udf-vortex-flow.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedBC.
--
-- This particular example defines the inviscid flow field
-- of a compressible vortex.

Rgas  = 287      -- J/kg.K
g     = 1.4      -- ratio of specific heats
-- radial limits of flow domain
r_i   = 1.0      -- metres
r_o   = 1.384
-- Set flow properties ar the inner radius.
p_i   = 100.0e3                  -- Pa
M_i   = 2.25
rho_i = 1.0                      -- kg/m**3
T_i   = p_i / (Rgas * rho_i)     -- K
a_i   = math.sqrt(g * Rgas * T_i)     -- m/s
u_i   = M_i * a_i                -- m/s

-- We'll use a bit of extra information to estimate
-- the locations of the ghost cells.
n = 40
dr = (r_o - r_i) / nnj

print("Set up inviscid vortex")
print("    p_i=", p_i, "M_i=", M_i, "rho_i=", rho_i, 
      "T_i=", T_i, "a_i=", a_i, "u_i=", u_i)

function vortex_flow(r)
   u   = u_i * r_i / r
   t1  = r_i / r
   t2  = 1.0 + 0.5 * (g - 1.0) * M_i * M_i * (1.0 - t1 * t1)
   rho = rho_i * math.pow( t2, 1.0/(g - 1.0) )
   p = p_i * math.pow( rho/rho_i, g )
   T = p / (rho * Rgas)
   return u, p, T
end

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- 
   -- We make an estimate of where the ghost cell is in space and
   -- then compute the vortex flow properties for that point.
   x = args.x
   y = args.y
   r = math.sqrt(x*x + y*y)
   theta = math.atan2(y, x)

   ghost1 = {}
   if which_boundary == NORTH then
      r_ghost1 = r + 0.5*dr
   else
      r_ghost1 = r - 0.5*dr
   end
   speed, p, T = vortex_flow(r_ghost1)
   ghost1.p = p
   ghost1.T = {}  -- temperatures as a table
   ghost1.T[0] = T
   ghost1.u = math.sin(theta) * speed
   ghost1.v = -math.cos(theta) * speed
   ghost1.w = 0.0
   ghost1.massf = {} -- mass fractions to be provided as a table
   ghost1.massf[0] = 1.0 -- mass fractions are indexed from 0 to nsp-1

   ghost2 = {}
   if which_boundary == NORTH then
      r_ghost2 = r + 1.5*dr
   else
      r_ghost2 = r - 1.5*dr
   end
   speed, p, T = vortex_flow(r_ghost2)
   ghost2.p = p
   ghost2.T = {}
   ghost2.T[0] = T
   ghost2.u = math.sin(theta) * speed
   ghost2.v = -math.cos(theta) * speed
   ghost2.w = 0.0
   ghost2.massf = {}
   ghost2.massf[0] = 1.0

   return ghost1, ghost2
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   x = args.x
   y = args.y
   r = math.sqrt(x*x + y*y)
   theta = math.atan2(y, x)
   speed, p, T = vortex(r)

   wall = {}
   wall.u = math.sin(theta) * speed
   wall.v = -math.cos(theta) * speed
   wall.w = 0.0
   wall.T_wall = T
   wall.massf = {}
   wall.massf[0] = 1.0
   return wall
end
