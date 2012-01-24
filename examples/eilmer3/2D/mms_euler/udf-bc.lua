-- Lua script for the south and west boundaries
-- of a Manufactured Solution which
-- treats Euler flow.
--
-- Author: Rowan J. Gollan
-- Date: 04-Feb-2008

M_PI = math.pi
cos = math.cos
sin = math.sin

L = 1.0
gam = 1.4

rho0 = 1.0
rhox = 0.15
rhoy = -0.1

uvel0 = 800.0
uvelx = 50.0
uvely = -30.0

vvel0 = 800.0
vvelx = -75.0
vvely = 40.0
wvel0 = 0.0

press0 = 1.0e5
pressx = 0.2e5
pressy = 0.5e5


function rho_function(x, y)
   rho = rho0 + rhox*sin((M_PI*x)/L) + rhoy*cos((M_PI*y)/(2.0*L))
   return rho;
end
function rho_south_bc(x) return rho_function(x, 0.0) end
function rho_west_bc(y) return rho_function(0.0, y) end


function pressure_function(x, y)
   p = press0 + pressx*cos((2.0*M_PI*x)/L) + pressy*sin((M_PI*y)/L)
   return p
end
function pressure_south_bc(x) return pressure_function(x, 0.0) end
function pressure_west_bc(y) return pressure_function(0.0, y) end


function u_function(x, y)
   u = uvel0 + uvelx*sin((3.0*M_PI*x)/(2.0*L)) + uvely*cos((3.0*M_PI*y)/(5.0*L))
   return u
end
function u_south_bc(x) return u_function(x, 0.0) end
function u_west_bc(y) return u_function(0.0, y) end


function v_function(x, y)
   v = vvel0 + vvelx*cos((M_PI*x)/(2.0*L)) + vvely*sin((2.0*M_PI*y)/(3.0*L))
   return v
end
function v_south_bc(x) return v_function(x, 0.0) end
function v_west_bc(y) return v_function(0.0, y) end


function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains {t, x, y, z, csX, csY, csZ, i, j, k, which_boundary}
   -- Set constant conditions across the whole boundary.
   x = args.x; y = args.y
   ghost = {}
   if args.which_boundary == SOUTH then
      ghost.p = pressure_south_bc(x) -- pressure, Pa
      rho = rho_south_bc(x)          -- density, kg/m^3
      ghost.u = u_south_bc(x)        -- x-velocity, m/s
      ghost.v = v_south_bc(x)        -- y-velocity, m/s
   else
      -- Assumed WEST and that we won't call this
      -- from any other boundary
      ghost.p = pressure_west_bc(y)  -- pressure, Pa
      rho = rho_west_bc(y)           -- density, kg/m^3
      ghost.u = u_west_bc(y)         -- x-velocity, m/s
      ghost.v = v_west_bc(y)         -- y-velocity, m/s
   end
   ghost.w = 0.0
   R = 287.1
   ghost.T = {}
   ghost.T[0] = ghost.p/(rho*R)      -- temperature, K
   ghost.massf = {}     -- mass fractions to be provided as a table
   ghost.massf[0] = 1.0 -- mass fractions are indexed from 0 to nsp-1
   ghost.Tvib = {}   -- vibrational temperatures also indexed from 0
   return ghost, ghost
end

function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- args contains {t, x, y, z, csX, csY, csZ, i, j, k, which_boundary}
   x = args.x; y = args.y
   wall = {}
   if args.which_boundary == SOUTH then
      wall.u = u_south_bc(x) 
      wall.v = v_south_bc(x)
      p = pressure_south_bc(x)
      rho = rho_south_bc(x)
   else
      wall.u = u_west_bc(y) 
      wall.v = v_west_bc(y)
      p = pressure_west_bc(y)
      rho = rho_west_bc(y)
   end
   wall.w = 0.0
   R = 287.1
   wall.T = {}
   wall.T[0] = p/(rho*R)
   wall.massf = {}
   wall.massf[0] = 1.0
   return wall
end
