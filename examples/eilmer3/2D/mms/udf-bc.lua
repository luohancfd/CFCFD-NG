-- Lua script for the boundaries of a Manufactured Solution
--
-- Author: Rowan J. Gollan
-- Date: 04-Feb-2008
-- Generalised by PJ, May-June-2011

pi = math.pi
cos = math.cos
sin = math.sin

L = 1.0
gam = 1.4

file = io.open("case.txt", "r")
case = file:read("*n")
file:close()

if case == 1 then
   -- Supersonic flow
   rho0=1.0; rhox=0.15; rhoy=-0.1; rhoxy=0.0; arhox=1.0; arhoy=0.5; arhoxy=0.0;
   u0=800.0; ux=50.0; uy=-30.0; uxy=0.0; aux=1.5; auy=0.6; auxy=0.0;
   v0=800.0; vx=-75.0; vy=40.0; vxy=0.0; avx=0.5; avy=2.0/3; avxy=0.0;
   p0=1.0e5; px=0.2e5; py=0.5e5; pxy=0.0; apx=2.0; apy=1.0; apxy=0.0
end

if case == 2 then
   -- Subsonic flow
   rho0=1.0; rhox=0.1; rhoy=0.15; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
   u0=70.0; ux=4.0; uy=-12.0; uxy=7.0; aux=5.0/3; auy=1.5; auxy=0.6;
   v0=90.0; vx=-20.0; vy=4.0; vxy=-11.0; avx=1.5; avy=1.0; avxy=0.9;
   p0=1.0e5; px=-0.3e5; py=0.2e5; pxy=-0.25e5; apx=1.0; apy=1.25; apxy=0.75
end

w0=0.0

function rho(x, y)
   return rho0 + rhox*sin(arhox*pi*x/L) + rhoy*cos(arhoy*pi*y/L)
          + rhoxy*cos(arhoxy*pi*x*y/(L*L))
end

function u(x, y)
   return u0 + ux*sin(aux*pi*x/L) + uy*cos(auy*pi*y/L)
          + uxy*cos(auxy*pi*x*y/(L*L))
end

function v(x, y)
   return v0 + vx*cos(avx*pi*x/L) + vy*sin((avy*pi*y)/L)
          + vxy*cos(avxy*pi*x*y/(L*L))
end

function p(x, y)
   return p0 + px*cos((apx*pi*x)/L) + py*sin(apy*pi*y/L)
          + pxy*sin(apxy*pi*x*y/(L*L))
end

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains {t, x, y, z, csX, csY, csZ, i, j, k, which_boundary}
   -- Set constant conditions across the whole boundary.
   x = args.x; y = args.y
   ghost = {}
   if args.which_boundary == NORTH then
      y = L
   elseif args.which_boundary == SOUTH then
      y = 0.0
   elseif args.which_boundary == EAST then
      x = L
   else
      -- WEST
      x = 0.0
   end
   ghost.p = p(x, y)        -- pressure, Pa
   ghost_rho = rho(x, y)    -- density, kg/m^3
   ghost.u = u(x, y)        -- x-velocity, m/s
   ghost.v = v(x, y)        -- y-velocity, m/s
   ghost.w = 0.0
   R = 287.1
   ghost.T = {}
   ghost.T[0] = ghost.p/(ghost_rho*R)      -- temperature, K
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
   if args.which_boundary == NORTH then
      y = L
   elseif args.which_boundary == SOUTH then
      y = 0.0
   elseif args.which_boundary == EAST then
      x = L
   else
      -- WEST
      x = 0.0
   end
   wall.u = u(x, y) 
   wall.v = v(x, y)
   wall_p = p(x, y)
   wall_rho = rho(x, y)
   wall.w = 0.0
   R = 287.1
   wall.T = {}
   wall.T[0] = wall_p/(wall_rho*R)
   wall.massf = {}
   wall.massf[0] = 1.0
   return wall
end
