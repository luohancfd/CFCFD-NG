-- Lua script for the boundaries of a Manufactured Solution
--
-- Author: Rowan J. Gollan
-- Date: 04-Feb-2008
-- Generalised by PJ, May-June-2011

pi = math.pi
cos = math.cos
sin = math.sin
exp = math.exp

L = 1.0
R = 287.0
gam = 1.4

file = io.open("case.txt", "r")
case = file:read("*n")
file:close()

if case == 1 or case == 3 then
   -- Supersonic flow
   rho0=1.0; rhox=0.15; rhoy=-0.1; rhoxy=0.0; arhox=1.0; arhoy=0.5; arhoxy=0.0;
   u0=800.0; ux=50.0; uy=-30.0; uxy=0.0; aux=1.5; auy=0.6; auxy=0.0;
   v0=800.0; vx=-75.0; vy=40.0; vxy=0.0; avx=0.5; avy=2.0/3; avxy=0.0;
   p0=1.0e5; px=0.2e5; py=0.5e5; pxy=0.0; apx=2.0; apy=1.0; apxy=0.0
end

if case == 2 or case == 4 then
   -- Subsonic flow
   rho0=1.0; rhox=0.1; rhoy=0.15; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
   u0=70.0; ux=4.0; uy=-12.0; uxy=7.0; aux=5.0/3; auy=1.5; auxy=0.6;
   v0=90.0; vx=-20.0; vy=4.0; vxy=-11.0; avx=1.5; avy=1.0; avxy=0.9;
   p0=1.0e5; px=-0.3e5; py=0.2e5; pxy=-0.25e5; apx=1.0; apy=1.25; apxy=0.75
end

w0=0.0

if case == 1 or case == 2 then
   function S(x, y) return 1.0 end
else
   function S(x, y)
      rsq = (x - L/2)^2 + (y - L/2)^2
      return exp(-16.0*rsq/(L*L))
   end
end

function rho(x, y)
   return rho0 + S(x,y)*rhox*sin(arhox*pi*x/L) + S(x,y)*rhoy*cos(arhoy*pi*y/L)
          + S(x,y)*rhoxy*cos(arhoxy*pi*x*y/(L*L))
end

function u(x, y)
   return u0 + S(x,y)*ux*sin(aux*pi*x/L) + S(x,y)*uy*cos(auy*pi*y/L)
          + S(x,y)*uxy*cos(auxy*pi*x*y/(L*L))
end

function v(x, y)
   return v0 + S(x,y)*vx*cos(avx*pi*x/L) + S(x,y)*vy*sin((avy*pi*y)/L)
          + S(x,y)*vxy*cos(avxy*pi*x*y/(L*L))
end

function p(x, y)
   return p0 + S(x,y)*px*cos((apx*pi*x)/L) + S(x,y)*py*sin(apy*pi*y/L)
          + S(x,y)*pxy*sin(apxy*pi*x*y/(L*L))
end

function fill_table(t, x, y)
   t.p = p(x, y)
   t_rho = rho(x, y)
   t.u = u(x, y)
   t.v = v(x, y)
   t.w = 0.0
   t.T = {}
   t.T[0] = t.p/(t_rho*R)      -- temperature, K
   t.massf = {}     -- mass fractions to be provided as a table
   t.massf[0] = 1.0 -- mass fractions are indexed from 0 to nsp-1
   t.Tvib = {}   -- vibrational temperatures also indexed from 0
   return t
end


function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains {t, x, y, z, csX, csY, csZ, i, j, k, which_boundary}
   -- Set constant conditions across the whole boundary.
   x = args.x; y = args.y
   i = args.i; j = args.j; k = args.k
   ghost1 = {}
   ghost2 = {}
   if args.which_boundary == NORTH then
      cell = sample_flow(block_id, i, j+1, k)
      ghost1 = fill_table(ghost1, cell.x, cell.y)
      cell = sample_flow(block_id, i, j+2, k)
      ghost2 = fill_table(ghost2, cell.x, cell.y)
   elseif args.which_boundary == SOUTH then
      cell = sample_flow(block_id, i, j-1, k)
      ghost1 = fill_table(ghost1, cell.x, cell.y)
      cell = sample_flow(block_id, i, j-2, k)
      ghost2 = fill_table(ghost2, cell.x, cell.y)
   elseif args.which_boundary == EAST then
      cell = sample_flow(block_id, i+1, j, k)
      ghost1 = fill_table(ghost1, cell.x, cell.y)
      cell = sample_flow(block_id, i+2, j, k)
      ghost2 = fill_table(ghost2, cell.x, cell.y)
   elseif args.which_boundary == WEST then
      cell = sample_flow(block_id, i-1, j, k)
      ghost1 = fill_table(ghost1, cell.x, cell.y)
      cell = sample_flow(block_id, i-2, j, k)
      ghost2 = fill_table(ghost2, cell.x, cell.y)
   elseif args.which_boundary == TOP then
      cell = sample_flow(block_id, i, j, k+1)
      ghost1 = fill_table(ghost1, cell.x, cell.y)
      cell = sample_flow(block_id, i, j, k+2)
      ghost2 = fill_table(ghost2, cell.x, cell.y)
   elseif args.which_boundary == BOTTOM then
      cell = sample_flow(block_id, i, j, k-1)
      ghost1 = fill_table(ghost1, cell.x, cell.y)
      cell = sample_flow(block_id, i, j, k-2)
      ghost2 = fill_table(ghost2, cell.x, cell.y)
   else
      print("ghost-cell: We just fell through a hole in the floor.")
   end
   return ghost1, ghost2
end

function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- args contains {t, x, y, z, csX, csY, csZ, i, j, k, which_boundary}
   x = args.x; y = args.y
   face = {}
   face.u = u(x, y) 
   face.v = v(x, y)
   face_p = p(x, y)
   face_rho = rho(x, y)
   face.w = 0.0
   face.T = {}
   face.T[0] = face_p/(face_rho*R)
   face.massf = {}
   face.massf[0] = 1.0
   return face
end
