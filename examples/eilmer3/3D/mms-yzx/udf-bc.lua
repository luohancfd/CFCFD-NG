-- Lua script for the boundaries of a Manufactured Solution
--
-- Author: Rowan J. Gollan
-- Date: 04-Feb-2008
-- Generalised by PJ, May-June-2011
-- yzx orientation, 07-Nov-2012

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
   rho0=1.0; rhoy=0.15; rhoz=-0.1; rhoyz=0.0; arhoy=1.0; arhoz=0.5; arhoyz=0.0;
   u0=800.0; uy=50.0; uz=-30.0; uyz=0.0; auy=1.5; auz=0.6; auyz=0.0;
   v0=800.0; vy=-75.0; vz=40.0; vyz=0.0; avy=0.5; avz=2.0/3; avyz=0.0;
   p0=1.0e5; py=0.2e5; pz=0.5e5; pyz=0.0; apy=2.0; apz=1.0; apyz=0.0
end

if case == 2 or case == 4 then
   -- Subsonic flow
   rho0=1.0; rhoy=0.1; rhoz=0.15; rhoyz=0.08; arhoy=0.75; arhoz=1.0; arhoyz=1.25;
   u0=70.0; uy=4.0; uz=-12.0; uyz=7.0; auy=5.0/3; auz=1.5; auyz=0.6;
   v0=90.0; vy=-20.0; vz=4.0; vyz=-11.0; avy=1.5; avz=1.0; avyz=0.9;
   p0=1.0e5; py=-0.3e5; pz=0.2e5; pyz=-0.25e5; apy=1.0; apz=1.25; apyz=0.75
end

w0=0.0

if case == 1 or case == 2 then
   function S(y, z) return 1.0 end
else
   function S(y, z)
      rsq = (y - L/2)^2 + (z - L/2)^2
      return exp(-16.0*rsq/(L*L))
   end
end

function rho(y, z)
   return rho0 + S(y,z)*rhoy*sin(arhoy*pi*y/L) + S(y,z)*rhoz*cos(arhoz*pi*z/L)
          + S(y,z)*rhoyz*cos(arhoyz*pi*y*z/(L*L))
end

function u(y, z)
   return u0 + S(y,z)*uy*sin(auy*pi*y/L) + S(y,z)*uz*cos(auz*pi*z/L)
          + S(y,z)*uyz*cos(auyz*pi*y*z/(L*L))
end

function v(y, z)
   return v0 + S(y,z)*vy*cos(avy*pi*y/L) + S(y,z)*vz*sin((avz*pi*z)/L)
          + S(y,z)*vyz*cos(avyz*pi*y*z/(L*L))
end

function p(y, z)
   return p0 + S(y,z)*py*cos((apy*pi*y)/L) + S(y,z)*pz*sin(apz*pi*z/L)
          + S(y,z)*pyz*sin(apyz*pi*y*z/(L*L))
end

function fill_table(t, y, z)
   t.p = p(y, z)
   t_rho = rho(y, z)
   t.v = u(y, z)
   t.w = v(y, z)
   t.u = 0.0
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
   y = args.y; z = args.z
   i = args.i; j = args.j; k = args.k
   ghost1 = {}
   ghost2 = {}
   if args.which_boundary == NORTH then
      cell = sample_flow(block_id, i, j+1, k)
      ghost1 = fill_table(ghost1, cell.y, cell.z)
      cell = sample_flow(block_id, i, j+2, k)
      ghost2 = fill_table(ghost2, cell.y, cell.z)
   elseif args.which_boundary == SOUTH then
      cell = sample_flow(block_id, i, j-1, k)
      ghost1 = fill_table(ghost1, cell.y, cell.z)
      cell = sample_flow(block_id, i, j-2, k)
      ghost2 = fill_table(ghost2, cell.y, cell.z)
   elseif args.which_boundary == EAST then
      cell = sample_flow(block_id, i+1, j, k)
      ghost1 = fill_table(ghost1, cell.y, cell.z)
      cell = sample_flow(block_id, i+2, j, k)
      ghost2 = fill_table(ghost2, cell.y, cell.z)
   elseif args.which_boundary == WEST then
      cell = sample_flow(block_id, i-1, j, k)
      ghost1 = fill_table(ghost1, cell.y, cell.z)
      cell = sample_flow(block_id, i-2, j, k)
      ghost2 = fill_table(ghost2, cell.y, cell.z)
   elseif args.which_boundary == TOP then
      cell = sample_flow(block_id, i, j, k+1)
      ghost1 = fill_table(ghost1, cell.y, cell.z)
      cell = sample_flow(block_id, i, j, k+2)
      ghost2 = fill_table(ghost2, cell.y, cell.z)
   elseif args.which_boundary == BOTTOM then
      cell = sample_flow(block_id, i, j, k-1)
      ghost1 = fill_table(ghost1, cell.y, cell.z)
      cell = sample_flow(block_id, i, j, k-2)
      ghost2 = fill_table(ghost2, cell.y, cell.z)
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
   y = args.y; z = args.z
   face = {}
   face.u = 0.0
   face.v = u(y, z) 
   face.w = v(y, z)
   face_p = p(y, z)
   face_rho = rho(y, z)
   face.T = {}
   face.T[0] = face_p/(face_rho*R)
   face.massf = {}
   face.massf[0] = 1.0
   return face
end
