-- Lua script for the boundaries of a Manufactured Solution
--
-- Author: Rowan J. Gollan
-- Date: 04-Feb-2008
-- Generalised by PJ, May-June-2011
-- Necessary terms have been added for turbulence model test by
-- Jianyong Wang, 13-July-2016 

pi = math.pi
cos = math.cos
sin = math.sin
exp = math.exp
max = math.max

L = 1.0
R = 287.0
gam = 1.4
PrT = 0.89
Cv = R/(gam-1)
Cp = gam*Cv

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
   -- Subsonic(turbulent) flow
   rho0=1.0; rhox=0.15; rhoy=-0.1; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
   u0=70.0; ux=7.0; uy=-8.0; uxy=5.5; aux=1.5; auy=1.5; auxy=0.6;
   v0=90.0; vx=-5.0; vy=10.0; vxy=-11.0; avx=1.5; avy=1.0; avxy=0.9;
   p0=1.0e5; px=0.2e5; py=0.175e5; pxy=-0.25e5; apx=1.0; apy=1.25; apxy=0.75;
   tke0=780.0; tkex=160.0; tkey=-120.0; tkexy=80.0; atkex=0.65; atkey=0.7; atkexy=0.8;
   omega0=150.0; omegax=-30.0; omegay=22.5; omegaxy=40.0; a_omegax=0.75; a_omegay=0.875; a_omegaxy=0.6
end

w0=0.0


function rho(x, y)
   return rho0 + rhox*cos(arhox*pi*x/L) + rhoy*sin(arhoy*pi*y/L)
          + rhoxy*cos(arhoxy*pi*x*y/(L*L))
end

function u(x, y)
   return u0 + ux*sin(aux*pi*x/L) + uy*cos(auy*pi*y/L)
          + uxy*cos(auxy*pi*x*y/(L*L))
end

function v(x, y)
   return v0 + vx*sin(avx*pi*x/L) + vy*cos((avy*pi*y)/L)
          + vxy*cos(avxy*pi*x*y/(L*L))
end

function p(x, y)
   return p0 + px*cos((apx*pi*x)/L) + py*sin(apy*pi*y/L)
          + pxy*sin(apxy*pi*x*y/(L*L))
end

function tke(x, y)
   return tke0 + tkex*cos((atkex*pi*x)/L) + tkey*sin(atkey*pi*y/L)
          + tkexy*cos(atkexy*pi*x*y/(L*L))
end

function omega(x, y)
   return omega0 + omegax*cos((a_omegax*pi*x)/L) + omegay*sin(a_omegay*pi*y/L)
          + omegaxy*cos(a_omegaxy*pi*x*y/(L*L))
end

function mu_t(x, y)
   return (-120.0*sin(0.7*pi*y) + 160.0*cos(0.65*pi*x) + 80.0*cos(0.8*pi*x*y) + 780.0)*(-0.1*sin(pi*y) + 0.15*cos(0.75*pi*x) + 0.08*cos(1.25*pi*x*y) + 1.0)/max(22.5*sin(0.875*pi*y) 
          - 30.0*cos(0.75*pi*x) + 40.0*cos(0.6*pi*x*y) + 150.0, 4.12478955692153*max(0.0, (9.9*pi*x*sin(0.9*pi*x*y) - 10.0*pi*sin(pi*y))^2 + (-3.3*pi*y*sin(0.6*pi*x*y) + 10.5*pi*cos(1.5*pi*x))^2 
          + (-3.3*pi*x*sin(0.6*pi*x*y) + 9.9*pi*y*sin(0.9*pi*x*y) + 12.0*pi*sin(1.5*pi*y) - 7.5*pi*cos(1.5*pi*x))*(-1.65*pi*x*sin(0.6*pi*x*y) + 4.95*pi*y*sin(0.9*pi*x*y) + 6.0*pi*sin(1.5*pi*y) 
          - 3.75*pi*cos(1.5*pi*x)) - (4.4*pi*x*sin(0.9*pi*x*y) - 1.46666666666667*pi*y*sin(0.6*pi*x*y) - 4.44444444444444*pi*sin(pi*y) + 4.66666666666667*pi*cos(1.5*pi*x))*(9.9*pi*x*sin(0.9*pi*x*y) 
          - 3.3*pi*y*sin(0.6*pi*x*y) - 10.0*pi*sin(pi*y) + 10.5*pi*cos(1.5*pi*x)))^0.5)
end

function fill_table(t, x, y)
   t.p = p(x, y)
   t_rho = rho(x, y)
   t.u = u(x, y)
   t.v = v(x, y)
   t.w = 0.0
   t.tke = tke(x, y)
   t.omega = omega(x, y)
   t.mu_t = mu_t(x, y)
   t.k_t = Cp*t.mu_t/PrT
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
   else -- WEST
      cell = sample_flow(block_id, i-1, j, k)
      ghost1 = fill_table(ghost1, cell.x, cell.y)
      cell = sample_flow(block_id, i-2, j, k)
      ghost2 = fill_table(ghost2, cell.x, cell.y)
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
   face.tke = tke(x, y)
   face.omega = omega(x, y)
   face.T = {}
   face.T[0] = face_p/(face_rho*R)
   face.massf = {}
   face.massf[0] = 1.0
   return face
end
