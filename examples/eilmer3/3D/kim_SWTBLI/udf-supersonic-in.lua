-- udf-supersonic-in.lua
--
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This script allows user to use an input file
-- that contains profiles of flow variables as
-- inflow conditions for a simulation (non-uniform
-- inflow boundary condition).
-- 
-- It has a added feature that allows for the
-- grid resolution of the input file and the 
-- boundary of the simulation not to be similar,
-- that is, the grid resolutions need not match
-- each other.
--
-- Written by Matthew McGilvray - December 2008
--   for mbcns2
--
-- Modified by Wilson Chan - 03 January 2009
--   for elmer3
--
-- Modified by Wilson Chan - 08 January 2009
--   corrected bug in ghost_cell function
--   that occurs on the last cell if a simulation
--   grid finer than the input grid is used (i.e.
--   when the last simulation cell is larger than
--   the last input cell.
--
-- Modified by Wilson Chan - 09 November 2010
--   Pulled the reading of input data out of the
--   "ghost cell" function
--
------------------------------------------------------------------------------------

-- Read in some Eilmer3 data
--
-- -- mbcns2 --
-- # x(m) y(m) rho(kg/m**3) u(m/s) v(m/s) e(J/kg) p1(Pa) M pitot(Pa) 
-- T(K) mu(Pa.s) mu_t(Pa.s) f[0](fraction) f[1](fraction) f[2](fraction) 
-- f[3](fraction) f[4](fraction) 
--
-- -- elmer3 --
-- # Variables: pos.x pos.y pos.z volume rho vel.x vel.y vel.z p a mu 
-- k[0] mu_t k_t S tke omega massf[0] e[0] T[0]
--
-- -- eilmer3 -- note that this list changes... so be careful
-- # Variables: 1:pos.x 2:pos.y 3:pos.z 4:volume 5:rho 6:vel.x 7:vel.y 
-- 8:vel.z 9:p 10:a 11:mu 12:k[0] 13:mu_t 14:k_t 15:S 16:tke 17:omega 
-- 18:massf[0] 19:e[0] 20:T[0] 21:M_local 22:pitot_p
-- 
-- Note : Setup for elmer3 & eilmer3 only allows for 1 species at the moment.

input_filename = "turb_flat_plate-x-206mm.dat"
io.input(input_filename)

-- Set up lists
x = {}; y = {}; z = {}; volume = {}; rho = {}; u = {}; v = {}; w = {};
p = {}; a = {}; mu = {}; k = {}; mu_t = {}; k_t = {}; S = {}; tke = {};
omega = {}; massf = {}; e = {}; T = {}; M = {}; pitot = {};

n = 1
while true do

   x[n], y[n], z[n], volume[n], rho[n], u[n], v[n], w[n], p[n], a[n], mu[n], k[n], mu_t[n], k_t[n], S[n], tke[n], omega[n], massf[n], e[n], T[n], M[n], pitot[n] = io.read("*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number","*number")

   if x[1]==nil then
      -- Just to be safe ... warn user if script is not reading any data at all.
      print("Warning: Input data seems incorrect! Check format of input data file!")
      break
   end
   if x[n]==nil then
      break
   end

n = n + 1
end


function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.

   n = 1
   while true do

      if args.y < y[1] then
      -- Mainly for situations when first simulation cell has a
      -- lower y-location than the first input cell.
        p_in = p[n];
        u_in = u[n];
        v_in = v[n];
        w_in = w[n];
        T_in = T[n];
        tke_in = tke[n];
        omega_in = omega[n];
        break
      end

      if y[n] == nil and y[n-1] < args.y then
      -- This 'if' statement basically looks for situations
      -- when the search for a correct y location has not located
      -- anything and that is when the final simulation cell has
      -- a higher y location than the final input cell.
        p_in = p[n-1];
        u_in = u[n-1];
        v_in = v[n-1];
        w_in = w[n-1];
        T_in = T[n-1];
        tke_in = tke[n-1];
        omega_in = omega[n-1];
        break
      end

      if y[n] > args.y then
      -- Operates on cells that don't fall under the previous
      -- two categories. This operation should be working on most
      -- of the cells.
        w1 = (args.y - y[n-1])/(y[n] - y[n-1]) -- w1 => weighting function 1 ?
        w2 = (y[n] - args.y)/(y[n] - y[n-1])   -- w2 => weighting function 2 ?
        p_in = w2*p[n-1] + w1*p[n];
        u_in = w2*u[n-1] + w1*u[n];
        v_in = w2*v[n-1] + w1*v[n];
        w_in = w2*w[n-1] + w1*w[n];
        T_in = w2*T[n-1] + w1*T[n];
        tke_in = w2*tke[n-1] + w1*tke[n];
        omega_in = w2*omega[n-1] + w1*omega[n];
        break
      end
      n = n + 1
   end
      
   ghost = {}
   ghost.p = p_in    -- pressure, Pa
   ghost.T = {}      -- temperatures, K (as a table)
   ghost.T[0] = T_in
   ghost.u = u_in    -- x-velocity, m/s
   ghost.v = v_in    -- y-velocity, m/s
   ghost.w = w_in    -- z-velocity, m/s
   ghost.tke = tke_in
   ghost.omega = omega_in
   ghost.massf = {}  -- mass fractions to be provided as a table
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
   wall.u = u_in
   wall.v = v_in
   wall.w = 0.0
   wall.T = {[0]=316.2,}
   wall.massf = {[0]=1.0,}
   return wall
end
