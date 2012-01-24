-- udf-supersonic-in.lua
-- Lua script for the user-defined functions 
-- called by the UserDefinedGhostCell BC.
--
-- This script is modified from the cone20-udf case 
-- and is used because we want to specify a supersonic
-- inflow condition that is variable (by a function) 
-- throughout the entire inflow boundary.
-- Wilson Chan 05-Nov-2008


function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- but we don't happen to use any of them.
   
   -- Some brief calculations to obtain inflow conditions at different
   -- axial locations.
   -- Start at cylinder leading edge at Mach 8.8 with an increasing
   -- axial mach no. gradient of 0.24/m (as specified in papers by
   -- Mallinson et al. and Boyce & Hillier.
   mach_inf = 8.8 + 0.24 * args.x
   -- Isentropic calculation, T_inf = T_0 / (1 + (g-1)/2*M^2)
   T_inf = 1150.0 / (1.0 + (1.4-1.0)/2.0 * mach_inf^2)
   -- Isentropic calculation, p_inf = p_0 / (T_0/T)^(g/(g-1))
   p_inf = 60.0e6 / (1150.0/T_inf)^3.5
   u_nominal = mach_inf * (1.4 * 297 * T_inf)^0.5

   -- Estimate turbulence quantities for free stream by specifying turbulence
   -- intensity and the ratio of turbulent-to-laminar viscosity
   --turbulence_intensity = 0.005
   --turb_lam_viscosity_ratio = 0.5
   --tke_inf = 1.5 * (turbulence_intensity * u_nominal)^2
   -- Sutherland's viscosity law for N2
   -- mu_inf =  mu_ref * (T / T_ref)**1.5 * (T_ref + S)/(T + S)
   ---mu_inf =  1.663e-5 * (T_inf / 273.0)^1.5 * (273.0 + 107.0)/(T_inf + 107.0)
   --mu_t_inf = turb_lam_viscosity_ratio * mu_inf
   --omega_inf = p_inf / (297.0 * T_inf) * tke_inf / mu_t_inf  
   tke_inf = 1.0e-12
   omega_inf = 1.0
   -- where p_inf/(R*T_inf) = rho_inf (Perfect gas equation)
   
   -- Input inflow conditions for each cell along boundary
   ghost = {}
   ghost.p = p_inf -- pressure, Pa
   ghost.T = {}  -- temperatures, K (as a table)
   ghost.T[0] = T_inf
   ghost.u = u_nominal * math.cos(math.rad(5.0 * args.y))
   ghost.v = u_nominal * math.sin(math.rad(5.0 * args.y))
   ghost.w = 0.0
   ghost.tke = tke_inf
   ghost.omega = omega_inf
   ghost.massf = {} -- mass fractions to be provided as a table
   ghost.massf[0] = 1.0 -- mass fractions are indexed from 0 to nsp-1
   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   -- but we don't happen to use any of them.
   
   -- Some brief calculations to obtain inflow conditions at different
   -- axial locations.
   -- Start at cylinder leading edge at Mach 8.8 with an increasing
   -- axial mach no. gradient of 0.24/m
   mach_inf = 8.8 + 0.24 * args.x
   -- Isentropic calculation, T_inf = T_0 / (1 + (g-1)/2*M^2)
   T_inf = 1150.0 / (1.0 + (1.4-1.0)/2.0 * mach_inf^2)
   u_nominal = mach_inf * (1.4 * 297 * T_inf)^0.5

   wall = {}
   wall.u = u_nominal * math.cos(math.rad(5.0 * args.y))
   wall.v = u_nominal * math.sin(math.rad(5.0 * args.y))
   wall.w = 0.0
   wall.T_wall = 295.0 
   wall.massf = {}
   wall.massf[0] = 1.0
   return wall
end
