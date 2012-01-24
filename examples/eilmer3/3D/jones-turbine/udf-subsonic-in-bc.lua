-- udf-subsonic-bc.lua
-- Lua script for the user-defined subsonic inflow at Radial Stator Inlet
-- called by the UserDefinedGhostCell BC.
-- PJ, 03-Nov-2008
-- Modified by M.P. van der Laan Nov-2009 for MP-model
-- Modified by PJ Petrie-Repar Aprv-2010 

-- input parameters:
T0 = 1056.48                                 -- Total temp in [K]
p0 = 513.15e3                               -- Total pressure [Pa]
--alpha = math.rad(25)                        -- Inflow angle [rad]
beta = math.rad(0) 			   -- Inflow beta [rad]	
alpha = math.rad(76)

-- constants and definitions:
R0 = 8.31451
M = 0.02897
R = R0/M                                   -- Gas constant for Argon [J/(kg.K)]
g = 1.4                                    -- Ratio of specific heats for Argon
Cp = g*R/(g-1)                             -- Specific-heat, constant volume [J/(kg.K)]
h0 = Cp*T0                                 -- Total enthalpy [J/kg]

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   -- Set constant conditions across the whole boundary.
   -- We assume that the inflow plane is a (y,z)-plane and that
   -- the velocity vector lies in the cylindrical surface normal
   -- to the inflow plane.
   x = args.x; y = args.y; z = args.z
   r = math.sqrt(x*x + y*y)   					-- r is radial coordinate in x,y-plane


   cell_flow = sample_flow(block_id, args.i, args.j, args.k) 	-- Adjacent cell properties

   vel_sq = cell_flow.u^2+cell_flow.v^2+cell_flow.w^2        	-- Square of inflow velocity [m^2/s^2]
   vel = math.sqrt(vel_sq)                   			-- Inflow velocity [m/s]
   M_sq = vel_sq/(h0-0.5*vel_sq)/(g-1)       			-- Square of Mach number [-]
   ratio = 1+0.5*(g-1)*M_sq                  			-- T0/T [-]
	
   vel_r = - vel * math.cos(alpha)
   vel_t = vel * math.sin(alpha)

--   print ("radius = ", r, "vel =", vel, "alpha =", alpha)  -- debug info
--   print ("vel_r = ", vel_r, "vel_t =", vel_t, "x y =", x, y)  -- debug info

   ghost = {}
   ghost.p = p0/math.pow(ratio,(g/(g-1)))    			-- Pressure [Pa]
   ghost.T = {}
   ghost.T[0] = (h0-0.5*vel_sq)/Cp           			-- Temperature [K]
   ghost.u = (-y * vel_t + x * vel_r) / r    -- x-velocity [m/s]
   ghost.v = (x * vel_t + y * vel_r) / r     -- y-velocity [m/s]
   ghost.w = 0.0			             -- z-velocity [m/s]
   ghost.massf = {}                          			-- Mass fractions to be provided as a table
   ghost.massf[0] = 1.0                      			-- Mass fractions are indexed from 0 to nsp-1  
--   print ("pres = ", ghost.p, "temp =", ghost.T[0])   -- debug info
--   print ("u = ", ghost.u, "v =", ghost.v)            -- debug info

   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
