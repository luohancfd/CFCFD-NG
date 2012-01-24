-- udf-mixing-plane-downstream-bc.lua
-- Lua script for the user-defined mixing plane testcase, downstream side.
-- called by the UserDefinedGhostCell BC.
-- M.P van der Laan, 09-Oct-2009

-- constants and definitions:
R0 = 8.31451
M = 0.039948
R = R0/M                                   -- Gas constant for Argon [J/(kg.K)]
g = 5/3                                    -- ratio of specific heats for Argon
Cp = g*R/(g-1)                             -- specific-heat, constant volume [J/(kg.K)]
omega = 0.0--2000				   -- rotor rotational speed [rad/s]

-- We want keep track of the number of TimeSteps. 
TimeStep = 0
count = 0

function get_cylindrical_flow(u, v, x, y)
   -- Function that transforms cartesian velocities to cylindrical velocities in xy plane.
   r = math.sqrt(x*x + y*y)

   u_theta = -y/r*u + x/r*v
   u_r     =  x/r*u + y/r*v
   return u_theta, u_r
end

function get_cartesian_flow(u_theta, u_r, x, y)
   -- Function that calculates cartesian velocities in xy plane.
   r = math.sqrt(x*x + y*y)

   u = -y/r*u_theta + x/r*u_r 
   v =  x/r*u_theta + y/r*u_r 
   return u, v
end

function get_relative_flow(alpha, beta, C, x, y)
   -- Function that calculates relative cartesian velocities for the moving rotor.
   r = math.sqrt(x*x + y*y)

   U = omega*r				  -- Rotor speed [rad/s]

   u_theta = C*math.sin(alpha)*math.cos(beta) - U
   u_r = C*math.sin(beta)

   w = C*math.cos(alpha)*math.cos(beta)
   u, v = get_cartesian_flow(u_theta, u_r, x, y)
   return u, v, w
end

function get_flow_upstream()
   -- Function for calculating density averaged alpha, beta, T0 and p0 in upstream flow cells.
   -- Averaging is done is for constant radii.
   i = nni+1
   T0_radial_rho_av = {}
   p0_radial_rho_av = {}
   alpha_radial_rho_av = {}
   beta_radial_rho_av = {}
   for k = 2, nnk+1 do
      T0_sum = 0
      p0_sum = 0
      alpha_sum = 0
      beta_sum = 0
      rho_sum = 0
      for j = 2, nnj+1 do 
         cell_flow = sample_flow(block_id-1-3, i, j, k)        	-- Get flow in upstream block 

         vel_sq = cell_flow.u^2+cell_flow.v^2+cell_flow.w^2    	-- Square of inflow velocity [m^2/s^2]
   	 vel = math.sqrt(vel_sq)                   		-- Inflow velocity [m/s]
   	 M_sq = vel_sq/cell_flow.a^2      			-- Square of Mach number [-]
   	 ratio = 1+0.5*(g-1)*M_sq				-- T0/T [-]

         T0 = ratio*cell_flow.T[0]                  		-- Total Temperature [K]
         p0 = cell_flow.p*math.pow(ratio,(g/(g-1)))    		-- Total pressure [Pa]

         u_theta, u_r = get_cylindrical_flow(cell_flow.u, cell_flow.v, cell_flow.x, cell_flow.y)  --u_theta_rel --> u_theta
         alpha = math.atan(u_theta/cell_flow.w)             	-- Inflow angle beta
         beta = math.asin(u_r/vel)	                        -- Inflow angle alpha


         rho = cell_flow.rho 

      	 alpha_sum = alpha_sum + alpha*rho                      -- Sum all alpha*rho for a constant radius
      	 beta_sum = beta_sum + beta*rho 			-- Sum all beta*rho for a constant radius
      	 T0_sum = T0_sum + T0*rho 				-- Sum all T0*rho for a constant radius
      	 p0_sum = p0_sum + p0*rho				-- Sum all p0*rho for a constant radius
         
	 rho_sum =  rho_sum + rho         	                -- Total Density of flow cells downstream mixing plane cross section, for a constant radius
      end
   alpha_radial_rho_av[k] = alpha_sum/rho_sum	                -- Density averaged alpha for a constant radius [rad]
   beta_radial_rho_av[k] = beta_sum/rho_sum 			-- Density averaged beta  for a constant radius [rad]
   T0_radial_rho_av[k] = T0_sum/rho_sum 			-- Density averaged total temperature  for a constant radius [K]
   p0_radial_rho_av[k] = p0_sum/rho_sum  			-- Density averaged total pressure for a constant radius [Pa]
   end
   return alpha_radial_rho_av, beta_radial_rho_av, T0_radial_rho_av, p0_radial_rho_av
end

function ghost_cell(args)
   -- Function that returns the flow states for the ghost cells 
   -- at the downstream side of the Mixing Plane.
   -- For use in the inviscid flux calculations.
   -- Set constant conditions across the whole boundary.
   -- We assume that the inflow plane is a (y,z)-plane.
   -- When this function is called we have the following available
   -- global variables: args{x,y,z,i,j,k}, block_id, nni, nnj, nnk,
   -- which_boundary.

   r = math.sqrt(args.x^2 + args.y^2)   
   k = args.k --Radial index

   -- Averaged T0 and p0 in upstream flow cells, calculated once per timestep
   count = count +1
   if count > nnj*nnk*2*TimeStep then
      TimeStep= TimeStep +1
      alpha, beta, T0, p0 = get_flow_upstream()
   end

   -- Calculate flow states for ghost cells using absolute variables 
   cell_flow = sample_flow(block_id, args.i, args.j, args.k)	-- Adjacent cell properties
   u_theta_rel, u_r = get_cylindrical_flow(cell_flow.u, cell_flow.v, cell_flow.x, cell_flow.y)
   u_theta_abs = u_theta_rel + omega*r 				-- Absolute tangetial velocity

   vel_sq = u_theta_abs^2+u_r^2+cell_flow.w^2   		-- Square of inflow velocity [m^2/s^2]
   vel = math.sqrt(vel_sq)                   			-- Inflow velocity [m/s]
   h0 = Cp*T0[k]						-- Enthalpy
   M_sq = vel_sq/(h0-0.5*vel_sq)/(g-1)				-- Square of Mach number [-]
   ratio = 1+0.5*(g-1)*M_sq                  			-- T0/T [-]

   ghost = {}
   ghost.p = p0[k]/math.pow(ratio,(g/(g-1)))    		-- Pressure [Pa]
   ghost.T = {}
   ghost.T[0] = T0[k]/ratio  					-- Temperature [K]
   ghost.u, ghost.v, ghost.w = get_relative_flow(alpha[k], beta[k], vel, cell_flow.x, cell_flow.y)
   ghost.massf = {}                          			-- Mass fractions to be provided as a table
   ghost.massf[0] = 1.0                      			-- Mass fractions are indexed from 0 to nsp-1
   return ghost, ghost
end

function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
