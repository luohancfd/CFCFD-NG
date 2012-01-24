-- udf-mixing-plane-upstream-bc.lua
-- Lua script for the user-defined mixing plane testcase, upstream side.
-- called by the UserDefinedGhostCell BC.
-- M.P van der Laan, 09-Oct-2009

-- Setup exit pressure BC
--pressure_exit = "gradient"  			-- Use forced exit pressure
pressure_exit = "constant"			-- Use constant exit pressure
p_exit_set = 324.64e3 				-- Exit pressure

-- constants and definitions as stated in lua table:
R0 = 8.31451
M = 0.02897
R = R0/M                                        -- Gas constant for Argon [J/(kg.K)]
g = 1.4                                    	-- Ratio of specific heats Argon [-]

-- We want keep track of the TimeStep value.
TimeStep = 0
count = 0

function get_average_exit_pressure()
   -- Function for calculating density averaged static pressure in the downstream mixing plane cross section.
   i_mp = nni+1
   rho_sum = 0
   p_sum = 0
   p_exit_rho_av = {}
   for k = 2, nnk+1 do
      for j = 2, nnj+1 do 
   	 cell_flow = sample_flow(block_id, i_mp, j, k)
         
         rho = cell_flow.rho 
         p = cell_flow.p
         
         p_sum = p_sum + p*rho	 		-- Sum all p*rho
	 rho_sum =  rho_sum + rho 		-- Total density of flow cells downstream mixing plane cross section, for a constant radius
      end
   end
   p_exit_rho_av = p_sum/rho_sum 	        -- Density averaged pressure for a constant radius [Pa]
   print ("average exit pressure = ", p_exit_rho_av )  -- debug info

   return p_exit_rho_av
end

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   -- Set constant conditions across the whole boundary.
   -- We assume that the inflow plane is a (y,z)-plane and that
   -- the velocity vector lies in the cylindrical surface normal
   -- to the inflow plane.
   -- When this function is called we have the following available
   -- global variables: args{x,y,z,i,j,k}, block_id, nni, nnj, nnk,
   -- which_boundary.

   if pressure_exit == "gradient" then				-- Use forced exit pressure
      switch0 = 0; switch1 = 1  						
   -- Average static exit pressure from outlet side downstream block, calculated once for each timestep  
      count = count + 1
      if count > nnj*nnk*2*TimeStep then
         TimeStep= TimeStep + 1
         p_exit_average = get_average_exit_pressure() 
         delta_p = p_exit_set - p_exit_average
      end
   elseif pressure_exit == "constant" then			-- Use constant exit pressure
      switch0 = 1; switch1 = 0; delta_p = 0					
   end

   cell_flow = sample_flow(block_id, args.i, args.j, args.k)  	-- Adjacent cell properties

   ghost = {}
   ghost.p = switch0*p_exit_set + switch1*(cell_flow.p + delta_p)
   ghost.T = {}
   ghost.T[0] = ghost.p/cell_flow.rho/R         		-- Temperature [K]
   ghost.u = cell_flow.u             				-- x-velocity [m/s]
   ghost.v = cell_flow.v   					-- y-velocity [m/s]
   ghost.w = cell_flow.w   					-- z-velocity [m/s]
   ghost.massf = {}                          			-- Mass fractions to be provided as a table
   ghost.massf[0] = 1.0                      			-- Mass fractions are indexed from 0 to nsp-1
   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
