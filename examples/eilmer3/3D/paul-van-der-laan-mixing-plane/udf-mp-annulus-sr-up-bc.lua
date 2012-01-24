-- udf-mixing-plane-upstream-bc.lua
-- Lua script for the user-defined mixing plane testcase, upstream side.
-- called by the UserDefinedGhostCell BC.
-- M.P van der Laan, 09-Oct-2009

-- constants and definitions as stated in lua table:
R0 = 8.31451
M = 0.039948
R = R0/M                                   -- Gas constant for Argon [J/(kg.K)]
g = 5/3                                    -- Ratio of specific heats Argon [-]

-- We want keep track of the number of TimeSteps.
TimeStep = 0
count = 0

function get_average_pressure()
   -- Function for calculating density averaged static pressure in the downstream mixing plane cross section.
   -- Averaging is done is for constant radii.
   i_mp = 2
   p_radial_rho_av = {}
   for k = 2, nnk+1 do
      rho_sum = 0
      p_sum = 0
      for j = 2, nnj+1 do 
   	 cell_flow = sample_flow(block_id+1+3, i_mp, j, k)
         
         rho = cell_flow.rho 
         p = cell_flow.p
         
         p_sum = p_sum + p*rho	 			-- Sum all p*rho
	 rho_sum =  rho_sum + rho			-- Total density of flow cells downstream mixing plane cross section, for a constant radius
      end
      p_radial_rho_av[k] = p_sum/rho_sum 	        -- Density averaged pressure for a constant radius [Pa]
   end
   return p_radial_rho_av
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

   -- Average static pressure from downstream block, calculated once for each timestep  
   count = count + 1
   k = args.k 	--radial_index
   if count > nnj*nnk*2*TimeStep then
      TimeStep= TimeStep +1
      pghost = get_average_pressure() 
   end

   cell_flow = sample_flow(block_id, args.i, args.j, args.k) -- Adjacent cell properties
   
   ghost = {}
   ghost.p = pghost[k]				-- Pressure [Pa]
   ghost.T = {}
   ghost.T[0] = pghost[k]/cell_flow.rho/R       -- Temperature [K]
   ghost.u = cell_flow.u             		-- x-velocity [m/s]
   ghost.v = cell_flow.v   			-- y-velocity [m/s]
   ghost.w = cell_flow.w   			-- z-velocity [m/s]
   ghost.massf = {}                          	-- Mass fractions to be provided as a table
   ghost.massf[0] = 1.0                      	-- Mass fractions are indexed from 0 to nsp-1
   return ghost, ghost
end


function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   return sample_flow(block_id, args.i, args.j, args.k)
end
