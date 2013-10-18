-- udf-process.lua
-- This file sets up functions that will be called
-- from the main time-stepping loop.

print("Hello from the set-up stage of udf-process.")
print("nblock=", nblock)

-- constants required for the porous block
K = 1.0e-6		-- permeability [1/m]
phi = 0.8		-- porosity [non-dimensional]
C_f = 1.5               -- inertial coefficient/Ergun constant

-- flow conditions of the fluid within the porous block

function at_timestep_start(args)
   -- Here, we could do things like run an external program
   -- and set up data files for use in the boundary conditions.
   --
   -- args contains {t, step}
   if (args.step % 1000) == 0 then
      print("At start of timestep ", args.step, " t=", args.t)
   end
   return
end

function at_timestep_end(args)
   if (args.step % 1000) == 0 then
      print("At end of timestep ", args.step, " t=", args.t)
   end
   return
end

function source_vector(args, cell)
   -- args contains t
   -- cell table contains most else

   Q = create_empty_gas_table()
   Q.p = cell.p
   Q.T[0] = cell.T[0]
   for i=0,nsp-1 do
      Q.massf[i] = cell.massf[i]
   end
   eval_thermo_state_pT(Q)
   eval_transport_coefficients(Q)
   mu = Q.mu 
   rho = cell.rho
   u = 0.0
   v = cell.v      -- or set the same as the high p inflow, if not working
   src = {}
   src.mass = 0.0
   src.momentum_z = 0.0
   if cell.x > 0.3 and cell.x < 0.6 and  -- change these dimensions if you change the geometry!
      cell.y > 0.08 and cell.y < 0.1 then
         src.momentum_x = -mu*u/(K*phi) - C_f*phi*rho*u*u/(K^0.5)
         src.momentum_y = -mu*v/(K*phi) - C_f*phi*rho*v*v/(K^0.5)   
      --src.momentum_x = -10000.0    -- dummy values for testing
      --src.momentum_y = -10000.0
   else
   src.momentum_x = 0.0
   src.momentum_y = 0.0
   end
   src.total_energy = 0.0
   src.energies = {[0]=0.0}
   src.species = {[0]=0.0,}
   return src
end

