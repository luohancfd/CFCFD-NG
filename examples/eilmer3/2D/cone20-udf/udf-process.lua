-- udf-process.lua
-- This file sets up functions that will be called
-- from the main time-stepping loop.

print("Hello from the set-up stage of udf-process.")
print("nblock=", nblock)

function at_timestep_start(args)
   -- Here, we could do things like run an external program
   -- and set up data files for use in the boundary conditions.
   --
   -- args contains {t, step}
   if (args.step % 100) == 0 then
      print("At start of timestep ", args.step, " t=", args.t)
   end
   return
end

function at_timestep_end(args)
   if (args.step % 100) == 0 then
      print("At end of timestep ", args.step, " t=", args.t)
   end
   return
end

function source_vector(args, cell_data)
   -- args contains t
   -- cell_data table contains most else
   src = {}
   src.mass = 0.0
   src.momemtum_x = 0.0
   src.momentum_y = 0.0
   src.momentum_z = 0.0
   if cell_data.x > 0.2 and cell_data.x < 0.4 and 
      cell_data.y > 0.2 and cell_data.y < 0.3 then
      src.total_energy = 100.0e+6  -- J/m**3
      src.energies = {[0]=100.0e6}
   else
      src.total_energy = 0.0
      src.energies = {[0]=0.0}
   end
   src.species = {[0]=0.0,}
   return src
end