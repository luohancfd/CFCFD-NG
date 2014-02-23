-- udf-process.lua
-- This file sets up functions that will be called
-- from the main time-stepping loop.

print("Hello from the set-up stage of udf-process.")
print("nblks=", nblks)

function at_timestep_start(args)
   if (args.step ~= 0) then
      -- do nothing, just leave
      return
   end
   -- For the 0th step only
   mass = 0.0
   for ib=0,(nblks-1) do
      imin = blks[ib].imin; imax = blks[ib].imax
      jmin = blks[ib].jmin; jmax = blks[ib].jmax
      blk_id = blks[ib].id
      for j=jmin,jmax do
         for i=imin,imax do
            cell = sample_flow(blk_id, i, j, k)
            -- We are only given p and T
            -- so need to compute density
            -- using gas model
            Q = create_empty_gas_table()
            Q.p = cell.p
            Q.T = cell.T
            for isp=0,(nsp-1) do Q.massf[isp] = cell.massf[isp] end
            eval_thermo_state_pT(Q)
            rho = Q.rho
            -- Now we can compute mass in cell using volume of cell
            mass = mass + rho*cell.vol
         end
      end
   end
   print("Mass (kg) of gas in domain: ", mass)
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