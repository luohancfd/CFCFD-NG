-- udf-source-template.lua
-- Lua template for the source terms of a Manufactured Solution.
--
-- PJ, 29-May-2011
--     07-Nov-2012 yzx version

-- dummy functions to keep eilmer3 happy
function at_timestep_start(args) return nil end
function at_timestep_end(args) return nil end

function source_vector(args, cell)
   src = {}
   y = cell.y
   z = cell.z
<insert-source-terms-here>
   src.mass = fmass
   src.momentum_z = fzmom
   src.momentum_y = fymom
   src.momentum_x = 0.0
   src.total_energy = fe
   src.species = {}
   src.species[0] = src.mass
   return src
end
