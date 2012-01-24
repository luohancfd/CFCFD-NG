-- Lua script for the source terms
-- of a Manufactured Solution which
-- treats Euler flow.
--
-- Author: Rowan J. Gollan
-- Date: 04-Feb-2008

-- dummy functions to keep eilmer3 happy

function at_timestep_start(args) return nil end
function at_timestep_end(args) return nil end

local T_i = 362.58
local alpha = 1000

-- Heaviside step function
local function H(T)
   if T >= T_i then
      return 1
   else
      return 0
   end
end

function source_vector(args, cell)
   src = {}
   src.mass = 0
   src.momentum_x = 0
   src.momentum_y = 0
   src.momentum_z = 0
   src.total_energy = 0
   src.species = {}
   src.species[0] = -alpha*cell.rho*cell.massf[0]*H(cell.T[0])
   src.species[1] = alpha*cell.rho*cell.massf[0]*H(cell.T[0])
   return src
end
