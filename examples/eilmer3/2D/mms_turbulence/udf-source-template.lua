-- udf-source-template.lua
-- Lua template for the source terms of a Manufactured Solution.
--
-- PJ, 29-May-2011
-- RJG, 06-Jun-2014
--   Declared maths functions as local
-- Jianyong Wang, 13-July-2016
--   To handle the min/max functions employed by Wilcox (2006) k-w
--   turbulence model and be compatible with file of "make_source_terms.py",
--   some functions are defined for generating source terms properly.

-- dummy functions to keep eilmer3 happy
function at_timestep_start(args) return nil end
function at_timestep_end(args) return nil end

function Heaviside(x)
   if x < 0.0 then 
      return 0.0
   else
      return 1.0
   end
end

function Min1(x, y)
   if y > x then
       return 0.125
   else
       return 0.0
   end
end

function Min2(x, y)
   if y > x then
       return x
   else
       return y
   end
end

local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi
local max = math.max
local sqrt = math.sqrt

function source_vector(args, cell)
   src = {}
   x = cell.x
   y = cell.y
<insert-source-terms-here>
   src.mass = fmass
   src.momentum_x = fxmom
   src.momentum_y = fymom
   src.momentum_z = 0.0
   src.total_energy = fe
   src.rtke = ftke
   src.romega = fomega
   src.species = {}
   src.species[0] = src.mass
   return src
end
