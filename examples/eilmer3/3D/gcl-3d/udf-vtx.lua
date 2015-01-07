-- Lua script for the vertex velocity
--
-- Author: Jason Kan Qin
-- Date: 13-Dec-2014

-- dummy functions to keep eilmer3 happy

function at_timestep_start(args) return nil end
function at_timestep_end(args) return nil end

local imin = 2
local jmin = 2
local kmin = 2
local imax = 7
local jmax = 7
local kmax = 7

function vtx_velocity(args)

   src = {}
   i = args.i
   j = args.j
   k = args.k
   
   if i == imin or i == imax or j == jmin or j == jmax or k == kmin or k == kmax  then      
       src.vel_x = 0.0
       src.vel_y = 0.0
       src.vel_z = 0.0
   else
       src.vel_x = math.random(-1,1)/50.0
       src.vel_y = math.random(-1,1)/50.0
       src.vel_z = math.random(-1,1)/50.0
       --src.vel_x = 0.02
       --src.vel_y = 0.04
       --src.vel_z = 0.02
   end
   
   return src
   
end
