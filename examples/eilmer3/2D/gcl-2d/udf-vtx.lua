-- Lua script for the vertex velocity
--
-- Author: Jason Kan Qin
-- Date: 13-Dec-2014

-- dummy functions to keep eilmer3 happy

function at_timestep_start(args) return nil end
function at_timestep_end(args) return nil end

function before_grid_update(args)
    if (args.step == 1) then
        print("HELLO")
        print("block_id= ", args.block_id)
    end
    return nil
end

local imin = 2
local jmin = 2
local imax = 12
local jmax = 12

function vtx_velocity(args)

   src = {}
   i = args.i
   j = args.j
   k = args.k   
   
   if i == imin or i == imax or j == jmin or j == jmax then
       
       src.vel_x = 0.0
       src.vel_y = 0.0
       src.vel_z = 0.0
   else
       --src.vel_x = 0.1
       --src.vel_y = 0.1
       src.vel_x = math.random(-2,2)/10.0
       src.vel_y = math.random(-2,2)/10.0
       src.vel_z = 0.0
   end
   
   return src
   
end
