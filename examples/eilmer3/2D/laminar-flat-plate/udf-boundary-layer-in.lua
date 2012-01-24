-- udf-boundary-layer-in.lua
-- This particular example defines the flow from an external data file.
-- PJ, 7-Mar-2011

-- read the y,p,T,u tables
dofile("profile.lua")

function find_nearest_profile_point(y_given)
   -- returns an index into the profile tables
   local d, min_d, nearest
   min_d = math.abs(y_given - y[0])
   nearest = 0
   for i=1,table.maxn(y) do
      d = math.abs(y_given - y[i])
      if d < min_d then
	 min_d = d; nearest = i
      end
   end	 
   return nearest
end

-- These tables will eventually contain indices into the profile tables
-- that indicate which profile point is closest to each cell and face. 
cell_index = {}
face_index = {}

function ghost_cell(args)
   -- Function that returns the flow states for a ghost cells.
   -- For use in the inviscid flux calculations.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   --
   -- Set conditions using data from the nearest profile point.
   local i
   if cell_index[args.j] == nil then
      cell_index[args.j] = find_nearest_profile_point(args.y)
   end
   i = cell_index[args.j]
   ghost = {}
   ghost.p = p[i]
   ghost.T = {}
   ghost.T[0] = T[i]
   ghost.u = u[i]
   ghost.v = 0.0
   ghost.w = 0.0
   ghost.massf = {}
   ghost.massf[0] = 1.0
   return ghost, ghost
end

function interface(args)
   -- Function that returns the conditions at the boundary 
   -- when viscous terms are active.
   --
   -- args contains t, x, y, z, csX, csY, csZ, i, j, k, which_boundary
   local i
   if face_index[args.j] == nil then
      face_index[args.j] = find_nearest_profile_point(args.y)
   end
   i = face_index[args.j]
   face = {}
   face.u = u[i]
   face.v = 0.0
   face.w = 0.0
   face.T_wall = T[i]
   face.massf = {}
   face.massf[0] = 1.0
   return face
end
