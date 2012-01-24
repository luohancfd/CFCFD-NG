-- This Lua script mimics the tests in:
-- geom_test.cxx and geom_test.py
--
-- The only difference is that some functionality
-- is not exposed to the Lua environment.  Some 
-- of this functionality is not supported by Lua
-- (eg. a += b: no '+=' in Lua) and some of it
-- has not needed to be exposed (as yet).
--
-- Author: R. J. Gollan
-- Date: 14-Mar-2008
-- 
-- based on work by P. A. Jacobs in:
-- geom_test.cxx
-- geom_test.py
--

-- Test if geometry objects functions are already available
if _G.Vector3 == nil then
   require "geometry"
   for k, v in pairs(geometry) do
      _G[k] = v
   end
end

local Nodes = {}
function Node(...)
   Nodes[#Nodes + 1] = Vector3(...)
   return Vector3(...)
end

function print_nodes()
   print("There are", #Nodes, "entries in Nodes table.")
   for i, node in ipairs(Nodes) do
      print(node)
   end
end

function main()
   print("Begin geom_test (in Lua)...")

   print("First, exercise the Vector3 class.")
   a = Vector3(1.1, 2.0, 3.0)
   print("initial : a=", a)
   -- In Lua, b = a would only make a reference to a
   b = a:copy()
   b = b + a
   print("b=a; b = b + a : b=", b)
   print("vabs(b)=", vabs(b))
   print("+a=", a, "; -b=", -b)
   c = a - b
   print("c=a-b: c=", c)
   c = 2.0*a
   print("c=2.0*a: c=", c)
   d = 2.0*a + Vector3(3.0)
   print("d=2.0*a+Vector3(3.0): d=", d)
   print("a == b : ", (a == b))
   print("c == d : ", (c == d))

   a = Vector3(math.sqrt(2.0), math.sqrt(2.0), 0.0)
   b = Vector3(-a.y, a.x, 0.0)
   print("a=", a, " b=", b, " cross(a,b)=", cross(a,b))
   print("unit(cross(a,b))=", unit(cross(a,b)))

end

main()
