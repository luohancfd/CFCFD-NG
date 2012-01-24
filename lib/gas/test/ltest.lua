#!/usr/bin/env lua

--
-- global data tables for module
--

----------------------------------------------------
-- Table to declare the allowable types for the units
-- being tested.
--
local allowable_types = {
   'function',
   'method'
}

----------------------------------------------------
-- Table to declare the allowable results for the units
-- being tested.
--
local allowable_results = {
   'passed',
   'failed'
}

----------------------------------------------------
-- Checks if a value exists in a table.
-- One of the neat aspects of Lua is the ability
-- to extend the standard library easily.
--
-- @param t table to be interrogated
-- @param val value to test for
-- @return Boolean indicating if the value is found.
--
-- @author Rowan J. Gollan
-- @version 29-May-2008
--
function table.has_value(t, val)
   for _,v in pairs(t) do
      if val == v then
	 return true
      end
   end
   return false
end

--
-- functions to read formatted file
--

---------------------------------------------------
-- Gathers up the title string from a "formatted"
-- test file.
--
-- @param t_str title string for a set of unit tests
-- @return side-effect: set global value "test_title"
-- 
-- @author Rowan J. Gollan
-- @version 29-May-2008
function Title(t_str)
   test_title = t_str
end


---------------------------------------------------
-- Gathers up the file name from a "formatted" test
-- output file.  Note this name should be the name of the
-- file which contains the implementation NOT the
-- name of the testing program.
--
-- @param f_str file name string
-- @return side-effect: set global value "test_file"
--
-- @author Rowan J. Gollan
-- @version 29-May-2008
function File(f_str)
   test_file = f_str
end

units = {}

---------------------------------------------------
-- Begins gathering data associated with a unit test.
-- 
-- @param test_name name of the unit test
-- @return side-effect: begins a new entry in the
--         global units table
--
-- @author Rowan J. Gollan
-- @version 29-May-2008
function Test(test_name)
   units[#units + 1] = { name=test_name }
end

---------------------------------------------------
-- Gathers the type (ie. function or method) tested
-- and adds it to the currently active unit.
--
-- Note: it was important to name this function 'Type'
-- as use of 'type' (all lowercase) would have clashed
-- with the built-in lua function 'type'.
--
-- @param type_str should be either 'function' or 'method'
-- @return side-effect: add type to active unit table.
--
-- @author Rowan J. Gollan
-- @version 29-May-2008
function Type(type_str)
   if table.has_value(allowable_types, type_str) then
      units[#units].type = type_str
   else
      print(string.format("The type: %s is unknown", type_str))
      print("The allowable types are: ")
      for _,v in pairs(allowable_types) do
	 print(v)
      end
   end
end

---------------------------------------------------
-- Gathers the class name to which a method belongs
-- if appropriate.
--
-- @param c_str the name of the class
-- @return side-effect: add class name to active unit table
--
-- @author Rowan J. Gollan
-- @version 29-May-2008
--
function Class(c_str)
   units[#units].class = c_str
end

---------------------------------------------------
-- Gathers the result of the test.
--
-- @param r_str string indicating the result
-- @return side-effect: add result to active unit table
--
-- @author Rowan J. Gollan
-- @version 29-May-2008
--
function Result(r_str)
   if table.has_value(allowable_results, r_str) then
      units[#units].result = r_str
   else
      print(string.format("The result: %s is unknown", type_str))
      print("The allowable results are: ")
      for _,v in pairs(allowable_results) do
	 print(v)
      end
   end
end

function main()
   if #arg ~= 1 then
      print("Usage: ltest input.ltest")
      os.exit(1)
   end

   dofile(arg[1])

   print("==================================================")
   print(test_title)
   print("--------------------------------------------------")
   print("Testing units in file:", test_file)

   ntests = #units
   passed = {}
   failed = {}
   for i,u in ipairs(units) do
      print(string.format("unit %d/%d:", i, ntests))
      if u.type == 'function' then
	 print(string.format("    function:  %s", u.name))
      else
	 print(string.format("    class:     %s", u.class))
	 print(string.format("    method:    %s", u.name))
      end
      print(string.format("    test:      %s", u.result))

      if u.result == "passed" then
	 passed[#passed + 1] = i
      else
	 failed[#failed + 1] = i
      end
   end

   print("--------------------------------------------------")
   print("Summary:")
   print(string.format("  no. of tests:   %d", ntests))
   print(string.format("  no. passed:     %d", #passed))
   print(string.format("  no. failed:     %d", #failed))
   if #failed >= 1 then
      print("  failed tests: ")
      for _,idx in ipairs(failed) do
	 print(string.format("    %d. %s", idx, units[idx].name))
      end
   end
   print("--------------------------------------------------")

end

main()
