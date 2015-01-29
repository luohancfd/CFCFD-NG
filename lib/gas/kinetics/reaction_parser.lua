-- Author: Rowan J. Gollan
-- Date: 13-Mar-2009 (Friday the 13th)
-- Place: NIA, Hampton, Virginia, USA
--
-- History:
--   20-Mar-2009 :- first put into production
--

require 'reac'

local validate_reaction = reac.validate_reaction
local parse_reaction_string = reac.parse_reaction_string
local transform_reaction = reac.transform_reaction

reactions = {}
scheme_t = {}
ode_t = {}

--%----------------------------------
--  Functions available to the user
--  for configuration purposes.
--%----------------------------------

function scheme(t)
   for k,v in pairs(t) do
      scheme_t[k] = v
   end
end

function ode_solver(t)
   for k,v in pairs(t) do
      ode_t[k] = v
   end
end

function reaction(t)
   assert(validate_reaction(t))
   reactions[#reactions+1] = t
end

function remove_reactions_by_label(label)
   -- Build a list with all reactions
   -- which match label
   indices = {}
   for i,r in ipairs(reactions) do
      if label == 'string' then
	 if r.label == label then
	    indices[#indices+1] = i
	 end
      else
	 for _,l in ipairs(label) do
	    if r.label == l then
	       indices[#indices+1] = i
	    end
	 end
      end
   end

   table.sort(indices, function(a,b) return a > b end)

   for _,i in ipairs(indices) do
      table.remove(reactions, i)
   end

end

function remove_reactions_by_number(n)
   if type(n) == 'number' then
      table.remove(reactions, n)
   else
      table.sort(n, function(a,b) return a > b end)
      for _,i in ipairs(n) do
	 table.remove(reactions, i)
      end
   end
end

function select_reactions_by_label(t)
   
   -- 1. Make a list of ALL available labels
   excluded_labels = {}
   for _,r in ipairs(reactions) do
      excluded_labels[#excluded_labels+1] = r.label or 'no label'
   end

   -- 2. Build a list of those labels we are keeping
   --    that is, selecting.
   --    We will remove the 'selected' ones from the
   --    the list of ALL labels
   indices_to_remove = {}

   if type(t) == 'string' then
      for i,e in ipairs(excluded_labels) do
	 if e == t then
	    indices_to_remove[#indices_to_remove+1] = i
	 end
      end
   else
      for _,l in ipairs(t) do
	 for i,e in ipairs(excluded_labels) do
	    if e == l then
	       indices_to_remove[#indices_to_remove+1] = i
	    end
	 end
      end
   end

   -- 3. Now remove the selected ones from the list
   --    What remains is what we are removing
   table.sort(indices_to_remove, function (a,b) return a > b end)

   for _,i in ipairs(indices_to_remove) do
      table.remove(excluded_labels, i)
   end

   -- 4. Now pass off the non-selected ones for removal
   remove_reactions_by_label(excluded_labels)

end

function select_reactions_by_number(n)
   -- If n is a single integer, make into a table
   -- with one entry
   if type(n) == 'number' then
      n = {n}
   end

   excluded_reactions = {}
   for i=1,#reactions do
      excluded_reactions[#excluded_reactions+1] = i
   end

   table.sort(n, function (a,b) return a > b end)

   for _,i in ipairs(n) do
      table.remove(excluded_reactions, i)
   end

   -- Now pass of all those not selected to the
   -- removal function
   remove_reactions_by_number(excluded_reactions)
      
end

--%---------------------------------------------
--  Validating functions
--%---------------------------------------------

local function check_scheme()
   scheme_t.update = scheme_t.update or 'chemical kinetic ODE'
   scheme_t.error_tolerance = scheme_t.error_tolerance or 0.1
   scheme_t.temperature_limits = scheme_t.temperature_limits or 
      {lower=300, upper=50000}
   scheme_t.temperature_limits.lower = scheme_t.temperature_limits.lower or 20
   scheme_t.temperature_limits.upper = scheme_t.temperature_limits.upper or 100000
end

local function check_ode()
   ode_t.step_routine = ode_t.step_routine or 'qss'
   ode_t.max_step_attempts = ode_t.max_step_attempts or 4
   ode_t.max_increase_factor = ode_t.max_increase_factor or 1.15
   ode_t.max_decrease_factor = ode_t.max_decrease_factor or 0.01
   ode_t.decrease_factor = ode_t.decrease_factor or 0.333
end

function main(config_file)
   dofile(config_file)

   check_scheme()
   check_ode()

   for i,r in ipairs(reactions) do
      r.number = i
      reactions[i] = transform_reaction(r, species, SUPPRESS_WARNINGS)
   end

end



