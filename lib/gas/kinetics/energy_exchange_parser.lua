-- Author: Rowan J. Gollan
-- Date: 30-Aug-2012
-- Place: Greenslopes, Queensland, Australia
--

require 'exch'

local parse_energy_exch_string = exch.parse_energy_exch_string
local validate_mechanism = exch.validate_mechanism
local transform_mechanism = exch.transform_mechanism

-- Set defaults for scheme and ode solver
scheme_t = {
   update = "energy exchange ODE",
   temperature_limits = {
      lower = 20.0,
      upper = 100000.0
   },
   error_tolerance = 0.000001
}
ode_t = {
   step_routine = 'rkf',
   max_step_attempts = 4,
   max_increase_factor = 1.15,
   max_decrease_factor = 0.01,
   decrease_factor = 0.333
}
-- Empty table to gather mechanisms
mechanisms = {}

---------------------------------------------
-- Functions available to the user for
-- configuration purposes.
---------------------------------------------

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

function mechanism(t)
   --assert(validate_mechanism(t))
   mechanisms[#mechanisms+1] = t
end

------------------------------------------
-- Main routine
------------------------------------------

function main(config_file)
   dofile(config_file)

   mechs = {}
   for i,m in ipairs(mechanisms) do
      local t = transform_mechanism(m, species, modes)
      for _,m2 in ipairs(t) do
	 mechs[#mechs+1] = m2
      end
   end
   
   -- Find the largest mode index
   imode_max = 0
   for _,m in ipairs(mechs) do
      if m.imode > imode_max then
         imode_max = m.imode
      end
   end
   
   -- Initialise the rates table
   rates = {}
   for ir=1, imode_max, 1 do
       rates[ir] = {}
       rates[ir].imode = ir
       rates[ir].mechanisms = {}
   end

   -- Gather up rates.
   for i,m in ipairs(mechs) do
      ir = m.imode
      rates[ir].mechanisms[#rates[ir].mechanisms+1] = m
   end
   
   --  Set imode and mechanisms for rates without mechanisms
   for i,r in ipairs(rates) do
       print(i)
       if not r then
          r = {}
          r.imode = i
          r.mechanisms = {}
       end
    end

end

