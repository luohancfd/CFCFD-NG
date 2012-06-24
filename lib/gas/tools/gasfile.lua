#!/usr/bin/env lua
-- Author: Rowan J. Gollan
-- Date: 28-July-2008
--
-- History:
--   30-Mar-2009 - added thermally_perfect_gas option
--   20-May-2009 - added Noble-Abel gas option
--   04-Aug-2009 - added van der Waals gas option (Brendan T. O'Flaherty)
--   30-Aug-2009 - added two temperature gas option (Daniel F. Potter)

require 'ideal_gas'
require 'thermally_perfect_gas'
require 'Noble_Abel_gas'
require 'van_der_Waals_gas'
require 'four_temperature_gas'
require 'three_temperature_gas'
require 'two_temperature_gas'
require 'one_temperature_gas'
require 'fully_coupled_one_temperature_gas'
require 'fully_coupled_two_temperature_gas'
require 'fully_coupled_four_temperature_gas'

local create_ideal_gas = ideal_gas.create_ideal_gas
local create_thermally_perfect_gas = thermally_perfect_gas.create_thermally_perfect_gas
local create_Noble_Abel_gas = Noble_Abel_gas.create_Noble_Abel_gas
local create_van_der_Waals_gas = van_der_Waals_gas.create_van_der_Waals_gas
local create_four_temperature_gas = four_temperature_gas.create_three_temperature_gas
local create_three_temperature_gas = three_temperature_gas.create_three_temperature_gas
local create_two_temperature_gas = two_temperature_gas.create_two_temperature_gas
local create_one_temperature_gas = one_temperature_gas.create_one_temperature_gas
local create_fully_coupled_one_temperature_gas = fully_coupled_one_temperature_gas.create_fully_coupled_one_temperature_gas
local create_fully_coupled_two_temperature_gas = fully_coupled_two_temperature_gas.create_fully_coupled_two_temperature_gas
local create_fully_coupled_four_temperature_gas = fully_coupled_four_temperature_gas.create_fully_coupled_four_temperature_gas

gas_models = {}
gas_models["ideal gas"] = create_ideal_gas
gas_models["thermally perfect gas"] = create_thermally_perfect_gas
gas_models["Noble-Abel gas"] = create_Noble_Abel_gas
gas_models["van der Waals gas"] = create_van_der_Waals_gas
gas_models["four temperature gas"] = create_four_temperature_gas
gas_models["three temperature gas"] = create_three_temperature_gas
gas_models["two temperature gas"] = create_two_temperature_gas
gas_models["one temperature gas"] = create_one_temperature_gas
gas_models["fully coupled one temperature gas"] = create_fully_coupled_one_temperature_gas
gas_models["fully coupled two temperature gas"] = create_fully_coupled_two_temperature_gas
gas_models["fully coupled four temperature gas"] = create_fully_coupled_four_temperature_gas

function print_usage()
   print "Usage gasfile:"
   print "> gasfile input output"
   print ""
   os.exit(1)
end

function main()
   if #arg ~= 2 then
      print_usage()
   end

   input = arg[1]
   output = arg[2]

   -- 1. Load input file, see what we can do with it.
   dofile(input)

   -- 2. Now try to create an appropriately detailed
   --    output file for use by the C++ code.
   f = assert(io.open(output, "w"))
   print(string.format("%s model selected", model))
   if gas_models[model] ~= nil then
      gas_models[model](species, f)
   else
      print(string.format("Model name: %s is not known or not yet implemented.", model))
      print "Available models are:"
      for k,v in pairs(gas_models) do
	 print(k)
      end
      print("No output file created: ", output)
      os.exit(1)
   end
   f:close()
end

main()
