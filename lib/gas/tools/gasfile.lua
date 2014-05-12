#!/usr/bin/env lua
-- Author: Rowan J. Gollan
-- Date: 28-July-2008
--
-- History:
--   30-Mar-2009 - added thermally_perfect_gas option
--   20-May-2009 - added Noble-Abel gas option
--   04-Aug-2009 - added van der Waals gas option (Brendan T. O'Flaherty)
--   30-Aug-2009 - added two temperature gas option (Daniel F. Potter)
--   12-Mar-2012 - added real gas Bender option (Peter Blyton)
--   18-Apr-2012 - added real gas MBWR option (Peter Blyton)
--   15-May-2012 - added real gas REFPROP option (Peter Blyton)

require 'tab_serialise'
require 'refprop'
require 'four_temperature_gas'
require 'three_temperature_gas'
require 'two_temperature_gas'
require 'one_temperature_gas'
require 'fully_coupled_one_temperature_gas'
require 'fully_coupled_two_temperature_gas'
require 'fully_coupled_four_temperature_gas'
require 'multi_temperature_gas'

local serialise = tab_serialise.serialise

-- For each gas model, list which properties
-- are required from the species database.
ideal_gas = {
   value_list = {"M", "gamma", "d", "e_zero", "q"},
   model_list = {"viscosity", "thermal_conductivity"},
   EOS = "perfect gas",
   TBM = "constant specific heats",
   mixing_rule = "Wilke",
   sound_speed = "equilibrium",
   diffusion_coeffs = "hard sphere"
}

simple_gas = {
   value_list = {"M", "gamma", "d", "e_zero", "q"},
   model_list = {"viscosity", "thermal_conductivity"},
   EOS = "simple gas",
   TBM = "constant specific heats",
   mixing_rule = "Wilke",
   sound_speed = "equilibrium",
   diffusion_coeffs = "hard sphere"
}

thermally_perfect_gas = {
   value_list = {"M", "d"},
   model_list = {"viscosity", "thermal_conductivity"},
   EOS = "perfect gas",
   TBM = "thermally perfect",
   mixing_rule = "Wilke",
   sound_speed = "equilibrium",
   diffusion_coeffs = "hard sphere"
}

Noble_Abel_gas = {
   value_list = {"M", "gamma", "d", "T_c", "p_c"},
   model_list = {"viscosity", "thermal_conductivity"},
   EOS = "Noble-Abel gas",
   TBM = "thermally real",
   mixing_rule = "Wilke",
   sound_speed = "equilibrium",
   diffusion_coeffs = "hard sphere"
}

van_der_Waals_gas = {
   value_list = {"M", "gamma", "d", "T_c", "p_c"},
   model_list = {"viscosity", "thermal_conductivity"},
   EOS = "van der Waals gas",
   TBM = "thermally real",
   mixing_rule = "Wilke",
   sound_speed = "equilibrium",
   diffusion_coeffs = "hard sphere"
}

real_gas_Bender = {
   value_list = {"M", "reference_state", "Cp0_coeffs",
		 "Bender_EOS_coeffs"},
   model_list = {"viscosity", "thermal_conductivity"},
   EOS = "Bender",
   TBM = "dense thermally real",
   mixing_rule = "Wilke",
   sound_speed = "equilibrium",
   diffusion_coeffs = "hard sphere"
}

real_gas_MBWR = {
   value_list = {"M", "rho_c", "reference_state",
		 "MBWR_EOS_coeffs", "Cp0_coeffs"},
   model_list = {"viscosity", "thermal_conductivity"},
   EOS = "MBWR",
   TBM = "dense thermally real",
   mixing_rule = "Wilke",
   sound_speed = "equilibrium",
   diffusion_coeffs = "hard sphere"
}

real_gas_REFPROP = {
   special_fn = refprop.create_real_gas_REFPROP
}

four_T_gas = {
   special_fn = four_temperature_gas.create_four_temperature_gas
}

three_T_gas = {
   special_fn = three_temperature_gas.create_three_temperature_gas
}

two_T_gas = {
   special_fn = two_temperature_gas.create_two_temperature_gas
}

one_T_gas = {
   special_fn = one_temperature_gas.create_one_temperature_gas
}

fully_coupled_one_T_gas = {
   special_fn = fully_coupled_one_temperature_gas.create_fully_coupled_one_temperature_gas
}

fully_coupled_two_T_gas = {
   special_fn = fully_coupled_two_temperature_gas.create_fully_coupled_two_temperature_gas
}

fully_coupled_four_T_gas = {
   special_fn = fully_coupled_four_temperature_gas.create_fully_coupled_four_temperature_gas
}

multi_T_gas = {
   special_fn = multi_temperature_gas.create_multi_temperature_gas
}

gas_models = {}
gas_models["ideal gas"] = ideal_gas
gas_models["simple gas"] = simple_gas
gas_models["thermally perfect gas"] = thermally_perfect_gas
gas_models["Noble-Abel gas"] = Noble_Abel_gas
gas_models["van der Waals gas"] = van_der_Waals_gas
gas_models["real gas Bender"] = real_gas_Bender
gas_models["real gas MBWR"] = real_gas_MBWR
gas_models["real gas REFPROP"] = real_gas_REFPROP
gas_models["four temperature gas"] = four_T_gas
gas_models["three temperature gas"] = three_T_gas
gas_models["two temperature gas"] = two_T_gas
gas_models["one temperature gas"] = one_T_gas
gas_models["fully coupled one temperature gas"] = fully_coupled_one_T_gas
gas_models["fully coupled two temperature gas"] = fully_coupled_two_T_gas
gas_models["fully coupled four temperature gas"] = fully_coupled_four_T_gas
gas_models["multi temperature gas"] = multi_T_gas

function print_usage()
   print "Usage gasfile:"
   print "> gasfile input output"
   print ""
   os.exit(1)
end
function list_available_species()
   e3bin = os.getenv("E3BIN") or os.getenv("HOME").."/e3bin"
   dir = e3bin.."/species"
   tmpname = os.tmpname()
   os.execute(string.format("ls -1 %s/*.lua > %s", dir, tmpname))
   tmpfile = assert(io.open(tmpname, "r"))
   species = {}
   for line in tmpfile:lines() do
      sp = string.match(line, "([%a%d_]+).lua")
      species[sp] = true
   end
   tmpfile:close()
   os.execute(string.format("rm %s", tmpname))
   return species
end

function create_gas_file(species, t, f)
   species_avail = list_available_species()
   
   f:write(string.format("-- Auto-generated by gasfile on: %s\n",
			 os.date("%d-%b-%Y %X")))
   f:write("model = 'composite gas'\n")
   f:write(string.format("equation_of_state = '%s'\n", t.EOS))
   f:write(string.format("thermal_behaviour = '%s'\n", t.TBM))
   f:write(string.format("mixing_rule = '%s'\n", t.mixing_rule))
   f:write(string.format("sound_speed = '%s'\n", t.sound_speed))
   if #species > 1 then
      f:write(string.format("diffusion_coefficients = '%s'\n", t.diffusion_coeffs))
   end
   f:write("ignore_mole_fraction = 1.0e-15\n")
   if _G.T_COLD then
      f:write(string.format("T_COLD = %f\n", _G.T_COLD))
   else
      f:write("T_COLD = 20.0\n")
   end

   species2 = {} -- converted species name (if necessary)
   f:write("species = {")
   for _,sp in ipairs(species) do
      -- Convert species names like N2+ --> N2_plus
      if string.match(sp, '+') then
          sp = string.gsub(sp, '+', '_plus')
      end
      if sp == 'e-' then
          sp = 'e_minus'
      end
      if not species_avail[sp] then
	 print(string.format("Species: %s cannot be found in the collection of species.\n", sp))
	 print("Check for an appropriate file in:\n")
	 e3bin = os.getenv("E3BIN") or os.getenv("HOME").."/e3bin"
	 dir = e3bin.."/species"
	 print("   ", dir)
	 print("Bailing out!\n")
	 os.exit(1)
      end
      f:write(string.format("'%s', ", sp))
      species2[#species2+1] = sp
   end
   f:write("}\n")
   f:write("\n")
   
   e3bin = os.getenv("E3BIN") or os.getenv("HOME").."/e3bin"
   dir = e3bin.."/species"
   default_file = dir.."/defaults.lua"
   dofile(default_file)

   for _,sp in ipairs(species2) do
      e3bin = os.getenv("E3BIN") or os.getenv("HOME").."/e3bin"
      dir = e3bin.."/species/"
      file = dir..sp..".lua"
      dofile(file)
      
      f:write(string.format("%s = {}\n", sp))

      key = sp..".atomic_constituents"
      f:write(string.format("%s = {\n", key))
      if _G[sp]['atomic_constituents'] then
	 for k,v in pairs(_G[sp]['atomic_constituents']) do
	    f:write(string.format("%s=%d, ", k,v))
	 end
      else
	 print("WARNING: No 'atomic_constituents' set for species: ", sp)
	 print("No constituents will be listed in the gas model file.")
      end
      f:write("\n}\n")

      key = sp..".charge"
      if _G[sp]['charge'] then
	 f:write(string.format("%s = %d\n", key, _G[sp]['charge']))
      else
	 print("WARNING: No 'charge' value set for species: ", sp)
	 print("A default value of ", default['charge'], " will be written")
	 print("to the gas model file.")
	 f:write(string.format("%s = %d\n", key, default['charge']))
      end

      for __,val in ipairs(t.value_list) do
	 var = sp.."."..val
	 f:write(string.format("%s = ", var))
	 if _G[sp][val] then
	    serialise(_G[sp][val], f)
	 else
	    serialise(default[val], f)
	 end
	 f:write("\n")
      end

      for __,model in ipairs(t.model_list) do
	 var = sp.."."..model
	 f:write(string.format("%s = ", var))
	 if _G[sp][model] then
	    serialise(_G[sp][model], f)
	 else
	    serialise(default[model], f)
	 end
	 f:write("\n")
      end
      
      models_w_CEA_coeffs = {'thermally perfect', 'thermally real'}
      for _,v in ipairs(models_w_CEA_coeffs) do
	 if t.TBM == v then
	    var = sp..".CEA_coeffs"
	    f:write(string.format("%s = ", var))
	    serialise(_G[sp].CEA_coeffs, f)
	    f:write("\n")
	 end
      end
   end
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
   if not model then
      print("A gas model must be specified in the input file. For example:")
      print("  model = 'ideal gas'")
      print("Bailing out!")
      os.exit(1)
   end
   f = assert(io.open(output, "w"))
   -- print(string.format("%s model selected", model))
   if gas_models[model] ~= nil then
      if gas_models[model].special_fn then
	 gas_models[model].special_fn(species, f)
      else
	 create_gas_file(species, gas_models[model], f)
      end
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
