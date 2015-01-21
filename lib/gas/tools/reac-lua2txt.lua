#!/usr/bin/env lua
-- Author: Rowan J. Gollan
-- Date: 22-Jan-2015
--

require 'reaction_parser'

local process = reaction_parser.main

reactions = reaction_parser.reactions
scheme_t = reaction_parser.scheme_t
ode_t = reaction_parser.ode_t

scheme = reaction_parser.scheme
ode_solver = reaction_parser.ode_solver
reaction = reaction_parser.reaction
remove_reactions_by_label = reaction_parser.remove_reactions_by_label
remove_reactions_by_number = reaction_parser.remove_reactions_by_number
select_reactions_by_label = reaction_parser.select_reactions_by_label
select_reactions_by_number = reaction_parser.select_reactions_by_number

function print_help()
   print "reac-lua2txt -- a tool to convert Lua chemistry files to a text description."
   print "Usage:"
   print "> reac-lua2txt gas-model.lua chemistry-inp.lua output.txt"
   print ""
   print "where:"
   print "  gas-model.lua -- name of a Lua gas model file."
   print "  chemistry-input.lua  -- name of Lua chemistry input file"
   print "  output.txt -- name of the output text file"
   print ""
   return
end

function main()
   if ( #arg ~= 3 ) then
      print("Incorrect number of arguments. Exactly 3 arguments expected.")
      print_help()
      os.exit(1)
   end
   
   -- 1. We are interested in the species table from
   --    the gas model file.
   gfile = arg[1]
   dofile(gfile)
   -- Make species file contain reverse lookup
   nsp = #species
   for isp,sp in ipairs(species) do
      species[sp] = isp - 1
   end
   species.size = nsp 

   -- 2. Next we parse the reactions file
   input = arg[2]
   process(input)

   -- 3. Now we can loop over the reactions writing out the text file
   f = assert(io.open(arg[3], 'w'))
   f:write(string.format("%d %d\n", #species, #reactions))
   f:write("#Reactant Stoichiometric Matrix\n")
   for ir,r in ipairs(reactions) do
      r_coeffs = {}
      for _,t in pairs(r.f_coeffs) do
	 r_coeffs[t[1]] = t[2]
      end
      r_str = ""
      for isp=0,nsp-1 do
	 if r_coeffs[isp] then
	    r_str = r_str .. string.format("%.1f ", r_coeffs[isp])
	 else
	    r_str = r_str .. "0.0 "
	 end
      end
      r_str = r_str .. "\n"
      f:write(r_str)
   end
   f:write("#Product Stoichiometric Matrix\n")
   for ir,r in ipairs(reactions) do
      p_coeffs = {}
      for _,t in pairs(r.b_coeffs) do
	 p_coeffs[t[1]] = t[2]
      end
      p_str = ""
      for isp=0,nsp-1 do
	 if p_coeffs[isp] then
	    p_str = p_str .. string.format("%.1f ", p_coeffs[isp])
	 else
	    p_str = p_str .. "0.0 "
	 end
      end
      p_str = p_str .. "\n"
      f:write(p_str)
   end
   f:write("#Third body reactions? No : 0, Yes : 1 or 2\n")
   for ir,r in ipairs(reactions) do
      if not r.third_body then
	 f:write("0\n")
      else
	 all_efficiencies_are_one = true
	 for isp=1,nsp do
	    if r.efficiencies[isp][2] ~= 1.0 then
	       all_efficiencies_are_one = false
	    end
	 end
	 if all_efficiencies_are_one then
	    f:write("1\n")
	 else
	    f:write("2 ")
	    for isp=1,nsp do
	       if r.efficiencies[isp] then
		  f:write(string.format("%12.6f ", r.efficiencies[isp][2]))
	       else
		  f:write(string.format("%12.6f ", 0.0))
	       end
	    end
	    f:write("\n")
	 end
      end
   end
   f:close()
end


main()

