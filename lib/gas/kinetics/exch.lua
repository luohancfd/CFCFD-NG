-- Author: Rowan J. Gollan
-- Date: 29-Aug-2012
-- Place: Greenslopes, Queensland, Australia
--

module(..., package.seeall)

require 'reac'
require 'lex_elems'

local transform_species_str = reac.transform_species_str

-- For the moment, assume translation mode is always 0
ITRANS = 0

-- lexical elements
-- get common elements from lex_elems.lua
for k,v in pairs(lex_elems) do
   _G[k] = v
end
-- List some specific elements for this module
local ListKw = lpeg.C(lpeg.P("*list")) * Space
local AllKw = lpeg.C(lpeg.P("*all")) * Space
local ExchType = lpeg.C(lpeg.S("VTE")) * Space -- letters designating exchange types
local ExchDesc = lpeg.C(lpeg.R("AZ")^1)

-- grammar
local Colliders = lpeg.V"Colliders"
local Exchange = lpeg.V"Exchange"
local Mechanism = lpeg.V"Mechanism"

local G = lpeg.P{ Mechanism,
		  Mechanism = lpeg.Ct(Colliders * Space * Colon * Space * Exchange),
		  Colliders = lpeg.Ct(Species * Space * (Tilde*Tilde) * Space *
					( Species + -- e.g. N2 ~~ O2
					 lpeg.Ct((Open*Species*(Comma*Species)^0*Close)) +  -- e.g. O2 ~~ (N2, H2)
					 (Open*ListKw*Close) + -- e.g. N2 ~~ (*list)
					 (Open*AllKw*Close) )), -- e.g. N2 ~~ (*all)
		  Exchange = lpeg.Ct(ExchType * Dash * ExchType*(Space*ExchDesc)^0)
		}
ThermMechG = Space * G * -1
ColliderG = lpeg.P{ Colliders,
		    Colliders = lpeg.Ct(Species * Space * (Tilde*Tilde) * Space *
					( Species + -- e.g. N2 ~~ O2
					   lpeg.Ct((Open*Species*(Comma*Species)^0*Close)) +  -- e.g. O2 ~~ (N2, H2)
					   (Open*ListKw*Close) + -- e.g. N2 ~~ (*list)
					   (Open*AllKw*Close) )) } -- e.g. N2 ~~ (*all) 
ExchangeG = lpeg.P{ Exchange,
		    Exchange = lpeg.Ct(ExchType * Dash * ExchType*(Space*ExchDesc)^0)}


function parse_energy_exch_string(s)
   local t = lpeg.match(ThermMechG, s)
   if t == nil then
      print("Problem parsing energy exchange string:")
      print(s)
      print("Please check input for energy exchange mechanisms.")
      print("Bailing out!")
      os.exit(1)
   end
   return t
end

local function find_mode_index(modes, mode_idx, e_str, s_str)
   -- 1. First based on species name
   local cname1
   if e_str == "V" then
      cname1 = s_str.."-vibration"
   elseif e_str == "E" then
      cname1 = "e_minus-translation"
   else
      print(string.format("Mode string %s" % e_str))
      print("is not known in energy exchange file.")
      print("Bailing out!")
      os.exit(1)
   end

   for k,t in pairs(modes) do
      for c,_ in pairs(t) do
	 if cname1 == c then
	    return mode_idx[k]
	 end
      end
   end

   if e_str == "E" then
      -- If we're here we did NOT find "e_minus-translation"
      print("Problem finding the mode 'e_minus-translation")
      print("as one of the thermal modes while trying to")
      print("setup the energy exchange mechanisms.")
      print("Bailing out!")
      os.exit(1)
   end

   -- 2. Search for 'all' in case this species is lumped with others
   -- We've already searched for the special case of e-
      
   local cname2 = "all-vibration"
   for k,t in pairs(modes) do
      for c,_ in pairs(t) do
	 if cname2 == c then
	    return mode_idx[k]
	 end
      end
   end

   -- 3. If we got this far, we did NOT find this mode.
   print("There was a problem determining the index of mode:")
   print(string.format("%s-%s", e_str, s_str))
   print("The following mode names were searched for:")
   print(cname1)
   print(cname2)
   os.exit(1)
end

function transform_mechanism(m, species, thermal_modes)
   local t = parse_energy_exch_string(m[1])
   print("t[1][1]= ", t[1][1])
   local p = transform_species_str(t[1][1]) -- species p as string
   local ip = species[p]
   local q = t[1][2] -- species q(s): could be string, list, or table
   local iqs = {}

   if type(q) == 'table' then
      for _,sq in ipairs(q) do
	 local sq_t = transform_species_str(sq)
	 iqs[#iqs+1] = species[sq]
      end
   elseif type(q) == 'string' then
      q_t = transform_species_str(q)
      if species[q_t] then
	 -- We found q_t as a species entry
	 iqs[#iqs+1] = species[q_t]
      elseif q == '*list' then
	 if not m.list then
	    print("Keyword '*list' has been used to designate colliders")
	    print("but a list entry is not available in the mechanism.")
	    print("Bailing out!")
	    os.exit(1)
	 end
	 for _,sq in ipairs(m.list) do
	    local sq_t = transform_species_str(sq)
	    iqs[#iqs+1] = species[sq]
	 end
      else
	 print("Keyword used for colliders is unknown: ", q)
	 print("Bailing out!")
	 os.exit(1)
      end
   else
      print("Don't know what to do with q of type: ", type(q))
      print("Bailing out!")
      os.exit(1)
   end

   if not m.rt then
      print("A relaxation time model 'rt={...}' has not been declared.")
      print("Bailing out!")
      os.exit(1)
   end

   -- User list the 'model' first
   m.rt.model = m.rt[1]

   local mechs = {}
   local m_type = tostring(t[2][1]..t[2][2])
   if t[2][3] then
      m_type = m_type.."-"..t[2][3]
   end
   local e_str = t[2][1]
   local s_str = p
   local imode = find_mode_index(thermal_modes, mode_idx, e_str, s_str)

   if ( imode == nil ) then
      print("A mode for thermal exchange listed in the energy exchange file")
      print("is NOT available as part of the gas model.")
      print("mode : ", mode)
      print("This was found when building exchange mechanism:")
      print(m[1])
      print("")
      print("Bailing out!")
      os.exit(1)
   end

   for i,q in ipairs(iqs) do
      mechs[#mechs+1] = {ip=ip, iq=q, type=m_type, relaxation_time=m.rt}
      if m_type == 'VT' then
	 mechs[#mechs].itrans = ITRANS
	 mechs[#mechs].imode = imode
      elseif m_type == "VV-THO" then
	 mechs[#mechs].itrans = ITRANS
	 mechs[#mechs].imode = imode
	 local qmode = "V_"..species[q]
	 mechs[#mechs].iTvq = thermal_modes[qmode]
      elseif m_type == "ET" then
	 mechs[#mechs].itrans = ITRANS
	 mechs[#mechs].imode = imode
      else
	 print("Mode type: ", m_type, " is not known.")
	 print("This occurred when trying to build energy exchange mechanism:")
	 print(m[1])
	 print("")
	 print("Bailing out!")
	 os.exit(1)
      end
   end

   return mechs

end


