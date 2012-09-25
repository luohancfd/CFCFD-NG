-- Author: Rowan J. Gollan
-- Date: 29-Aug-2012
-- Place: Greenslopes, Queensland, Australia
--

module(..., package.seeall)

-- For the moment, assume translation mode is always 0
ITRANS = 0

-- lexical elements
local Space = lpeg.S(" \n\t")^0
local Number = lpeg.R("09")^1
local Underscore = lpeg.S("_")
local Element = ((lpeg.R("AZ") * lpeg.R("az")^0) + lpeg.P("e"))
local Solid = lpeg.P("S")
local ElecLevel = lpeg.R("%w")^(-3) -- %w matches alphanumeric characters
                                    -- ^(-3) says to match at most 3 occurrences
local PM = lpeg.S("+-")
Species = lpeg.C(((Element * Number^0)^1 * PM^0)^1 * (Underscore * (Solid + ElecLevel))^0) * Space
local Tilde = lpeg.P("~")
local Dash = lpeg.P("-") * Space
local Comma = lpeg.P(",") * Space
local Colon = lpeg.P(":") * Space
local Open = "(" * Space
local Close = ")" * Space
local ListKw = lpeg.C(lpeg.P("*list")) * Space
local ExchType = lpeg.C(lpeg.S("VTE")) * Space -- letters designating exchange types

-- grammar
local Colliders = lpeg.V"Colliders"
local Exchange = lpeg.V"Exchange"
local Mechanism = lpeg.V"Mechanism"

local G = lpeg.P{ Mechanism,
		  Mechanism = lpeg.Ct(Colliders * Colon * Exchange),
		  Colliders = lpeg.Ct(Species * Space * (Tilde*Tilde) * Space *
					( Species + -- e.g. N2 ~~ O2
					 lpeg.Ct((Open*Species*Space*(Comma*Species)^0*Close)) +  -- e.g. O2 ~~ (N2, H2)
					 (Open*ListKw*Close) )),
		  Exchange = lpeg.Ct(ExchType * Dash * ExchType)
		}
ThermMechG = Space * G * -1
ColliderG = lpeg.P{ Colliders,
		    Colliders = lpeg.Ct(Species * Space * (Tilde*Tilde) * Space *
					( Species + -- e.g. N2 ~~ O2
					   lpeg.Ct((Open*Species*Space*(Comma*Species)^0*Close)) +  -- e.g. O2 ~~ (N2, H2)
					   (Open*ListKw*Close) )) } -- e.g. N2 + (*list)
ExchangeG = lpeg.P{ Exchange,
		    Exchange = lpeg.Ct(ExchType * Dash * ExchType)}


function parse_energy_exch_string(s)
   local t = lpeg.match(ThermMechG, s)
   return t
end

function transform_mechanism(m, species, thermal_modes)
   local t = parse_energy_exch_string(m[1])
   local ip = species[t[1][1]]
   local q = t[1][2]
   local iqs = {}

   if type(q) == 'table' then
      for _,sq in ipairs(q) do
	 iqs[#iqs+1] = species[sq]
      end
   elseif type(q) == 'string' then
      if species[q] then
	 -- We found q as a species entry
	 iqs[#iqs+1] = species[q]
      elseif q == '*list' then
	 if not m.list then
	    print("Keyword '*list' has been used to designate colliders")
	    print("but a list entry is not available in the mechanism.")
	    print("Bailing out!")
	    os.exit(1)
	 end
	 for _,sq in ipairs(m.list) do
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

   for i,q in ipairs(iqs) do
      mechs[#mechs+1] = {ip=ip, iq=q, type=m_type, relaxation_time=m.rt}
      if m_type == 'VT' then
	 mechs[#mechs].itrans = ITRANS
	 -- Need to do this more intelligently
	 mechs[#mechs].imode = 1
      end
   end

   return mechs

end


