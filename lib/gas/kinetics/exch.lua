-- Author: Rowan J. Gollan
-- Date: 29-Aug-2012
-- Place: Greenslopes, Queensland, Australia
--

module(..., package.seeall)

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

			
		-- * Space * ( Species + (Open Species^1 * (Comma * Species)^0 * Comma^0 Close) +  (Open ListKw Close) ) * Space,					       
