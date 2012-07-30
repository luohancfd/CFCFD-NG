require 'lpeg'
require 'reac'

local SG = reac.SpeciesG
local S2G = reac.Species2G

local t1 = "N"
local t2 = "N2"
local t3 = "Co"
local t4 = "N2+"
local t5 = "N2O2-"
local t6 = "N2_S"
local t7 = "N2_X"

tests = {t1, t2, t3, t4, t5, t6, t7}

for _,t in ipairs(tests) do
   print(lpeg.match(SG, t))
end

for _,t in ipairs(tests) do
   tt = lpeg.match(S2G, t)
   for k,v in pairs(tt) do
      print(k, v)
      if type(v) == 'table' then
	 for k1,v1 in pairs(v) do
	    print(k1, v1)
	 end
      end
   end

end


