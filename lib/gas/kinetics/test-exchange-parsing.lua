require 'lpeg'
require 'exch'

print("Testing 'collider' grammar.")
local CG = exch.ColliderG

local t1 = "N2 ~~ O2"
local t2 = "N2 ~~ (N2, O2, H2, N2)"
local t3 = "N2 ~~ (*list)"
local t4 = "N2 ~~ ( N2_4So, N2_X, CO_X, CO2 )"

local tests = {t1, t2, t3, t4}

for _,t in ipairs(tests) do
   result = lpeg.match(CG, t)
   if result then
      for k,v in pairs(result) do
	 print(k,v)
	 if type(v) == 'table' then
	    for __,sp in ipairs(v) do
	       print(sp)
	    end
	 end
      end
      print("")
   else
      print("No match for: ", t)
      print("")
   end
end
print("Done.\n")

--[=[
print("Testing 'exchange' grammar.")
local EG = exch.ExchangeG

t1 = "V-T"
t2 = "V - V"
t3 = "E-T"
t4 = "T-W"
t5 = "V-V THO"

tests = {t1, t2, t3, t4, t5}
for _,t in ipairs(tests) do
   result = lpeg.match(EG, t)
   if result then
      for i,et in ipairs(result) do
	 print(i, et)
      end
      print("")
   else
      print("No match for: ", t)
      print("")
   end
end
print("Done.\n")

print("Testing complete 'mechanism' grammar.")
local MG = exch.ThermMechG

t1 = "N2 ~~ O2 : V-V"
t2 = "N2 ~~ (O2, H2) : V-T"
t3 = "O2 ~~ (*list) : E-T"
t4 = "N2 ~~ (O2, N2) : V-V THO"

tests = {t1, t2, t3, t4}

for _,t in ipairs(tests) do
   result = lpeg.match(MG, t)
   if result then
      for it,tab in ipairs(result) do
	 print("tab: ", it)
	 for i, e in ipairs(tab) do
	    print(i, e)
	 end
      end
      print("")
   else
      print("No match for: ", t)
      print("")
   end
end
print("Done.\n")


   --]=]
