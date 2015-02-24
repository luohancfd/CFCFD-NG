-- test.lua
-- An example of driving the wrapped classes and functions
-- used in preparing a flow simulation.

-- For the moment, just use the example from Rowan's luaflowstate_demo.
fs = FlowState:new{p=1.0e5, T=300.0, u=1000.0, v=200.0}
fsTab = fs:toTable{}
for k,v in pairs(fsTab) do
    print(k,v)
    if ( k == 'gas' ) then
       for k1,v1 in pairs(v) do
           print(k1, v1)
           if ( k1 == 'massf' ) then
               print("massf[1]= ", v1[1])
           end
       end
    end
end
