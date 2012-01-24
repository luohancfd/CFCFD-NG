-- diaphragm_2.lua
-- Secondary diaphragm boundary condition.
-- Comments are in diaphragm_1.lua.

p_burst = 200.0e3
is_burst = false
t_trigger = 0.0e-6
t_hold = 10.0e-6
dt_burst = 2.0e-5
d_tube = 0.085
y_max = d_tube/2
y_current = -1.0
dofile("diaphragm_common.lua")
