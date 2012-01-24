-- diaphragm_1.lua
-- User-defined parameters controlling the primary diaphragm opening.
-- The working functions are actually in diaphragm_common.lua.

p_burst = 500.0e3
-- diaphragm bursts once a single cell pressure difference
-- on interface exceeds this rupture pressure. 

is_burst = true
-- =false, initially diaphragm is not burst;
-- =true, diaphragm begins to rupture immediately.

t_trigger = 0.0e-6
-- This is the trigger time for diaphragm.
-- This value only used to initialise variable.
t_hold = 0.0e-6
-- Hold time for diaphragm after first cell triggers rupture.
-- This is a delay before it begins opening.
dt_burst = 1.0e-4
-- Time taken for diaphragm to fully open up after initiation of rupture. 
-- This ensures area of opening increases linearly with time,
-- not radius (i.e. for axisymmetric analysis). 
-- Opens according to formula y_current=((t-t_trigger)/dt_burst*ymax^2),
-- where y==radius.
-- Change this for 2D non-axisymmetric simulation!

d_tube = 0.085
-- Tube diameter (assuming axisymmetric analysis).
y_max = d_tube/2
-- Axisymmetric, so remember to divide by 2.
-- This is vertical position of diaphragm opening ending point. 

y_current=-1.0
-- current position of diaphragm opening; -1.0 means closed.
-- Note, currently this only works for a vertically aligned diaphragm. 

dofile("diaphragm_common.lua")
