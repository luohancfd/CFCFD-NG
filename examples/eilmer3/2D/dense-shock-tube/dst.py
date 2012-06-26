"""
dst.py -- A dense real-gas shock tube example. Adapted from the classic-shock-tube
example and gasdyn example.

Further work is required to compare this simulation to that of:

Argrow, BM (1996). 'Computational Analyisis of Dense Gas Shock Tube Flow',
Shock Waves, vol 6 num 4, 241-248, Springer

.. Author: Peter Blyton
.. Version: 26/06/2012
"""

gdata.title = "Dense gas shock tube with CO2."
select_gas_model(fname="REFPROP-lut-CO2.FLD.lua.gz")

driver = FlowCondition(p=6.0e6, T=1000.0)
driven = FlowCondition(p=30.0e3, T=300.0)

a = Node(0.5, 0.0, label="a"); b = Node(0.5, 0.1, label="b")
c = Node(0.0, 0.1, label="c"); d = Node(0.0, 0.0,  label="d")
e = Node(1.0, 0.0, label="e"); f = Node(1.0, 0.1, label="f")
south0 = Line(d, a); south1 = Line(a, e) # lower boundary along x-axis
north0 = Line(c, b); north1 = Line(b, f) # upper boundary
# left-end, diaphragm, right-end
west0 = Line(d, c); east0west1 = Line(a, b); east1 = Line(e, f) 

# Define the blocks, boundary conditions and set the discretisation.
blk_0 = Block2D(make_patch(north0, east0west1, south0, west0), 
                nni=400, nnj=2,
                fill_condition=driver, label="driver")
blk_1 = Block2D(make_patch(north1, east1, south1, east0west1), 
                nni=400, nnj=2,
                fill_condition=driven, label="driven")
identify_block_connections()

# Some simulation parameters
gdata.flux_calc = ADAPTIVE
gdata.max_time = 100.0e-6
gdata.max_step = 8000
gdata.dt = 1.0e-9

sketch.xaxis(0.0, 1.0, 0.5, -0.05)
sketch.yaxis(0.0, 0.1, 0.1, -0.05)
sketch.window(0.0, 0.0, 1.0, 0.1, 0.02, 0.02, 0.17, 0.035)
