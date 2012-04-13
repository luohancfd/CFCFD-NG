# High-performance shock tube with helium driving air 
# in a constant-diameter tube.  The temperatures in the air
# are high enough to induce strong thermochemical effects.
#
# Adapted from examples/mbcns2/sod2/, examples/eilmer3/2D/sod/He-air/.
# Authors: PAJ and RJG
# Date: 24-Mar-2012

gdata.title = "High-performance shock tube with helium driving air."
# Combine a LUT air model with a composite gas of pure helium.
create_gas_file(model="ideal gas", species=['He',], 
                fname="LUT-plus-He.lua", lut_file="cea-lut-air.lua.gz")
species_list = select_gas_model(fname="LUT-plus-He.lua")
print "species_list=", species_list

helium = FlowCondition(p=30.0e6, T=3000, massf={'He':1.0})
air = FlowCondition(p=30.0e3, T=300.0, massf={'LUT':1.0})

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
                fill_condition=helium, label="driver")
blk_1 = Block2D(make_patch(north1, east1, south1, east0west1), 
                nni=400, nnj=2,
                fill_condition=air, label="driven")
identify_block_connections()

# Some simulation parameters
gdata.flux_calc = ADAPTIVE
gdata.max_time = 100.0e-6
gdata.max_step = 8000
gdata.dt = 1.0e-9

sketch.xaxis(0.0, 1.0, 0.5, -0.05)
sketch.yaxis(0.0, 0.1, 0.1, -0.05)
sketch.window(0.0, 0.0, 1.0, 0.1, 0.02, 0.02, 0.17, 0.035)


