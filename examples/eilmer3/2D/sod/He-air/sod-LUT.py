# Sod's shock tube problem with air from LUT, and He treated as ideal
#
# Adapted from examples/mbcns2/sod2
#
# Author: Rowan J. Gollan
# Date: 16-Feb-2009

gdata.title = "One-dimensional shock tube with helium driving air."
gdata.dimensions = 2
gdata.stringent_cfl = 1

create_gas_file(model="ideal gas", species=['He',], 
                fname="LUT-plus-He.lua", lut_file="cea-lut-air.lua.gz")
species_list = select_gas_model(fname="LUT-plus-He.lua")
print "species_list=", species_list

helium = FlowCondition(p=1.0e5, u=0.0, v=0.0, T=348.4, massf=[0.0, 1.0])
air = FlowCondition(p=1.0e4, u=0.0, v=0.0, T=278.7, massf=[1.0, 0.0])

# Set up two quadrilaterals in the (x,y)-plane by first defining
# the corner nodes, then the lines between those corners.
# The labelling is not significant; it is just to make the MetaPost
# picture look the same as that produced by the Tcl scriptit program.
a = Node(0.5, 0.0, label="a")
b = Node(0.5, 0.1, label="b")
c = Node(0.0, 0.1, label="c")
d = Node(0.0, 0.0,  label="d")
e = Node(1.0, 0.0, label="e")
f = Node(1.0, 0.1, label="f")

south0 = Line(d, a); south1 = Line(a, e) # lower boundary along x-axis
north0 = Line(c, b); north1 = Line(b, f) # upper boundary
west0 = Line(d, c); east0west1 = Line(a, b); east1 = Line(e, f)

# Define the blocks, boundary conditions and set the discretisation.
nx = 50; ny = 2
blk_0 = Block2D(make_patch(north0, east0west1, south0, west0), nni=nx, nnj=ny,
                fill_condition=helium, label="driver")
blk_1 = Block2D(make_patch(north1, east1, south1, east0west1), nni=nx, nnj=ny,
                fill_condition=air, label="driven")
identify_block_connections()

# Some simulation parameters
gdata.flux_calc = AUSMDV
gdata.max_time = 0.4e-3
gdata.max_step = 600
gdata.dt = 1.0e-6


