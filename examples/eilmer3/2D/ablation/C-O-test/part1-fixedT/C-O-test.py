## \file C-O-test.py
## \brief oxidation of a spherical graphite model
##        testing reactions at the wall.
## part1: fixed T wall, no wall reactions. Generate a
##        solution to start reacting wall case from.
## \author EJF, 23-Sep-2013

from math import cos, sin, tan, sqrt, pi
from cfpylib.grid.shock_layer_surface import *

job_title = "Porous test case."
print job_title

gdata.title = job_title
gdata.axisymmetric_flag = 1

#
# 1. Setup the gas model
#
species = select_gas_model(model='two temperature gas', species=['C','O2','O','CO','C2','C3','C_plus','O_plus','e_minus'])
set_reaction_update("Park2001-C-O-flowfield.lua", reacting_flag=1)
gm = get_gas_model_ptr()
massf_inf = [ 0.0 ] * gm.get_number_of_species()
massf_inf[species.index('O2')] = 1.0

#
# 2. Define flow conditions
#
inflow = FlowCondition(p=500.0, u=5000.0, v=0.0, T=2500.0, massf=massf_inf)
initial = FlowCondition(p=5.0, u=0.0, v=0.0, T=300.0, massf=massf_inf)

#
# 3. Define the geometry
#
a = Node( 0.045,    0.0,       label="a")
b = Node( 0.0,      0.0,       label="b")
c = Node( 0.013181, 0.031820,  label="c")
d = Node( 0.045,    0.045,     label="d")
e = Node( 0.0675,   0.038972,  label="e")
f = Node( -0.020,    0.0,       label="f")
g = Node( -0.020,    0.050625,  label="g")
h = Node( -0.016875, 0.106875,  label="h")
i = Node( 0.045,    0.135,     label="i")
j = Node( 0.07875,  0.095625,  label="j")
k = Node( 0.084375, 0.0675,    label="k")

bc = Arc(b, c, a, "ab")
cd = Arc(c, d, a, "cd")
de = Arc(d, e, a, "de")
east = Polyline([bc, cd, de], "east") # wall
west = Bezier([f, g, h, i], "west") # inflow
south = Line(f, b, "south") # symmetry axis
north = Bezier([i, j, k, e], "north") # outflow
#
# 4. Define the blocks, boundary conditions and set the discretisation
#
nnx = 40; nny= 40
nbx = 2; nby = 2

blk_0 = SuperBlock2D(make_patch(north,east,south,west), nni=nnx, nnj=nny, nbi=nbx, nbj=nby,
                bc_list=[ExtrapolateOutBC(), FixedTBC(3500.0), SlipWallBC(), SupInBC(inflow)],
                fill_condition=initial, label="BLOCK-0")
# N = extr, E = extr, S = wall, W = inflow UserDefinedBC(filename="udf-char-mass-in.lua",sets_conv_flux=1,is_wall=1)
identify_block_connections()

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.diffusion_flag = 1
gdata.diffusion_model = "ConstantLewisNumber"
gdata.max_time = 5.0e-4  # seconds
gdata.max_step = 100000
gdata.dt = 1.0e-10
gdata.dt_max = 1.0e-7
gdata.dt_plot = gdata.max_time/20.0
gdata.cfl = 0.25

sketch.xaxis(0.0, 1.0, 0.2, -0.05)
sketch.yaxis(0.0, 1.0, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)

