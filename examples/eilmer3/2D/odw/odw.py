# odw.py
#
# A Python input file to describe the oblique
# detonation wave used as a verification case.
#
# Thie Python file prepared by...
# Rowan J Gollan
# and adjusted for the new geometry spec by PJ (Aug-06)
# 17-May-2009: updated for Eilmer3 by RJG

gdata.title = "The oblique detonation wave verification case."
select_gas_model(fname="binary-gas.lua")

inflow = FlowCondition(p=86.1e3, u=964.302, v=0.0, T=[300.0], massf=[1.0, 0.0])
initial = FlowCondition(p=28.7e3, u=0.0, v=0.0, T=[300.0], massf=[1.0, 0.0])

#
# Geometry
xmin = -0.25
xmax = 1.75
ymin = 0.0
ymax = 2.0

nnx = 64
nny = 64

from oblique_detonation import *
from math import pi
od =  ObliqueDetonation( pi/4.0, 300.0, 3.0, 1.0)
wall = od.create_wall_spline(0.0, xmax, 1.0e-5)

a = Node( xmin, 0.0, label="a")
b = Node( 0.0, 0.0, label="b")
c = Node( wall.eval(1.0).x, wall.eval(1.0).y, label="c" )
d = Node( xmin, ymax, label="d")
e = Node( 0.0, ymax, label="e")
f = Node( xmax, ymax, label="f")

south0 = Line(a, b)
west0 = Line(a, d)
south1 = wall
# bc = Line(b, c) # straighten the wall
east0west1 = Line(b, e)
east1 = Line(c, f)
north0 = Line(d, e)
north1 = Line(e, f)

nnx0 = int( 0.125 * nnx )
nnx1 = nnx - int( 0.125 * nnx )

blk_0 = Block2D(
    psurf=make_patch(north0, east0west1, south0, west0),
    nni=nnx0, nnj=nny,
    bc_list=[ExtrapolateOutBC(), AdjacentBC(), SlipWallBC(), SupInBC(inflow)],
    fill_condition=inflow,
    label="blk-0"
    )
blk_1 = Block2D(
    psurf=make_patch(north1, east1, south1, east0west1),
    nni=nnx1, nnj=nny,
    bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(), SlipWallBC(), AdjacentBC()],
    fill_condition=inflow,
    label="blk-1"
    )
identify_block_connections()

# Simulate the reaction between reactants
# to form products by giving an appropriate
# user-defined source vector

gdata.udf_file = "udf-source.lua"
gdata.udf_source_vector_flag = 1

# Do a little more setting of global data.
gdata.flux_calc = AUSMDV
gdata.max_time = 1.0e-2  # seconds
gdata.max_step = 30000
gdata.dt = 1.0e-6
gdata.dt_plot = gdata.max_time/20.0
gdata.dt_history = 10.0e-5

# Values to make the SVG look good
sketch.xaxis(-0.25, 1.75, 0.25, -0.06)
sketch.yaxis(0.0, 2.0, 0.25, -0.06)
sketch.window(-0.25, 0.0, 1.75, 2.0, 0.05, 0.05, 0.17, 0.17)
