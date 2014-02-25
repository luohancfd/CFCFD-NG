# Sod's shock tube problem with a mixture of ideal gases
#
# Description
# -----------
# A mixture of two species N2 and O2 with mass fractions: 0.7778 and 0.2222
# respectively is present. This is roughly 80% N2 and 20% O2 by volume which is
# representative of room temperature air.
#
# This test case has been adapted from an old mb_cns test case:
#  sod_mix/test_case_1.sit
#
# Author: Rowan J. Gollan
# Date: 30-July-2008
#
# 25-Feb-2014 PJ temperature monitoring examples

gdata.title = "One-dimensional shock tube with (fake) air driving (fake) air."
gdata.dimensions = 2
gdata.stringent_cfl = 1
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.10

select_gas_model(model='ideal gas',
                 species=['N2', 'O2'])

high_p = FlowCondition(p=1.0e5, u=0.0, v=0.0, T=348.4, massf=[0.7778, 0.2222])
low_p = FlowCondition(p=1.0e4, u=0.0, v=0.0, T=278.7, massf=[0.7778, 0.2222])

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
                fill_condition=high_p, label="driver")
blk_1 = Block2D(make_patch(north1, east1, south1, east0west1), nni=nx, nnj=ny,
                fill_condition=low_p, label="driven", mcell_list=[(nx/2,0),])
identify_block_connections()

# Some simulation parameters
gdata.flux_calc = AUSMDV
gdata.max_time = 0.6e-3
gdata.max_step = 600
gdata.dt = 1.0e-6

# Exercise the halt of large temperature-change code
# for Elise and Dan's radiation-coupled calculations.
MonitorLocation(0.55,0.0)
MonitorLocation(0.25,0.1)
gdata.halt_on_large_flow_change = True
gdata.tolerance_in_T = 20.0


