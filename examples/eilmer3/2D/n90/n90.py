# n90.py
# Updated version of n90.sit (an mb_cns example)
#
# R.J.Gollan
# Updated on 12-Mar-2008
# PJ - Elmer3 port, July 2008
# RJG - updated for new kinetics library, Nov 2008
# PJ - different ways to set mass fractions, July 2009

gdata.title = "Cylinder in Mach 10 nitrogen flow"

# Gas model selection
species_list = select_gas_model(model='thermally perfect gas', 
                                species=['N2', 'N'])
print "species_list=", species_list
set_reaction_scheme("nitrogen-2sp-2r.lua",reacting_flag=1)

# Flow conditions
# mf = [1.0, 0.0]
# mf = 1.0
# mf = 1
mf = {'N2':1.0}
inflow = FlowCondition(p=500.0, u=5000.0, v=0.0, T=700.0, massf=mf)
initial = FlowCondition(p=5.0, u=0.0, v=0.0, T=300.0, massf=mf)
print "inflow=", inflow
print "initial=", initial

# Geometry
a = Node( 0.045,    0.0,       label="a")
b = Node( 0.0,      0.0,       label="b")
c = Node( 0.013181, 0.031820,  label="c")
d = Node( 0.045,    0.045,     label="d")
e = Node( 0.0675,   0.038972,  label="e")
f = Node(-0.020,    0.0,       label="f")
g = Node(-0.020,    0.050625,  label="g")
h = Node(-0.016875, 0.106875,  label="h")
i = Node( 0.045,    0.135,     label="i")
j = Node( 0.07875,  0.095625,  label="j")
k = Node( 0.084375, 0.0675,    label="k")

bc = Arc(b, c, a, "ab")
cd = Arc(c, d, a, "cd")
de = Arc(d, e, a, "de")
east = Polyline([bc, cd, de], "east")
west = Bezier([f, g, h, i], "west")
south = Line(f, b, "south")
north = Bezier([i, j, k, e], "north")

# Block setup
NNR = 60
NNT = 40

blk = Block2D(make_patch(north, east, south, west),
              nni=NNR, nnj=NNT,
              fill_condition=initial)
blk.set_BC(WEST, SUP_IN, inflow)
blk.set_BC(NORTH, EXTRAPOLATE_OUT)

# Simulation parameters
gdata.flux_calc = ADAPTIVE
gdata.max_time = 100.0e-6
gdata.max_step = 40000
gdata.dt = 1.0e-8
gdata.cfl = 0.5
gdata.dt_plot = 20.0e-6

sketch.xaxis(-0.02, 0.10, 0.02, -0.005)
sketch.yaxis(0.0, 0.14, 0.02, -0.005)
sketch.window(-0.02, 0.0, 0.10, 0.12, 0.05, 0.05, 0.17, 0.17)

