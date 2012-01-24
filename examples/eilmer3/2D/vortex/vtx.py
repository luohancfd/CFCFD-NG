# file: vtx.py
# PJ, 14-Dec-2006
#     01-Feb-2010 ported to Eilmer3
gdata.title = "Inviscid supersonic vortex -- flow in a bend."

# Geometry
R_inner = 1.0
R_outer = 1.384
a = Node(0.0, 0.0)
b = Node(0.0, R_inner)
c = Node(0.0, R_outer)
d = Node(R_inner, 0.0)
e = Node(R_outer, 0.0)
north0 = Arc(c, e, a)
east0 = Line(d, e)
south0 = Arc(b, d, a)
west0 = Line(b, c)

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
# The following flow condition is not really important because
# the actual data will be taken from the user-defined boundaries.
initial = FlowCondition(p=1000.0, u=0.0, v=0.0, T=348.43)

blk_0 = Block2D(psurf=make_patch(north0, east0, south0, west0),
                fill_condition=initial,
                nni=80, nnj=40,
                bc_list=[UserDefinedBC("udf-vortex-flow.lua"), 
                         ExtrapolateOutBC(),
                         UserDefinedBC("udf-vortex-flow.lua"),
                         UserDefinedBC("udf-vortex-flow.lua")],
                label="Duct")

# Simulation-control information
gdata.flux_calc = ADAPTIVE
gdata.max_time = 20.0e-3
gdata.max_step = 6000
gdata.dt = 1.0e-6
gdata.dt_plot = 5.0e-3

# Some hints to scale and place the SVG layout figure.
sketch.xaxis(0.0, 1.5, 0.5, -0.1)
sketch.yaxis(0.0, 1.5, 0.5, -0.1)
sketch.window(0.0, 0.0, 1.5, 1.5, 0.05, 0.05, 0.17, 0.17)

