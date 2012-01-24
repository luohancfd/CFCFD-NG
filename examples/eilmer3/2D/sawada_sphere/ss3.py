# file: ss3.py
#
# Sphere in equilibrium air modelling Case 3 from
#     K. Sawada & E. Dendou (2001)
#     Validation of hypersonic chemical equilibrium flow calculations
#     using ballistic-range data.
#     Shock Waves (2001) Vol. 11, pp 43--51
#
# Experimental shock stand-off distance is 2.59mm
# Sawada & Dendou CFD value:               2.56mm
#
# This script derived from rbody, 22-Jan-2004.
# and the Python version: ss3.py, 04-Apr-2005, 10-Aug-2006, 27-Nov-2006
# PJ
#
# The grid is a bit wasteful because the shock lies close to
# the body for equilibrium air, however, this grid layout 
# (as used in rbody) allows us to play with perfect-gas models
# without hitting the inflow boundary with the shock.
#
# Updated: 12-Nov-2008 by RJG for use in Elmer3


# The following JOB name is used to build file names at the end.
JOB = "ss3"

# Radius of body
R = 31.8e-3                      # m
T_body = 296.0                   # surface T, not relevant for inviscid flow
body_type = "sphere"             # choose between "cylinder" and "sphere"

# Free-stream flow definition
p_inf = 20.0e3                   # Pa
T_inf = 296.0                    # degrees K
u_inf = 4.68e3                   # flow speed, m/s

# For equilibrium chemistry, use the look-up-table (which has
# been previously created).
print "About to select gas model."
select_gas_model(fname='cea-lut-air.lua.gz')
print "Gas model selection: done."

# Define flow conditions
inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf)
initial = FlowCondition(p=0.3*p_inf, u=0.0, v=0.0, T=T_inf)

# Job-control information
do_viscous = 0                   # flag for viscous/inviscid calc
nn = 60                          # grid resolution, both ix and iy
t_final = 10.0 * R / u_inf       # allow time to settle at nose
t_plot = t_final / 1.0           # plot only once
TitleText = "Blunt Body " + JOB + ": R=" + str(R) + ", gas='equilibrium air'" + \
            ", p=" + str(p_inf) + ", v=" + str(u_inf) + ", T=" + str(T_inf) + \
            ", viscous=" + str(do_viscous)
gdata.title = TitleText
gdata.case_id = 0
if do_viscous:
    gdata.viscous_flag = 1
    gdata.viscous_delay = t_plot
if body_type == "sphere":
    gdata.axisymmetric_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = t_final
gdata.max_step = 400000
gdata.dt = 1.0e-8
gdata.cfl = 0.40
gdata.dt_plot = t_plot
gdata.dt_history = 1.0e-6
 
# Begin geometry details...
# Note that mbcns_prep.py has already imported the math module.
deg2rad = math.pi / 180.0
alpha1 = 135.0 * deg2rad
alpha2 = 50.8 * deg2rad
# The node coordinates are scaled with the body radius.
# The labels are not required but make the MetaPost plot
# look a little like the plot produced by scriptit.tcl.
a = Node(0.0, 0.0, label="a")
b = Node(-1.0 * R, 0.0, label="b")
c = Node(math.cos(alpha1) * R, math.sin(alpha1) * R, label="c")
d = Node(0.0, R, label="d")
e = Node(math.cos(alpha2) * R, math.sin(alpha2) * R, label="e")
f = Node(1.4 * R, 1.5 * R, label="f")
g = Node(1.5 * R, 2.5 * R, label="g")
h = Node(1.5 * R, 3.5 * R, label="h")
i = Node(-1.5 * R, 0.0, label="i")
j = Node(-1.5 * R, 1.5 * R, label="j")
k = Node(-1.0 * R, 2.8 * R, label="k")

east0  = Polyline([Arc(b, c, a), Arc(c, d, a), Arc(d, e, a)])
north0 = Bezier([e, f, g, h,]); north0.reverse()
south0 = Line(i, b)
west0  = Bezier([i, j, k, h,])

print "ss3: block to be defined."
cluster_functions = [RobertsClusterFunction(0, 1, 1.2),
                     RobertsClusterFunction(1, 0, 1.1),
                     RobertsClusterFunction(0, 1, 1.2),
                     RobertsClusterFunction(1, 0, 1.1)]
boundary_conditions = [ExtrapolateOutBC(), FixedTBC(T_body),
                       SlipWallBC(), SupInBC(inflow)]

blk_0 = Block2D(psurf=make_patch(north0, east0, south0, west0),
                fill_condition=initial,
                nni=nn, nnj=nn,
                cf_list=cluster_functions,
                bc_list=boundary_conditions,
                label="BLOCK-0", hcell_list=[(nn,1)])

# Some hints to scale and place the sketch.
# If you change the radius, you'll probably have to adjust the axes.
sketch.xaxis( -0.060, 0.050, 0.020, -0.010)
sketch.yaxis(  0.0,   0.110, 0.020, 0.0)
sketch.window(-1.5*R, 0.0, 1.5*R, 3.0*R, 0.05, 0.05, 0.15, 0.15)
