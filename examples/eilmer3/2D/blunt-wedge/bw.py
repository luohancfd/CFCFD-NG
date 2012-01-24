# bw.py
# MECH4470/CFD Exercise: Hypersonic flow over a blunt wedge.
# PJ, 07-Dec-2006
#     31-Jan-2010 ported to Eilmer3

from math import sqrt, sin, cos, tan, pi

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Free stream
g_gas = 1.4                     # Ideal Air
R_gas = 287.0 
M_inf = 5.0                     # Specified Mach number
p_inf = 100.0e3                 # Select a static pressure
T_inf = 100.0                   # and a temperature
a_inf = sqrt(T_inf * R_gas * g_gas) # determine sound speed
u_inf = M_inf * a_inf                    # and velocity
# Also, handy to know dynamic pressure for nondimensionalization
# of the pressures and drag forces.
q_inf = 0.5 * g_gas * p_inf * M_inf * M_inf
print "Free-stream velocity, u_inf=", u_inf
print "    static  pressure, p_inf=", p_inf
print "    dynamic pressure, q_inf=", q_inf
free_stream = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf)
# For transient simulation, we start with a low pressure.
initial = FlowCondition(p=1000.0, u=0.0, v=0.0, T=100.0)

# Geometry
Rn = 10.0e-3                    # radius of cylindrical nose
xEnd = 8.0 * Rn                 # downstream extent of wedge
alpha = 10.0 / 180.0 * pi       # angle of wedge wrt free stream
delta = 10.0e-3                 # offset for inflow boundary

# First, specify surface of cylinder and wedge
a = Node(0.0, 0.0, label='a') # Centre of curvature for nose
b = Node(-Rn, 0.0, label='b')
c = Node(-Rn*sin(alpha), Rn*cos(alpha), label='c')
bc = Arc(b, c, a)
# Down-stream end of wedge
d = Node(xEnd, c.y+(xEnd-c.x)*tan(alpha), label='d')
print "height at end of plate yd=", d.y
cd =Line(c, d)

# Outer-edge of flow domain has to contain the shock layer
# Allow sufficient for shock stand-off at the stagnation line.
R2 = Rn + delta
e = Node(-R2, 0.0, label='e')
# The shock angle for a 10 degree ramp with sharp leading edge
# is 20 degrees (read from NACA 1135, chart 2),
# however, the blunt nose displaces the shock a long way out
# so we allow some more space.
# We need to set the boundary high enough to avoid the shock
R3 = Rn + 2.0 * delta
f = Node(-R3*sin(alpha), R3*cos(alpha), label='f')
# Now, put in intermediate control points so that we can use
# cubic Bezier curve for the inflow boundary around the nose
# and a straight line downstream of point f.
e1 = Node(e.x, delta, label='e1')
alpha2 = 40.0 / 180.0 * pi
f1 = Node(f.x-delta*cos(alpha2), f.y-delta*sin(alpha2), label='f1')
ef = Bezier([e, e1, f1, f])
g = Node(xEnd, f.y+(xEnd-f.x)*tan(alpha2), label='g')
fg = Line(f,g)

# Define straight-line segments between surface and outer boundary.
eb = Line(e, b); fc = Line(f, c); dg = Line(d, g)

# Define the blocks using the path segments.
# Note that the EAST face of region0 wraps around the nose and
# that the NORTH face of region0 is adjacent to the WEST face
# of region1.
region0 = make_patch(fc, bc, eb, ef)
cf = fc.copy(); cf.reverse() # common boundary but opposite sense
region1 = make_patch(fg, dg, cd, cf)
cluster0 = RobertsClusterFunction(0, 1, 1.2)
cluster1 = RobertsClusterFunction(1, 0, 1.2)
nni0 = 40
nnj0 = 40
nni1 = 100
blk_0 = Block2D(region0, nni=nni0, nnj=nnj0,
                cf_list=[cluster0,None,cluster0,None],
                fill_condition=initial,
                xforce_list=[0,1,0,0])
blk_1 = Block2D(region1, nni=nni1, nnj=nnj0,
                cf_list=[None,cluster1,None,cluster1],
                fill_condition=initial,
                xforce_list=[0,0,1,0])
identify_block_connections()
blk_0.bc_list[WEST] = SupInBC(free_stream)
blk_1.bc_list[NORTH] = SupInBC(free_stream)
blk_1.bc_list[EAST] = ExtrapolateOutBC()


# We can set individual attributes of the global data object.
job_title = "Blunt Wedge Rn=" + str(Rn)
job_title += (" q_inf=%12.3e" % q_inf) + (" d.y=%10.5f" % d.y)
print job_title
gdata.title = job_title
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 40.0 * Rn / u_inf
print "Final time=", gdata.max_time
gdata.max_step = 5000
gdata.dt = 1.0e-8
gdata.dt_plot = gdata.max_time
gdata.dt_history = gdata.max_time / 100.0
HistoryLocation(b.x-0.001, b.y) # just in front of the stagnation point

sketch.xaxis(-0.020, 0.080, 0.020, -0.004)
sketch.yaxis(0.0, 0.100, 0.020, -0.004)
sketch.window(-0.02, 0.0, 0.08, 0.10, 0.05, 0.05, 0.17, 0.17)




