# inf_cyl.py
# Rowan J Gollan & PJ
# 10-May-2005, 19-June-2005
#
# This file can be used to simulate the test case reported by:
# Giordano, et al. (1997)
# Vibrationally Relaxing Flow of N2 past an Infinite Cylinder
# Journal of Thermophysics and Heat Transfer, vol 11, no 1, pp 27 - 35
#
# Description:
# Pure N2 flow (chemically inert)
# Mach = 6.5
# T_inf = 300 K
# p_inf = 50 - 500 Pa
#
# Updated: 20-Feb-2007
# This update brings the script inline with mbcns2
#
# Updated: 17-Feb-2012
# This update brings the script up-to-date for use with eilmer3

Rc = 1.0 # cylinder radius

from math import sqrt
from cfpylib.gasdyn.billig import x_from_y

# Setup simulation
gdata.title = "Inifinite cylinder in N2 flow (M=6.5)"

# Inflow conditions
M_inf = 6.5
T_inf = 300.0
p_inf = 50.0

# Set the gas model and
# use this to compute some other flow properties

select_gas_model(fname='gas-model.lua')

gmodel = get_gas_model_ptr()
gd = Gas_data(gmodel)
gd.p = p_inf
gd.T[0] = T_inf
gd.T[1] = T_inf
gd.massf[0] = 1.0
gmodel.eval_thermo_state_pT(gd)

u_inf = M_inf * gd.a

# Thermal energy exchange model
set_energy_exchange_update("N2-TV.lua")

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=[T_inf, T_inf], massf=[1.0,])
initial = FlowCondition(p=p_inf/3.0, u=0.0, v=0.0, T=[T_inf, T_inf], massf=[1.0,])

# Build the geometry from the bottom-up, starting with nodes...
# Scale it with the cylinder radius.

# In preparation for defining nodes, generate a few sample points
# along the expected shock position
# (which is estimated via Billig's correlation).
xys = []
for y in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]:
    x = x_from_y(y, M_inf, theta=0.0, axi=0, R_nose=Rc)
    xys.append((x,y)) # a new coordinate pair
    print "x=", x, "y=", y
    
a = Node(-Rc, 0.0, label="a")
b = Node( 0.0, Rc, label="b")
c = Node( 0.0, 0.0, label="c")
d = [] # will use a list to keep the nodes for the shock boundary
for x, y in xys:
    # the outer boundary should be a little further than the shock itself
    d.append(Node(-1.1*x, 1.1*y, label="d"))

# ...then lines, arcs, etc, that will make up the block boundaries.
axis = Line(d[0], a) # first-point of shock to nose of cylinder
cylinder = Arc(a, b, c)
shock = Spline(d)
top = Line(d[-1], b)  # top-point of shock to top of cylinder

# Specify the boundary discretization and conditions...
nnx = 30
nny = 3 * nnx

# ...and finally, assemble the block from its boundary faces.
block_0 = Block2D(make_patch(top, cylinder, axis, shock),
                  nni=nnx, nnj=nny,
                  bc_list=[ExtrapolateOutBC(), SlipWallBC(), SlipWallBC(), SupInBC(inflow)],
                  fill_condition=initial,
                  hcell_list=[(nnx, 1)])

# simulation control
gdata.flux_calc = AUSMDV
gdata.max_time = 10.0e-3 # should be large enough to allow steady flow
gdata.max_step = 30000
gdata.cfl = 0.25
gdata.stringent_cfl = 1
gdata.dt = 1.0e-8
gdata.dt_history = 1.0e-5
gdata.dt_plot = 1.0e-3

# Metapost setting to get a decent arrangement on A4 paper.
sketch.scales(0.05, 0.05)
sketch.origin(0.15, 0.05)
sketch.xaxis(-2.0, 0.0, 0.5, -0.1)
sketch.yaxis(0.0, 3.0, 0.5, 0.0)

