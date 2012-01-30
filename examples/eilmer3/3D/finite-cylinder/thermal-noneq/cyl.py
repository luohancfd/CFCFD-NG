# \file cyl.py
#
# This geometry is a set of three blocks describing a quarter-cylinder
# of finite length in supersonic flow.
#
# PJ, 20-Jun-2005, 04-Dec-2005 increase number of blocks along cylinder axis
#     06-Feb-2006 new geometry objects
#     19-Aug-2009 Eilmer3 port
#     23-Jan-2010 SuperBlock3D and use of MPI code cor comparison
# RJG, 02-Apr-2007 new reacting gas spec.
# DFP, 08-Dec-2011 port to thermal nonequilibrium

gdata.dimensions = 3
D = 15.0e-3  # Diameter of cylinder, metres
L = 2.0 * D  # (axial) length of cylinder

# Gas model used in the simulation.
select_gas_model(model='two temperature gas', species=['N2','N','N2_plus','N_plus','e_minus'])
set_reaction_scheme("nitrogen-5sp-6r.lua",reacting_flag=1)
set_energy_exchange_update("TV-TE_exchange.lua")
mf = {'N2':1.0}

# Free-stream properties
T_inf = 3000.0 # degrees K
p_inf = 2000.0 # Pa
u_inf = 10000.0 # m/s
gdata.title = "Cylinder L/D=%g in N2 at u=%g m/s." % (L/D, u_inf)
print "title=", gdata.title

# Flow conditions for fill and boundary conditions.
inflowCond = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=mf)
initialCond = FlowCondition(p=p_inf/3.0, u=0.0, v=0.0, T=300.0, massf=mf)

# Geometry is built from the bottom up.
Rc = D/2.0  # cylinder radius
    
# Define a few key nodes.
a = Node(-Rc, 0.0, 0.0, label="a")  # stagnation point on the cylinder
b = Node( 0.0, Rc, 0.0, label="b")  # top of cylinder
c = Node( 0.0, 0.0, 0.0, label="c") # centre of curvature

# In order to have a grid that fits reasonably close the the shock,
# use Billig's shock shape correlation to generate
# a few sample points along the expected shock position.
from math import sqrt
from cfpylib.gasdyn.billig import x_from_y
# ideal N2 properties used for shock shape estimate
R_N2 = 296.8
gamma_N2 = 1.4
a_inf = sqrt(gamma_N2 * R_N2 * T_inf)
M_inf = u_inf / a_inf
print "M_inf=", M_inf
xys = []
for y in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]:
    x = x_from_y(y*Rc, M_inf, theta=0.0, axi=0, R_nose=Rc)
    xys.append((x,y*Rc)) # a new coordinate pair
    print "x=", x, "y=", y

# Scale the Billig distances, depending on the expected behaviour
# relative to the gamma=1.4 ideal gas.
if gdata.reacting_flag == 1:
    b_scale = 0.87  # for finite-rate chemistry
else:
    b_scale = 1.1   # for ideal (frozen-chemistry) gas
d = [] # will use a list to keep the nodes for the shock boundary
for x, y in xys:
    # the outer boundary should be a little further than the shock itself
    d.append(Node(-b_scale*x, b_scale*y, 0.0, label="d"))
print "front of grid: d[0]=", d[0]

# Extent of the cylinder in the z-direction to end face.
c2 = c.clone(); c2.translate(0.0, 0.0, L/2.0)
e = d[0].clone().translate(0.0, 0.0, L/2.0)
f = a.clone().translate(0.0, 0.0, L/2.0)
g = Node(-Rc/2.0, 0.0, L/2.0)
h = Node(0.0, Rc/2.0, L/2.0)
i = Node(0.0, Rc, L/2.0)
# the domain is extended beyond the end of the cylinder
j = e.clone().translate(0.0, 0.0, Rc)
k = f.clone().translate(0.0, 0.0, Rc)

# ...then lines, arcs, etc, that will make up the domain-end face.
xaxis = Line(d[0], a) # first-point of shock to nose of cylinder
cylinder = Arc(a, b, c)
shock = Spline(d)
outlet = Line(d[-1], b)  # top-point of shock to top of cylinder
domain_end_face = CoonsPatch(xaxis, outlet, shock, cylinder)

# ...lines along which we shall extrude the domain-end face
yaxis0 = Line(d[0], e)
yaxis1 = Line(e, j)

# End-face of cylinder
xaxis = Line(f, g)
cylinder = Arc(f, i, c2)
inner = Arc(g, h, c2)
outlet = Line(i, h)
cyl_end_face = CoonsPatch(xaxis, outlet, cylinder, inner)
yaxis2 = Line(f, k)

# Third, set up the blocks from the geometric and flow elements.
nr = 20            # radial discretization
nc = int(1.5 * nr) # circumferential discretization
na = int(L/D * nc) # axial discretization along the cylinder
na1 = nc           # axial discretization off the end of the cylinder
nr2 = int(nr/2)    # radial discretization toward the cylinder axis

# The volume constructor extrudes the end-face along the axis in the k-direction.
# We want to divide the over-cylinder block up to make reasonable use of the
# cluster computer.
blk0 = SuperBlock3D(label="over-cylinder", nni=nr, nnj=nc, nnk=na, nbk=int(L/D),
                    parametric_volume=WireFrameVolume(domain_end_face,yaxis0,"k"),
                    fill_condition=initialCond)
for blk in blk0.blks[0][0]:
    # We work along the line of blocks in the k-direction
    blk.set_BC("WEST", "SUP_IN", inflow_condition=inflowCond)
    blk.set_BC("NORTH", "SUP_OUT")

blk1 = Block3D(label="outside-cylinder", nni=nr, nnj=nc, nnk=na1,
               parametric_volume=WireFrameVolume(domain_end_face,yaxis1,"k"),
               fill_condition=initialCond)
blk1.set_BC("WEST", "SUP_IN", inflow_condition=inflowCond)
blk1.set_BC("NORTH", "SUP_OUT")

blk2 = Block3D(label="beside-cylinder", nni=nr2, nnj=nc, nnk=na1,
               parametric_volume=WireFrameVolume(cyl_end_face,yaxis2,"k"),
               fill_condition=initialCond)
blk2.set_BC("EAST", "SUP_OUT")
blk2.set_BC("NORTH", "SUP_OUT")

identify_block_connections()

# Finally, Other simulation control parameters. ----------
gdata.viscous_flag = 0
gdata.flux_calc = ADAPTIVE
gdata.interpolation_type = "pT"
gdata.t_order = 1
gdata.x_order = 2
gdata.max_time = Rc/u_inf * 30
gdata.max_step = 40000
gdata.dt = 1.0e-10
gdata.cfl = 0.5
gdata.dt_history = 1.0e-5
gdata.dt_plot = gdata.max_time/10
gdata.print_count = 1

print "Total number of blocks=", len(blk0.blks)+2

