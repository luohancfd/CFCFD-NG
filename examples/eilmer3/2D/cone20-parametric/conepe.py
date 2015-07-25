# conepe.py
# Simple job-specification making use of parametric capabilities.
# PJ, 25-Jul-2015 -- adapted from the classic cone20, extended

# We can set individual attributes of the global data object.
gdata.dimensions = 2
gdata.axisymmetric_flag = 1

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)
Tinf = 1103.0 
a = sqrt(1.4*287.1*Tinf) # sound speed in free stream
M = 1.5
ux = M * a 
print "ux=", ux
inflow  = FlowCondition(p=95.84e3, u=ux, v=0.0, T=Tinf)

# Define cone/flow-domain geometry
theta = 32 # cone half-angle, degrees
L = 0.8 # axial length of cone, metres
rbase =  L * math.tan(math.radians(theta))
x0 = 0.2 # upstream distance to cone tip
H = 3.0 # height of flow domain, metres

gdata.title = "Mach %.1f flow over a %.1f-degree cone." % (M, theta)
print gdata.title

# Set up two quadrilaterals in the (x,y)-plane by first defining
# the corner nodes, then the lines between those corners.
a = Node(0.0, 0.0, label="A")
b = Node(x0, 0.0, label="B")
c = Node(x0+L, rbase, label="C")
d = Node(x0+L, H, label="D")
e = Node(x0, H, label="E")
f = Node(0.0, H, label="F")
ab = Line(a, b); bc = Line(b, c) # lower boundary including cone surface
fe = Line(f, e); ed = Line(e, d) # upper boundary
af = Line(a, f); be = Line(b, e); cd = Line(c, d) # vertical lines

# Define the blocks, with particular discretisation.
dx = 1.0/40
nx0 = int(x0/dx); nx1 = int(L/dx); ny = int(H/dx)
blk_0 = Block2D(make_patch(fe, be, ab, af), nni=nx0, nnj=ny,
                fill_condition=initial, label="BLOCK-0")
blk_1 = Block2D(make_patch(ed, cd, bc, be, "AO"), nni=nx1, nnj=ny,
                fill_condition=initial, label="BLOCK-1",
                hcell_list=[(9,0)], xforce_list=[0,0,1,0])
# Extend the flow domain
xend = x0+2*L
blk_2 = Block2D(CoonsPatch(c,Vector(xend,rbase/2),Vector(xend,H),d),
                nni=nx1, nnj=ny,
                fill_condition=initial, label="BLOCK-2")

# Set boundary conditions.
identify_block_connections()
blk_0.bc_list[WEST] = SupInBC(inflow, label="inflow-boundary")
blk_2.bc_list[EAST] = ExtrapolateOutBC(label="outflow-boundary")

# Do a little more setting of global data.
gdata.max_time = 30.0e-3  # seconds
gdata.max_step = 15000
gdata.dt = 1.0e-6
gdata.dt_plot = 1.5e-3
gdata.dt_history = 10.0e-5

sketch.xaxis(0.0, 2.0, 0.5, -0.05)
sketch.yaxis(0.0, 2.0, 0.5, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)
