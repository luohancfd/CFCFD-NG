## \file ram2.py
## \brief Ramjet design from Birte Haker -- this time with ReactionZones
## \author PJ, 25-Mar-2007, 02-Jun-2007
## \version 22-Apr-2011 (updated for Eilmer3)
##

from math import asin, sin, cos, pi
def radians(degrees): return degrees/180.0*pi
def polar_vector(degrees):
    return Vector(cos(radians(degrees)),sin(radians(degrees)))

gdata.title = "Ramjet in Mach 1.8 Flight -- Reacting Gas."
print gdata.title

# Gas model consists of a reactant and a product gas.
select_gas_model(fname='a-b-gas-model.lua')
set_reaction_scheme('single-step-reaction.lua', reacting_flag=1)

# Define flow conditions
p_inf = 90.0e3  # Pa
T_inf = 282     # degrees K
a_inf = sqrt(1.4 * 287.0 * T_inf)
M_inf = 1.8
u_inf = M_inf * a_inf
print "flight speed=", u_inf, "m/s"
fuel_air_mix = [0.5, 0.5, 0.0]
initial = FlowCondition(p=p_inf/3.0, u=0.0, v=0.0, T=T_inf, 
                        massf=fuel_air_mix)
inflow  = FlowCondition(p=p_inf, u=M_inf*a_inf, v=0.0, T=T_inf, 
                        massf=fuel_air_mix)

# Set up the geometry in the (x,y)-plane be first defining key parameters,
# then defining points and paths using those parameters.
mm = 1.0e-3         # metres per mm
L0 = 70.0 * mm      # distance upstream of inlet
Y0 = 300.0 * mm     # outer edge of flow domain
Y1 = 55.0/2 * mm    # inlet
Y2 = 63.0 * mm      # maximum width of body
L1 = 120.0 * mm     # forebody length
L2 = 300.0 * mm     # combustor length
L3 = 30.0 * mm      # nozzle length
L4 = 100.0 * mm     # trailing body
# Center-body
R2 = 24.0 * mm      # bluff nose
R3 = R2 + 20.0 * mm # block thickness around nose
Y3 = 21.0 * mm      # centrebody diameter

# Derived parameters.
R1 = (L1 * L1 + (Y2 - Y1)**2) / (2.0 * (Y2 - Y1))  # raduis of curvature of forebody
theta1 = asin(L1 / R1)
print "R1=", R1, "theta1=", theta1

# Key points on (x,y)-plane
# Outside of cowl
a = Node(0.0, Y1, label="a")
b = Node(L1, Y2, label="b")
c = Node(L1, Y2-R1, label="c")  # centre of curvature
z0 = Node(-L0, 0.0); z1 = Node(-L0, 2*Y1/3.0)
z2 = Node(-L0, Y1); z3 = Node(-L0, Y0)
z4 = Node(0.0, 0.0); z5 = Node(0.0, 2*Y1/3.0); z6 = Node(0.0, Y0)
z7 = Node(L1, Y0); z8 = Node(L1+L2+L3, Y0); z9 = z8 + Vector(L4,0.0)
b2 = b + Vector(L2+L3,0.0); z10 = b2 + Vector(L4,0.0)
# Inside of cowl.
g = b - Vector(0.0,5.0*mm)
h = 0.33*a + 0.67*g  # part way along the line
f = g + Vector(L2-5.0*mm, 0.0)
e = f + Vector(0.0, -5.0*mm)
# d is defined after q

# Centrebody
i = Node(L1-5.0*mm, 0.0)
j = i - Vector(R2, 0.0)
k = i + R2 * polar_vector(135.0)
l = i + Vector(0.0, R2)
m = Node(L1, R2)
n = Node(L1, Y3)
o = Node(L1+L2, Y3)
p = o + Vector(0.0, 5.0*mm)
q = p + 20.0 * mm * polar_vector(30.0)
r = q + 5.0 * mm * polar_vector(10.0)
s = q + 20.0 * mm * polar_vector(5.0)
t = Node(z10.x, s.y)

# Finish throat of nozzle
d = Node(q.x, e.y)

# Intermediate points through combustor
z11 = i - Vector(R3, 0.0)
z12 = i + R3 * polar_vector(135.0)
z13 = i + Vector(0.0, R3)
z14 = Vector(L1, R3)

def simple_quad_patch(ne, se, sw, nw):
    "Make a quadrilateral from its corners."
    south = Line(sw, se); north = Line(nw, ne)
    west = Line(sw, nw); east = Line(se, ne)
    return make_patch(north, east, south, west)

# Define the blocks, boundary conditions and set the discretisation.
nx0 = 30; ny0 = 20; ny1 = 10; ny2 = 70
nx1 = 50; nx2 = int(nx1*0.6); nx3 = 4*nx1; ny3 = 20; ny4 = 4
nx4 = 40

# Region upstream of inlet.
blk_0 = Block2D(simple_quad_patch(z5, z4, z0, z1), nni=nx0, nnj=ny0,
                fill_condition=initial, label="[0]")
blk_1 = Block2D(simple_quad_patch(a, z5, z1, z2), nni=nx0, nnj=ny1,
                fill_condition=initial, label="[1]")
blk_2 = Block2D(simple_quad_patch(z6, a, z2, z3), nni=nx0, nnj=ny2,
                fill_condition=initial, label="[2]")

# Blocks around the outer cowl.
blk_3 = Block2D(make_patch(Line(z6,z7),Line(b,z7),Arc(a,b,c),Line(a,z6)),
                nni=nx1+nx2, nnj=ny2, fill_condition=initial, label="[3]")
blk_4 = Block2D(simple_quad_patch(z8, b2, b, z7), nni=nx3, nnj=ny2,
                fill_condition=initial, label="[4]")

# Blocks inside entrance.
blk_5 = Block2D(make_patch(Line(z5,z12),Arc(z11,z12,i),Line(z4,z11),Line(z4,z5)),
                nni=nx1, nnj=ny0, fill_condition=initial, label="[5]")
blk_6 = Block2D(simple_quad_patch(h, z12, z5, a), nni=nx1, nnj=ny1,
                fill_condition=initial, label="[6]")
blk_7 = Block2D(make_patch(Arc(z11,z12,i),Line(k,z12),Arc(j,k,i),Line(j,z11)),
                nni=ny0, nnj=ny3, fill_condition=initial, label="[7]")
south8 = Polyline([Arc(k,l,i),Line(l,m)])
south9north8 = Polyline([Arc(z12,z13,i),Line(z13,z14)])
blk_8 = Block2D(make_patch(south9north8,Line(m,z14),south8,Line(k,z12)),
                nni=nx2, nnj=ny3, fill_condition=initial, label="[8]")
blk_9 = Block2D(make_patch(Line(h,g),Line(z14,g),south9north8,Line(z12,h)),
                nni=nx2, nnj=ny1, fill_condition=initial, label="[9]")

# Long part of combustor
blk_10 =  Block2D(simple_quad_patch(p, o, n, m), nni=nx3, nnj=ny4,
                fill_condition=initial, label="[10]")
blk_11 =  Block2D(simple_quad_patch(e, p, m, z14), nni=nx3, nnj=ny3,
                fill_condition=initial, label="[11]")
blk_12 =  Block2D(simple_quad_patch(f, e, z14, g), nni=nx3, nnj=ny1,
                fill_condition=initial, label="[12]")

# Converging-diverging nozzle
north13 = Polyline([Line(e,d),Line(d,b2)])
south13 = Polyline([Line(p,q),Line(q,r),Line(r,s)])
blk_13 = Block2D(make_patch(north13,Line(s,b2),south13,Line(p,e)),
                 nni=nx4, nnj=ny3, fill_condition=initial, label="[13]")

# Blocks after nozzle.
blk_14 = Block2D(simple_quad_patch(z10, t, s, b2), nni=nx3, nnj=ny3,
                 fill_condition=initial, label="[14]")
blk_15 = Block2D(simple_quad_patch(z9, z10, b2, z8), nni=nx3, nnj=ny2,
                 fill_condition=initial, label="[15]")

identify_block_connections()

# Boundary conditions updated.
blk_0.bc_list[WEST] = SupInBC(inflow)
blk_1.bc_list[WEST] = SupInBC(inflow)
blk_2.bc_list[WEST] = SupInBC(inflow)

blk_14.bc_list[EAST] = ExtrapolateOutBC()
blk_15.bc_list[EAST] = ExtrapolateOutBC()

# Do a little more setting of global data.
gdata.axisymmetric_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 5.0e-3  # seconds
print "maximum simulation time=", gdata.max_time, "seconds"
gdata.max_step = 80000
gdata.dt = 1.0e-8
gdata.dt_plot = gdata.max_time / 20.0  # so that we get a few images
gdata.dt_history = 10.0e-5

# We'll let reaction occur in the downstream end of the combustor
# and into the nozzle.
ReactionZone(Vector(L1+L2/4.0,Y3), Vector(L1+L2+L3,Y2))
# Also, we want the flow through the combustor to establish before
# allowing reactions to become active.
gdata.reaction_time_start = 2.5e-3

HistoryLocation(10.0*mm, 2.0*mm)  # just inside the inlet

sketch.xaxis(-0.10, 0.5, 0.1, -0.012)
sketch.yaxis(0.0, 0.25, 0.1, -0.004)
sketch.window(-0.1, 0.0, 0.4, 0.2, 0.0, 0.0, 0.6, 0.25)
