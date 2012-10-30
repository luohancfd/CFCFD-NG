## \file hayabusa.py
## \brief Simulating the JAXA Hayabusa sample return capsule
## \author DFP, 30-Oct-2012
## 
## Part1: Inviscid solution on a coarse grid

from cfpylib.gasdyn.billig import x_from_y, y_from_x
from cfpylib.nm.zero_solvers import bisection
from math import cos, sin, tan, sqrt, pi
from cfpylib.flow.shock_layer_surface import ShockLayerSurface

job_title = "JAXA Hayabusa sample return capsule."
print job_title

gdata.title = job_title
gdata.axisymmetric_flag = 1

#
# 1. Setup the gas model
#
species = select_gas_model(model='two temperature gas', species=['N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus'])
set_reaction_update("Park93-s03-AIC-EIIC.lua")
set_energy_exchange_update("TV-TE.lua")
gm = get_gas_model_ptr()
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()

#
# 2. Define flow conditions
#
rho_inf = 1.645e-4
T_inf = 233.25
u_inf = 11.6e3
massf_inf = [ 0.0 ] * gm.get_number_of_species()
massf_inf[species.index('N2')] = 0.767
massf_inf[species.index('O2')] = 0.233

# do some calculations to get pressure, Mach number and total mass-flux
Q = Gas_data(gm)
for itm in range(ntm):
    Q.T[itm] = T_inf
Q.rho = rho_inf
for isp,massf in enumerate(massf_inf):
    Q.massf[isp] = massf
gm.eval_thermo_state_rhoT(Q)
M_inf = u_inf / Q.a
p_inf = Q.p
print "p_inf = %0.1f, M_inf = %0.1f" % ( p_inf, M_inf )
inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=[T_inf]*ntm, massf=massf_inf)
initial = FlowCondition(p=p_inf/10.0, u=0.0, v=0.0, T=[T_inf]*ntm, massf=massf_inf)

#
# 3. Define the geometry
#
Rn = 200.0e-3
D = 400.0e-3
theta = 45.0 * pi / 180.0
L = 200.0e-3
global bx_scale, by_scale
bx_scale = 0.95
by_scale = 0.95

# origin
o = Node(0.0,0.0)

# surface nodes
global c,d
a = Node(-Rn,0.0, label='a')
b = Node(-Rn*cos(theta),Rn*sin(theta), label='b')
c = Node( b.x + ( D/2.0-b.y)/tan(theta), D/2.0, label='c')
d = Node( 0.0, c.y - abs(c.x), label='d')

# inflow boundary nodes
x_limit = c.x
inflow_nodes = []
np = 32
y_top = by_scale * y_from_x(-x_limit/bx_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
dy = y_top / ( np - 1 )
for iy in range(np):
    y = dy * iy
    x = - bx_scale * x_from_y(y/by_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
    inflow_nodes.append( Node(x,y) )

# find intersection of surface normal with inflow boundary at the top most point
global inflow_spline
inflow_spline = Spline(inflow_nodes)
def zero_func(y):
    dc_line = Line( d, Node( d.x+(c.x-d.x)*2.0, d.y+(c.y-d.y)*2.0 ) )
    t = (y-d.y)/((c.y-d.y)*2.0)
    return inflow_spline.eval_from_y(y).x - dc_line.eval(t).x
y_int = bisection( zero_func, by=0.0, uy=D )
x_int = inflow_spline.eval_from_y(y_int).x

# split the inflow spline
west_nodes = []
for node in inflow_nodes:
    if node.y < y_int: west_nodes.append(node)
west_nodes.append( Node( x_int, y_int ) )

# curves - block0
west0 = Spline(west_nodes)
east0 = Polyline( [Arc(a,b,o),Line(b,c)] )
#
# 4. Define the blocks, boundary conditions and set the discretisation
#
nnx = 30; nny=30
nbx = 2; nby = 2

blk_0 = SuperBlock2D(psurf=ShockLayerSurface(east0, west0),
		     fill_condition=initial,
		     nni=nnx, nnj=nny,
		     nbi=nbx, nbj=nby,
		     bc_list=[ExtrapolateOutBC(), SlipWallBC(), SlipWallBC(), SupInBC(inflow)],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-0", hcell_list=[(nnx,1)])
 
#
# 5. Simulation control parameters 
#
gdata.viscous_flag = 0
gdata.flux_calc = ADAPTIVE
gdata.max_time = Rn * 5 / u_inf    # 5 body lengths
gdata.max_step = 230000
gdata.dt = 1.0e-8
gdata.reaction_time_start = Rn * 0.1 /u_inf
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1 / u_inf    # 5 solutions
gdata.cfl = 0.5
gdata.cfl_count = 1
gdata.print_count = 20

#
# 6. svg sketch parameters
#
sketch.scales(0.3/Rn, 0.3/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-1.5*Rn, 0.0, 0.25*Rn, -0.1*Rn)
sketch.yaxis(0.0, 1.5*Rn, 0.25*Rn, 0.0)
