## \file   hayabusa.py
## \brief  Simulating the subscale Hayabusa model in 9km/s air.
## \author DFP
## \date   23-Nov-2009
## \date   09-May-2011: Simplified for inclusion in the examples by PJ
## The intent is to display the use of Dan's two-temperature gas model.
## The geometry description has been reduced to including just the 
## windward surface of the aeroshell.

gdata.title = "JAXA Hayabusa sample return capsule (subscale)."

species = select_gas_model(model='two temperature gas', 
                           species=['N2', 'N2_plus', 'NO', 'NO_plus', 
                                    'O2', 'O2_plus', 'N', 'N_plus',
                                    'O', 'O_plus', 'e_minus'])
set_reaction_update("Park93-s03-AIC-EIIC.lua")
set_energy_exchange_update("TV-TE.lua")
gm = get_gas_model_ptr()
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()

# Define free-stream flow conditions relevant to Mary D'Souza's experiment.
T_wall = [300.0, 300.0]
rho_inf = 1.733230e-03
T_inf = [1.069042e+03, 1.202479e+03]
u_inf = 9.121299e+03
massf_inf = [7.455342e-01, 5.161281e-17, 2.777519e-04, 9.223930e-05,
             7.806664e-03, 1.176211e-09, 2.156932e-02, 3.532030e-13,
             2.247197e-01, 6.071886e-08, 5.812489e-21,]
# Use the gas model to get pressure, Mach number and total mass-flux.
Q = Gas_data(gm)
for itm in range(ntm):
    Q.T[itm] = T_inf[itm]
Q.rho = rho_inf
for isp,massf in enumerate(massf_inf):
    Q.massf[isp] = massf
gm.eval_thermo_state_rhoT(Q)
global M_inf
M_inf = u_inf / Q.a
p_inf = Q.p
print "p_inf = ", p_inf
print "Freestream velocity corresponding to Mach %0.2f is %0.1f" % ( M_inf, u_inf )
print "total enthalpy = ", gm.total_enthalpy(Q) + 0.5 * u_inf**2
# Print molar fractions for info.
M = vectord(nsp)
molef = vectord(nsp)
for isp in range(nsp):
    M[isp] = gm.molecular_weight(isp)
convert_massf2molef(Q.massf,M,molef)
for isp in range(nsp):
    print "%10s = %12.6e" % (species[isp], molef[isp])

inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=massf_inf)
initial = FlowCondition(p=p_inf/10.0, u=0.0, v=0.0, T=T_inf, massf=massf_inf)

# Geometry construction.
# Setup a quarter sphere with the specified radius.
from math import cos, sin, tan, sqrt, pi
scale = 1.0 / 10.0        # This model may be subscale.
Rn = 200.0e-3 * scale     # nose radius
D = 400.0e-3 * scale      # diameter of heat shield
theta = 45.0 * pi / 180.0 # angle of skirt

origin = Node(0.0, 0.0)
a = Node(-Rn, 0.0, label='a')  # nose of aeroshell
b = Node(-Rn*cos(theta), Rn*sin(theta), label='b')  # transition to straight
c = Node( b.x + (D/2.0-b.y)/tan(theta), D/2.0, label='c')  # outer edge of aeroshell
d = Node( 0.0, c.y - abs(c.x), label='d')

# Generate inflow boundary points by scaling Billig's shock shape.
from cfpylib.gasdyn.billig import x_from_y, y_from_x
bx_scale = by_scale = 0.95
np = 20
y_top = by_scale * y_from_x(0.0, M_inf, theta=0.0, axi=1, R_nose=Rn)
dy = y_top / (np - 1)
inflow_points = []
for iy in range(np):
    y = dy * iy
    x = - bx_scale * x_from_y(y/by_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
    inflow_points.append(Vector(x,y))
inflow_spline = Spline(inflow_points)

# The inflow boundary will be a spline through the Billig points
# but truncated by the line rising from the end of 
# the straight-edge of the aeroshell.
from cfpylib.nm.zero_solvers import bisection
def intersection_func(y):
    dc_line = Line(d, Vector(d.x+(c.x-d.x)*2.0, d.y+(c.y-d.y)*2.0))
    t = (y-d.y)/((c.y-d.y)*2.0)
    return inflow_spline.eval_from_y(y).x - dc_line.eval(t).x
y_int = bisection(intersection_func, by=0.0, uy=D)
x_int = inflow_spline.eval_from_y(y_int).x
# Truncate the inflow spline.
west0_nodes = [point for point in inflow_points if point.y < y_int]
west0_nodes.append(Node(x_int, y_int))

# Bounding paths are now assembled from the nodes.
north0 = Line(west0_nodes[-1],c)                 # outflow boundary
east0  = Polyline([Arc(a,b,origin), Line(b,c)])  # aeroshell surface itself
south0 = Line(west0_nodes[0],a)                  # axis of symmetry
west0  = Spline(west0_nodes)       # Billig shock shape for inflow boundary

# Discretize the region and attach the boundary conditions.
nni = nnj = 30
bc_list = [ExtrapolateOutBC(), FixedTBC(T_wall[0]), SlipWallBC(), SupInBC(inflow)]
blk = SuperBlock2D(psurf=make_patch(north0, east0, south0, west0),
                   fill_condition=initial,
                   nni=nni, nnj=nnj, nbi=3, nbj=3, cf_list=[None,]*4,
                   bc_list=bc_list, wc_bc_list=[NonCatalyticWBC()]*4,
                   label="BLK", hcell_list=[(nni,1)])

# Do a little more setting of global data to control the simulation.
gdata.axisymmetric_flag = 1
gdata.reaction_time_start = Rn * 1 / u_inf
gdata.viscous_flag = 0
gdata.flux_calc = ADAPTIVE
gdata.max_time = Rn * 10 / u_inf    # 10 body lengths
gdata.max_step = 230000
gdata.dt = 1.0e-11
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1.0 / u_inf    # a solution every body length
gdata.cfl = 0.5
gdata.cfl_count = 1
gdata.print_count = 20

# Adjust the SVG rendering so that if may be used for documentation.
sketch.xaxis(-1.25*Rn, 0.0, 0.25*Rn, -0.1*Rn)
sketch.yaxis(0.0, 1.25*Rn, 0.25*Rn, 0.0)
sketch.window(-1.25*Rn, 0.0, 0.0, 1.25*Rn, 0.05, 0.05, 0.17, 0.17)
