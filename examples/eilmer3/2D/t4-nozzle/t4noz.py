# t4noz.py
# Supersonic part of the T4 nozzle.
# Which particular nozzle is determined by the contour file.
#
# PJ, 14-Mar-2011 adapted from Rainer Kirchhartz 2009 simulation files.
#---------------------------------------------------------------------
gdata.title = "Flow through the T4 nozzle."

# Shock tube conditions define flow condition at the nozzle throat.
# This throat condition is used as the input to the simulation of
# the nozzle expansion.
from estcj import reflected_shock_tube_calculation as rstc
gasName = 'air'
# Katsu's T4 shot 9378
T1 = 300.0   # test-gas fill temperature, K
p1 = 160.0e3 # test-gas fill pressure, Pa
Vs = 2707.0  # incident shock speed, m/s
pe = 38.0e6  # equilibrium pressure in nozzle supply region, Pa
estc_result = rstc(gasName, p1, T1, Vs, pe, area_ratio=10.0, task='stn')
throat = estc_result['state6']
throat_V = estc_result['V6']
print "Flow condition at nozzle throat:"
throat.write_state(sys.stdout)

# Estimate turbulence quantities for free stream
# by specifying the intensity as 0.05 and estimating the
# turbulence viscosity as 100 times the laminar viscosity.
throat_tke = 1.5 * (throat_V * 0.05)**2
throat_mu_t = 100.0 * throat.mu
throat_omega = throat.rho * throat_tke / throat_mu_t
print "Inflow turbulence: tke=", throat_tke, "omega=", throat_omega

select_gas_model(fname="cea-lut-air.lua.gz" )
inflow = FlowCondition(p=throat.p, u=throat_V, v=0.0, T=throat.T, 
		       massf=[1.,], tke=throat_tke, omega=throat_omega)
initial = FlowCondition(p=130.0, u=0.0, v=0.0, T=300.0, 
			massf=[1.,], tke=0.0, omega=1.0)

# Read in contour file for the particular nozzle.
x_c = []; y_c = []
fp = open("contour-t4-m4.data", 'r')
for line in fp.readlines():
    items = line.strip().split()
    if len(items) == 0: continue  # skip blank lines
    if items[0] == '#': continue  # skip header lines
    x_c.append(float(items[0]))
    y_c.append(float(items[1]))
area_ratio = (y_c[-1]/y_c[0])**2
print "nozzle area ratio=", area_ratio

# Define the nodes, which make up the nozzle contour in the simulation. 
dlist = [Vector(x_c[i],y_c[i]) for i in range(len(x_c))]
# Define nodes for rest of nozzle geometry.
# The throat is a bit of straight section before onset of expansion.
throat_centre = Node(-y_c[0],0.0) 
throat_wall = Node(-y_c[0],y_c[0])
# Start expansion at x=0.
divergence_centre = Node(0.0,0.0)
divergence_wall = Node(0.0,y_c[0])
end_of_cone = Node(x_c[1],y_c[1])
nozzle_end_centre = Node(x_c[-1],0.0)
nozzle_end_wall = Node(x_c[-1],y_c[-1])

# Define the lines making up the patch
n_south = Polyline([Line(throat_centre, divergence_centre), 
		    Line(divergence_centre, nozzle_end_centre)])
n_north = Polyline([Line(throat_wall, divergence_wall), 
		    Line(divergence_wall, end_of_cone), 
		    Spline(dlist[1:len(dlist)])])
n_west = Line(throat_centre, throat_wall)
n_east = Line(nozzle_end_centre, nozzle_end_wall)
expansion_region = make_patch(n_north, n_east, n_south, n_west)

# Boundary conditions
Twall = 300
bc_list = [FixedTBC(Twall), ExtrapolateOutBC(), SlipWallBC(), SupInBC(inflow)]

# Define the clustering and grid resolution.
x_clust = RobertsClusterFunction(1, 1, 1.05)
y_clust = RobertsClusterFunction(0, 1, 1.05)
cf_list = [x_clust, y_clust, x_clust, y_clust]
nozzle_x = 300
nozzle_y = 50

# The nozzle is defined as a linear sequence of blocks
# to allow space-marching in the simulation.
nbi = 20
nozzle_blk = SuperBlock2D(expansion_region, nni=300, nnj=50, nbi=nbi, nbj=1, 
			  bc_list=bc_list, cf_list = cf_list,
			  fill_condition=initial, label="nozzle")

# Attributes for the space-marching simulation.
gdata.dimensions = 2
gdata.sequence_blocks = 1
gdata.axisymmetric_flag = 1
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 1.50e-3  # seconds
gdata.max_step = 800000
gdata.dt = 0.2e-9
gdata.cfl = 0.4
gdata.cfl_count = 5
gdata.dt_plot = 0.04e-3
gdata.dt_history = 10.0e-6
gdata.turbulence_flag = 0
gdata.turbulence_model = 'baldwin_lomax'

sketch.xaxis(0.0, 1.2, 0.1, 0.)
sketch.yaxis(0.0, 0.4, 0.1, 0.)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.015, 0.015, 0.15, 0.15)
