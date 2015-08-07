# t4-m4-nozzle.py
# Runs an Eilmer simulation of the T4 M4 nozzle, without the use of 
# NENZFr. This script is built from NENZFr-generated files and allows 
# users to run nozzle simulations that require functionalities outside 
# the standard ones provided by NENZFr. The flow in the T4 M4 nozzle
# is assumed to be laminar and in chemical equilibrium. The inflow 
# conditions are equivalent to that at Mach 7 flight at q of 50 kPa. 
#
# Wilson Chan, 06-Aug-2015
#---------------------------------------------------------------------
from estcj import reflected_shock_tube_calculation as rstc
import cfpylib.gasdyn.cea2_gas as cea2

gdata.title = "Laminar equilibrium flow through the T4 Mach 4 nozzle."

# Run estcj to get flow conditions at the nozzle throat.
estc_result = rstc(cea2, "air", 131e3, 300.0, 1644.0, 8.32e6, 
                   None, area_ratio=27.0, task='stn')
throat = estc_result['state6']
throat_V = estc_result['V6']
print "Flow velocity at nozzle throat = {0} m/s.".format(throat_V)
print "Flow condition at nozzle throat:"
throat.write_state(sys.stdout)  

# Set inflow conditions at throat.
select_gas_model(fname="cea-lut-air.lua.gz" )
inflow = FlowCondition(p=throat.p, u=throat_V, v=0.0, T=throat.T, 
                       massf=[1.0,], label="inflow")

# Read in coordinates that define the contour of the nozzle. For the T4 M4
# nozzle, all coordinate pairs except the first are Bezier control points.
fp = open("Bezier-control-pts-t4-m4.data", 'r')
x_c = []; y_c = []
for line in fp.readlines():
    items = line.strip().split()
    if len(items) == 0: continue  # skip blank lines
    if items[0] == '#': continue  # skip header lines
    x_c.append(float(items[0]))
    y_c.append(float(items[1]))
bezCtrlPts = [Vector(x_c[i],y_c[i]) for i in range(1, len(x_c))]

# Specify dimensions for the throat and nozzle blocks.
Lthr = 0.0125   # Length of throat block
Lnoz = x_c[-1]  # Length of nozzle block

# Create nodes for the throat block (constant-area section of nozzle).
throat_axis = Node(-Lthr, 0.0)
throat_wall = Node(-Lthr, y_c[0])

# Create nodes for the nozzle block (diverging section of nozzle). 
divergence_axis = Node(0.0, 0.0)
divergence_wall = Node(0.0, y_c[0])
nozzle_end_axis = Node(x_c[-1], 0.0)
nozzle_end_wall = Node(x_c[-1], y_c[-1])

# Define paths and make patch for the throat block.
t_south = Line(throat_axis, divergence_axis)
t_north = Line(throat_wall, divergence_wall)
t_west = Line(throat_axis, throat_wall)
t_east = Line(divergence_axis, divergence_wall)
throat_region = make_patch(t_north, t_east, t_south, t_west)

# Define paths and make patch for the nozzle block.
n_north = Polyline([Line(divergence_wall,bezCtrlPts[0]), Bezier(bezCtrlPts)])
n_west = Line(divergence_axis, divergence_wall)
n_east = Line(nozzle_end_axis, nozzle_end_wall)
n_south = Line(divergence_axis, nozzle_end_axis)
expansion_region = ExpandingChannelPatch(n_south, n_north, n_west, n_east)
    
# Define grid clustering for throat and nozzle blocks.
bx = 1.1     # clustering in the axial direction
by = 1.02    # clustering in the radial direction
bx_throat = 2.01183994581  # axial clustering for throat block
x_clust = RobertsClusterFunction(1, 1, bx)
y_clust = RobertsClusterFunction(0, 1, by)
x_clust_throat = RobertsClusterFunction(1, 1, bx_throat)

# Create blocks.
nnj = 40; nbj = 2
throat_blk = SuperBlock2D(throat_region, nni=38, nnj=nnj, nbi=1, nbj=nbj,
                          bc_list=[FixedTBC(300.0), ExtrapolateOutBC(), 
                                   SlipWallBC(), SupInBC(inflow)],
                          cf_list=[x_clust_throat, y_clust, x_clust_throat, y_clust],
                          fill_condition=inflow, label="throat")

nozzle_blk = SuperBlock2D(expansion_region, nni=600, nnj=nnj, nbi=30, nbj=nbj,
                          bc_list=[FixedTBC(300.0), ExtrapolateOutBC(),
                                   SlipWallBC(), SupInBC(inflow)],
                          cf_list=[x_clust, y_clust, x_clust, y_clust],
                          fill_condition=inflow, label="nozzle")

identify_block_connections()

# Global attributes for the simulation.
gdata.dimensions = 2
gdata.axisymmetric_flag = 1
gdata.viscous_flag = 1
gdata.turbulence_model = 'None'
gdata.flux_calc = ADAPTIVE
gdata.max_time = 0.001  # seconds
gdata.max_step = 800000
gdata.dt = 1.0e-11
gdata.cfl = 0.4
gdata.cfl_count = 5
gdata.dt_plot = 1.0e6   # Set as large value for block-marching mode
gdata.dt_history = 10.0e-6

