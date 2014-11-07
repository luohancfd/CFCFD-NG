# rfp.py
gdata.title = "Rarefied flow over a flat plate."
# PJ, 2014-Nov-06
# To try out the JumpWallBC.

print gdata.title
gdata.viscous_flag = 1
gdata.max_time = 100.0e-6  # will allow several flow lengths   
gdata.dt_plot =  gdata.max_time / 10
gdata.dt_history = 1.0e-6
gdata.max_step = 3000000
gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-12  # only an initial guess - the simulation will take over

select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Schetz's worked example 5-1
p_inf = 10.0  # Pa
u_inf = 715.0   # m/s
T_inf = 50.0    # degrees K

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,])
initial = FlowCondition(p=p_inf/3, u=0.0, v=0.0, T=T_inf*3.0, massf=[1.0,])

# Geometry of plate and flow domain.
mm = 1.0e-3 # m
L = 20.0*mm   # Length of flat plate
H = 0.5 * L   # Height of flow domain

#                 p11
#    flow=>    -/- |
#           -/-    |
#       p01        |
# flow=> |         |
#       p00-------p10 ----> x
#            wall
# 
p01 = Node(0.0, H/4.0, label='p01'); p11 = Node(L, H, label='p11')
p00 = Node(0.0, 0.0, label='p00'); p10 = Node(L, 0.0, label='p10')

# Split the domain into a grid of blocks so that we can make use
# of our cluster computer.
factor = 1
ni = int(40 * factor); nj = int(40 * factor)
if 1:
    south_boundary = FixedTBC(300.0)
else:
    south_boundary = JumpWallBC(300, 1.0)

blk = SuperBlock2D(AOPatch(p00,p10,p11,p01), 
                   nni=ni, nnj=nj, nbi=2, nbj=2, 
                   fill_condition=initial,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(1,0, 1.016),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(1,0, 1.05)],
                   bc_list=[SupInBC(inflow), ExtrapolateOutBC(),
                            south_boundary, SupInBC(inflow)])

# Set to make a nice 2D rendering of the blocks.
sketch.xaxis(0.0, 0.025, 0.005, -0.0001)
sketch.yaxis(0.0, 0.025, 0.005, -0.0001)
sketch.window(0.0, 0.0, 0.025, 0.025, 0.05, 0.05, 0.17, 0.17)

