# Micro-combustion
# 5 mm x 0.3 mm channel
# methane/air, V = 40 cm/s, Phi = 1.0
# Xin Kang

gdata.title = "full channel simulation"
gdata.dimensions  = 2
gdata.axisymmetric_flag = 0
gdata.viscous_flag = 1


# Gas Model Set-up

select_gas_model(fname='thermally-perfect-drm19.lua')
set_reaction_scheme("drm19.lua",reacting_flag =1)
gdata.diffusion_model = 'ConstantLewisNumber'
gdata.diffusion_lewis_number = 1.0
gdata.diffusion_flag = 1


# Flow Conditions

molef = {'N2':7.52, 'O2':2.0, 'CH4':1.0}
gmodel = get_gas_model_ptr()
massf=gmodel.to_massf(molef)
P_exit = 1.01325e5 #Pa, 1 atm
T0 = 300.0 #K, total inlet stagnation temperature
u0 = 0.4 #m/s, expected combustor inlet velocity
inflow = FlowCondition(p=P_exit,T=T0,u=u0,v=0,massf=massf)


# Geometry

L = 5.0e-3 #m, length of the channel
h = 0.3e-3 #m, full channel simulation

a = Node(0.0,0.0)
b = Node(0.0,h)
c = Node(L,0.0)
d = Node(L,h)

ab = Line(a,b)
ac = Line(a,c)
cd = Line(c,d)
bd = Line(b,d)


# Block Configuration
# Ensure only one blocks along y axis to facilitate udf lua function

nxcells0 = 390
nycells0 = 23
nbi0 = 64
nbj0 = 1


blk0 = SuperBlock2D(make_patch(bd, cd, ac, ab),
                    nni = nxcells0, nnj = nycells0, nbi=nbi0, nbj=nbj0,
                    bc_list = [UserDefinedBC(filename="udf-wall.lua", is_wall=1),\
                               FixedPOutBC(P_exit),\
                               UserDefinedBC(filename="udf-wall.lua", is_wall=1),\
                               UserDefinedBC(filename="udf-massflux-in.lua")],
                    fill_condition = inflow,
                    cf_list = [None, None, None, None])

# Make Block connections
identify_block_connections()


# IgnitionZone

point0 = Vector3(0.5*L, 0.0)
point1 = Vector3(point0.x + 0.05*L, h)
IgnitionZone(2000.0, point0, point1)
gdata.ignition_time_stop = 1.0e-3 # s


# History locations

HistoryLocation(0.0, 0.0)
HistoryLocation(0.0, h/2.0)
HistoryLocation(0.0, h)
HistoryLocation(L/2.0, 0.0)
HistoryLocation(L/2.0, h/2.0)
HistoryLocation(L/2.0, h)
HistoryLocation(L, 0.0)
HistoryLocation(L, h/2.0)
HistoryLocation(L, h)


# Simulation Parameters

gdata.flux_calc = AUSM_PLUS_UP
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 0.3
gdata.max_time = 60.0e-3  # seconds
gdata.max_step = 20000000
gdata.dt = 3.0e-11
gdata.dt_plot = gdata.max_time/600.0
gdata.dt_history = 1.0e-7
