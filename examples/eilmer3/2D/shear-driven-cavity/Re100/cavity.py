# cavity.py
# The Ghia-Ghia-Shin shear-driven cavity, with near-incompressible flow.
# Peter J. 17-Mar-2014, 26-May-2014

gdata.title = "Shear-driven cavity flow Ghia et al. 1982"

select_gas_model(model='ideal gas', species=['air'])
Re = 100
mu = 1.847e-5 # Pa.s
u_lid = 1.0 # m/s
L = 1.0 # metre
rho = Re * mu / (u_lid * L) # kg/m**3
T = 300.0 # degrees K
print "rho=", rho, "kg/m**3"
p = rho * 287.1 * T  # Pa
print "p=", p, "Pa"
initial_flow = FlowCondition(p=p, T=T, u=0.0, v=0.0)

p00 = Vector(0.0,0.0); p10 = Vector(1.0,0.0)
p01 = Vector(0.0,1.0); p11 = Vector(1.0,1.0)
nn = 64

blk = SuperBlock2D(CoonsPatch(p00,p10,p11,p01), 
                   nni=nn, nnj=nn, nbi=2, nbj=2,
                   cf_list=4*[RobertsClusterFunction(1,1,1.1),],
                   fill_condition=initial_flow,
                   bc_list=[MovingWallBC(v_trans=[u_lid,0.0,0.0]),
                            AdiabaticBC(), AdiabaticBC(), AdiabaticBC()])
identify_block_connections()

gdata.viscous_flag = 1
gdata.gasdynamic_update_scheme = "euler"
gdata.flux_calc = ADAPTIVE
gdata.max_time = 20.0 # seconds; this is a slow flow
gdata.max_step = 5000000
gdata.dt = 1.0e-6
gdata.dt_plot = 500.0e-3
gdata.dt_history = 10.0e-3
# The following scales provide a reasonable picture.
sketch.xaxis(0.0, 1.0, 0.2, -0.1)
sketch.yaxis(0.0, 1.0, 0.2, -0.1)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.15, 0.15)

