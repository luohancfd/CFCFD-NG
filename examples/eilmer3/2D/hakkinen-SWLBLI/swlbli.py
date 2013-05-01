# swlbli.py
# PJ, 01-May-2013
# Model of Hakkinen et al's 1959 experiment. 

gdata.title = "Shock-wave laminar-boundary-layer interaction."
print gdata.title
# Conditions to match those of Figure 6: pf/p0=1.4, Re_shock=2.96e5
p_inf = 6205.0 # Pa
u_inf = 514.0 # m/s
T_inf = 164.4 # degree K

# Accept defaults giving perfect air (R=287 J/kg.K, gamma=1.4)
select_gas_model(model='ideal gas', species=['air'])
inflow  = FlowCondition(p=p_inf, u=u_inf, T=T_inf)

mm = 1.0e-3 # metres per mm
L1 = 10.0*mm; L2 = 90.0*mm; L3 = 67*mm
H1 = 37.36*mm
alpha = 3.09*math.pi/180.0 # angle of inviscid shock generator
tan_alpha = math.tan(alpha)
a0 = Vector(-L1, 0.0); a1 = a0+Vector(0.0,H1) # leading edge of shock generator
b0 = Vector(0.0, 0.0); b1 = b0+Vector(0.0,H1-L1*tan_alpha) # start plate
c0 = Vector(L3, 0.0); c1 = c0+Vector(0.0,H1-(L1+L2)*tan_alpha) # end shock generator
d0 = Vector(L2, 0.0); d1 = d0+Vector(0.0,H1) # end plate

# The following lists are in order [N, E, S, W]
rcf = RobertsClusterFunction(1,1,1.1)
ni0 = 20; nj0 = 80 # We'll scale discretization off these values

inlet = SuperBlock2D(CoonsPatch(a0,b0,b1,a1),
                     nni=ni0, nnj=nj0, nbi=1, nbj=2,
                     bc_list=[SlipWallBC(),None,SlipWallBC(),SupInBC(inflow)],
                     cf_list=[None,rcf,None,rcf],
                     fill_condition=inflow, label="in")
plate1 = SuperBlock2D(CoonsPatch(b0,c0,c1,b1),
                      nni=ni0*7, nnj=nj0, nbi=7, nbj=2,
                      bc_list=[SlipWallBC(),None,AdiabaticBC(),None],
                      cf_list=[None,rcf,None,rcf],
                      fill_condition=inflow, label="p1")
plate2 = SuperBlock2D(CoonsPatch(c0,d0,d1,c1),
                      nni=ni0*2, nnj=nj0, nbi=2, nbj=2,
                      bc_list=[SlipWallBC(),FixedPOutBC(6205.0),AdiabaticBC(),None],
                      cf_list=[None,rcf,None,rcf],
                      fill_condition=inflow, label="p2")
identify_block_connections()

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.0
gdata.max_time = 5.0*L2/u_inf  # in flow lengths
gdata.max_step = 200000
gdata.dt = 1.0e-8
gdata.dt_plot = gdata.max_time/10

sketch.xaxis(-0.020, 0.100, 0.020, -0.010)
sketch.yaxis(0.000, 0.040, 0.020, -0.004)
sketch.window(-0.02, 0.0, 0.10, 0.12, 0.05, 0.05, 0.25, 0.25)
