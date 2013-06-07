# cyl-flare.py
# PJ, 11-May-2013, 29-May-2013
# Model of the CUBRC hollow cylinder with extended-flare experiment. 

gdata.title = "Hollow cylinder with extended flare."
print gdata.title
gdata.axisymmetric_flag = 1
# Conditions to match those reported for CUBRC Run 14
p_inf = 31.88 # Pa
u_inf = 2304.0 # m/s
T_inf = 120.4 # degree K
T_vib = 2467.0 # degrees K (but we will ignore for ideal-gas)
select_gas_model(model='ideal gas', species=['N2'])
inflow  = FlowCondition(p=p_inf, u=u_inf, T=T_inf)
initial = FlowCondition(p=p_inf/5, u=0, T=T_inf)
T_wall = 295.2 # degree K

mm = 1.0e-3 # metres per mm
L1 = 101.7*mm # cylinder length
L2 = 220.0*mm # distance to end of flare
R1 = 32.5*mm
alpha = 30.0*math.pi/180.0 # angle of flare
tan_alpha = math.tan(alpha)
a0 = Vector(0.0, R1); a1 = a0+Vector(0.0,5*mm) # leading edge of cylinder
b0 = Vector(L1, R1); b1 = b0+Vector(-5*mm,20*mm) # start flare
c0 = Vector(L2, R1+tan_alpha*(L2-L1)); c1 = c0+Vector(0.0,25*mm) # end flare

rcfx = RobertsClusterFunction(1,0,1.2)
rcfy = RobertsClusterFunction(1,0,1.1)
ni0 = 200; nj0 = 80 # We'll scale discretization off these values
factor = 1
ni0 *= factor; nj0 *= factor

cyl = SuperBlock2D(CoonsPatch(a0,b0,b1,a1),
                   nni=ni0, nnj=nj0, nbi=6, nbj=2,
                   bc_list=[SupInBC(inflow),None,
                            FixedTBC(T_wall),SupInBC(inflow)],
                   cf_list=[rcfx,rcfy,rcfx,rcfy],
                   fill_condition=inflow, label="cyl")
flare = SuperBlock2D(CoonsPatch(b0,c0,c1,b1),
                     nni=ni0, nnj=nj0, nbi=6, nbj=2,
                     bc_list=[SupInBC(inflow),FixedPOutBC(p_inf/5),
                              FixedTBC(T_wall),None],
                     cf_list=[None,rcfy,None,rcfy],
                     fill_condition=initial, label="fl")
identify_block_connections()

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.0
# The settling of the separation bubble will probably dominate.
gdata.max_time = 2.5e-3 # long enough, looking at earlier simulations
gdata.max_step = 2000000
gdata.dt = 1.0e-9
gdata.dt_plot = 0.25e-3

sketch.xaxis(0.0, 0.250, 0.05, -0.010)
sketch.yaxis(0.0, 0.250, 0.05, -0.010)
sketch.window(0.0, 0.0, 0.250, 0.250, 0.05, 0.05, 0.25, 0.25)
