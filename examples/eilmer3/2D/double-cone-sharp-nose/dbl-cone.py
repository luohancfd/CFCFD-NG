# dbl-cone.py
# PJ, 12-June-2013
# Model of the CUBRC double-cone with sharp nose. 

gdata.title = "Double-cone, sharp nose."
print gdata.title
gdata.axisymmetric_flag = 1
# Conditions to match those reported for CUBRC Run 35
p_inf = 18.55 # Pa
u_inf = 2576.0 # m/s
T_inf = 102.2 # degree K
T_vib = 2711.0 # degrees K (but we will ignore for ideal-gas)
select_gas_model(model='ideal gas', species=['N2'])
inflow  = FlowCondition(p=p_inf, u=u_inf, T=T_inf)
initial = FlowCondition(p=p_inf/5, u=0, T=T_inf)
T_wall = 295.8 # degree K

mm = 1.0e-3 # metres per mm
a0 = Vector(0.0, 0.0)
a1 = Vector(0.0,5*mm) # leading edge of domain
b0 = Vector(92.08*mm,42.94*mm) # junction between cones
b1 = Vector(76*mm,61*mm) # out in the free stream
c0 = Vector(153.69*mm,130.925*mm) # downstream-edge of second cone
c1 = Vector(124*mm,181*mm) # out in the free stream
d0 = Vector(193.68*mm,130.925*mm) # down-stream edge of domain
d1 = Vector(193.68*mm,181*mm)

rcfx = RobertsClusterFunction(1,0,1.2)
rcfy = RobertsClusterFunction(1,0,1.1)
ni0 = 120; nj0 = 40 # We'll scale discretization off these values
factor = 2
ni0 *= factor; nj0 *= factor

cone1 = SuperBlock2D(CoonsPatch(a0,b0,b1,a1),
                     nni=ni0, nnj=nj0, nbi=6, nbj=2,
                     bc_list=[SupInBC(inflow),None,
                              FixedTBC(T_wall),SupInBC(inflow)],
                     cf_list=[rcfx,rcfy,rcfx,rcfy],
                     fill_condition=inflow, label="cone1")
cone2 = SuperBlock2D(CoonsPatch(b0,c0,c1,b1),
                     nni=ni0, nnj=nj0, nbi=6, nbj=2,
                     bc_list=[SupInBC(inflow),None,
                              FixedTBC(T_wall),None],
                     cf_list=[None,rcfy,None,rcfy],
                     fill_condition=initial, label="cone2")
cone3 = SuperBlock2D(CoonsPatch(c0,d0,d1,c1),
                     nni=int(ni0/2), nnj=nj0, nbi=2, nbj=2,
                     bc_list=[SupInBC(inflow),FixedPOutBC(p_inf/5),
                              FixedTBC(T_wall),None],
                     cf_list=[None,rcfy,None,rcfy],
                     fill_condition=initial, label="out")
identify_block_connections()

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.0
# The settling of the separation bubble will probably dominate.
gdata.max_time = 5.0e-3 # long enough, maybe
gdata.max_step = 4000000
gdata.dt = 1.0e-9
gdata.dt_plot = 0.25e-3

sketch.xaxis(0.0, 0.250, 0.05, -0.010)
sketch.yaxis(0.0, 0.250, 0.05, -0.010)
sketch.window(0.0, 0.0, 0.250, 0.250, 0.05, 0.05, 0.25, 0.25)
