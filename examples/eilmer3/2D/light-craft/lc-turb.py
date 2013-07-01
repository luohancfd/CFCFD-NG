# lc.py
# PJ, 22-June-2013, 28-June-2013 turbulent
# David Buttsworth and Con Doolan's Light Craft. 

gdata.title = "Light Craft."
print gdata.title
gdata.axisymmetric_flag = 1
# Conditions to match the TUSQ operating condition.
p_inf = 701.5 # Pa
u_inf = 989.8 # m/s
T_inf = 72.46 # degree K
p_init = 1000.0 # Pa, DRB estimate
T_wall = 300.0 # degree K

# Estimate turbulence quantities for free stream
# by specifying the intensity as 0.01 and estimating the
# turbulence viscosity as 10 times the laminar viscosity.
tke_inf = 1.5 * (u_inf * 0.01)**2
rho_inf = p_inf / (287.0 * T_inf)
from cfpylib.gasdyn.sutherland import mu
mu_t_inf = 10.0 * mu(T_inf, 'Air')
omega_inf = rho_inf * tke_inf / mu_t_inf
gdata.turbulence_model = "k_omega"
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

select_gas_model(model='ideal gas', species=['air'])
inflow  = FlowCondition(p=p_inf, u=u_inf, T=T_inf, tke=tke_inf, omega=omega_inf)
initial = FlowCondition(p=p_init, u=0, T=T_wall, tke=0.0, omega=1.0)

mm = 1.0e-3 # metres per mm

# Conical forebody
a0 = Vector(0.0, 0.0)
a1 = Vector(0.0,10.0)*mm # leading edge of domain
y_outer = 130.0 # mm outer radius of simulated domain
a2 = Vector(0.0,y_outer)*mm
alpha = 0.5*24.57*math.pi/180.0 # cone half-angle
f0 = Vector(1.0,math.tan(alpha))*261.7*mm
f1 = Vector(261.7,74.9)*mm
f2 = Vector(261.7,y_outer)*mm

rcfx0 = RobertsClusterFunction(1,1,1.2)
rcfx1 = RobertsClusterFunction(1,0,1.2)
rcfx2 = RobertsClusterFunction(1,0,1.1)
rcfy0 = RobertsClusterFunction(1,0,1.1)
rcfy1 = RobertsClusterFunction(1,1,1.1)
rcfy2 = RobertsClusterFunction(0,1,1.1)
ni0 = 240; nj0 = 20 # We'll scale discretization off these values
factor = 2.0
ni0 = int(factor*ni0); nj0 = int(factor*nj0)

blk_0 = SuperBlock2D(CoonsPatch(a0,f0,f1,a1),
                     nni=ni0, nnj=nj0, nbi=6, nbj=1,
                     bc_list=[None,None,
                              FixedTBC(T_wall),SupInBC(inflow)],
                     cf_list=[rcfx0,rcfy1,rcfx0,rcfy0],
                     fill_condition=initial, label="cone1")
blk_1 = SuperBlock2D(CoonsPatch(a1,f1,f2,a2),
                     nni=ni0, nnj=int(nj0*1.5), nbi=6, nbj=1,
                     bc_list=[SupInBC(inflow),None,
                              None,SupInBC(inflow)],
                     cf_list=[rcfx0,rcfy0,rcfx0,rcfy0],
                     fill_condition=initial, label="cone2")

b = Vector(1.0,math.tan(alpha))*285.7*mm 
c_test = Vector(333.3,55.0)*mm # to check against trimmed spline
d_test = Vector(347.57,40.0)*mm # to check against trimmed spline
e = Vector(347.57,0.0)*mm
h = Vector(333.3,65.4)*mm
g0 = Vector(381.0,0.0)*mm
g1 = Vector(381.0,40.0)*mm
g2 = Vector(381.0,70.0)*mm
g3 = Vector(381.0,83.3)*mm
g4 = Vector(381.0,y_outer)*mm

# Cowl, first spline, combustor
spl3 = Spline([Vector(261.7+0,74.9)*mm,
               Vector(261.7+10.3,74.0)*mm,
               Vector(261.7+20.5,72.9)*mm,
               Vector(261.7+30.8,71.5)*mm,
               Vector(261.7+41.0,70.0)*mm,
               Vector(261.7+51.2,68.5)*mm,
               Vector(261.7+61.4,67.0)*mm,
               Vector(261.7+71.6,65.4)*mm])
# Cowl, second spline, expansion
spl4 = Spline([Vector(261.7+71.6,65.4)*mm,
               Vector(261.7+77.7,69.6)*mm,
               Vector(261.7+83.9,73.3)*mm,
               Vector(261.7+90.5,76.4)*mm,
               Vector(261.7+97.3,79.0)*mm,
               Vector(261.7+104.3,81.0)*mm,
               Vector(261.7+111.7,82.4)*mm,
               Vector(261.7+119.3,83.3)*mm])
# Body, spline 1, combustor
spl1 = Spline([Vector(285.7,62.2)*mm,
               Vector(292.1,61.8)*mm,
               Vector(299.0,60.8)*mm,
               Vector(305.9,59.7)*mm,
               Vector(312.8,58.6)*mm,
               Vector(319.6,57.4)*mm,
               Vector(326.5,56.2)*mm,
               Vector(333.3,55.0)*mm])
# Body, spline 2, expansion
spl2 = Spline([Vector(333.3,55.0)*mm,
               Vector(347.57,40.0)*mm,
               Vector(350.0,37.8)*mm,
               Vector(370.0,23.79)*mm,
               Vector(390.0,14.2)*mm,
               Vector(410.0,7.73)*mm,
               Vector(430.0,3.51)*mm,
               Vector(450.0,1.07)*mm,
               Vector(476.2,0.0)*mm])
# On the body, we want only a little bit of spline 2,
# just along to the position x=347.57mm y=40.0mm.
spl2b = spl2.copy(direction=-1) # going up the west face
def my_error(t, path=spl2b): return abs(path.eval(t).y - 0.040)
from cfpylib.nm import zero_solvers
t_trim = zero_solvers.secant(my_error, 0.6, 0.61)
if t_trim == 'FAIL':
    print "Failed to solve intersection with spline 2"
    sys.exit()
spl2b.t0 = t_trim
c = spl2b.eval(1.0)
d = spl2b.eval(0.0)
assert vabs(c - c_test) < 1.0e-6
assert vabs(d - d_test) < 1.0e-6

# Combustor
blk_2 = SuperBlock2D(CoonsPatch(Polyline([Line(f0,b),spl1]), spl3,
                                Line(f0,f1), Line(c,h)),
                     nni=int(ni0/2), nnj=nj0, nbi=3, nbj=1,
                     bc_list=[FixedTBC(T_wall),None,
                              FixedTBC(T_wall),None],
                     cf_list=[None,rcfy1,None,rcfy1],
                     fill_condition=initial, label="comb1")

# Outside cowl
blk_3 = SuperBlock2D(CoonsPatch(f1,g3,g4,f2),
                     nni=int(ni0/2), nnj=int(nj0*1.5), nbi=2, nbj=1,
                     bc_list=[SupInBC(inflow),None,
                              FixedTBC(T_wall),None],
                     cf_list=[None,rcfy0,None,rcfy0],
                     fill_condition=initial, label="comb2")

# Tail and expansion
blk_4 = SuperBlock2D(CoonsPatch(e,g0,g1,d),
                     nni=int(ni0/4), nnj=nj0, nbi=1, nbj=1,
                     bc_list=[None,None,
                              SlipWallBC(),FixedTBC(T_wall)],
                     cf_list=[rcfx1,None,rcfx1,None],
                     fill_condition=initial, label="tail0")
blk_5 = SuperBlock2D(CoonsPatch(Line(d,g1), Line(c,g2),
                                spl2b, Line(g1,g2)),
                     nni=int(ni0/4), nnj=nj0, nbi=1, nbj=1,
                     bc_list=[None,None,
                              None,FixedTBC(T_wall)],
                     cf_list=[rcfx1,None,rcfx1,rcfy2],
                     fill_condition=initial, label="tail1")
blk_6 = SuperBlock2D(CoonsPatch(Line(c,g2), spl4,
                                Line(c,h), Line(g2,g3)),
                     nni=int(ni0/4), nnj=int(nj0), nbi=1, nbj=1,
                     bc_list=[FixedTBC(T_wall),None,
                              None,None],
                     cf_list=[rcfx1,rcfy2,rcfx1,rcfy1],
                     fill_condition=initial, label="tail2")

x_exit = 600.0
k0 = Vector(x_exit,0.0)*mm
k1 = Vector(x_exit,40.0)*mm
k2 = Vector(x_exit,70.0)*mm
k3 = Vector(x_exit,83.3)*mm
k4 = Vector(x_exit,y_outer)*mm

# Downstream region
blk_7 = SuperBlock2D(CoonsPatch(g0,k0,k1,g1),
                     nni=int(ni0/4), nnj=nj0, nbi=1, nbj=1,
                     bc_list=[None,FixedPOutBC(p_init),
                              SlipWallBC(),None],
                     cf_list=[rcfx2,None,rcfx2,None],
                     fill_condition=initial, label="exit0")
blk_8 = SuperBlock2D(CoonsPatch(g1,k1,k2,g2),
                     nni=int(ni0/4), nnj=nj0, nbi=1, nbj=1,
                     bc_list=[None,FixedPOutBC(p_init),
                              None,None],
                     cf_list=[rcfx2,None,rcfx2,None],
                     fill_condition=initial, label="exit1")
blk_9 = SuperBlock2D(CoonsPatch(g2,k2,k3,g3),
                     nni=int(ni0/4), nnj=int(nj0), nbi=1, nbj=1,
                     bc_list=[None,FixedPOutBC(p_init),
                              None,None],
                     cf_list=[rcfx2,None,rcfx2,rcfy2],
                     fill_condition=initial, label="exit2")
blk_10 = SuperBlock2D(CoonsPatch(g3,k3,k4,g4),
                     nni=int(ni0/4), nnj=int(nj0*1.5), nbi=1, nbj=1,
                     bc_list=[SupInBC(inflow),FixedPOutBC(p_init),
                              None,None],
                     cf_list=[rcfx2,None,rcfx2,rcfy0],
                     fill_condition=initial, label="exit3")

identify_block_connections()

# Add the transducer locations as HistoryLocations
trans1 = Vector(1.0,math.tan(alpha))*231.62*mm
HistoryLocation(trans1.x, trans1.y, trans1.z, label='t1p0')
HistoryLocation(trans1.x, trans1.y, trans1.z, i_offset=1, label='t1p1')
HistoryLocation(trans1.x, trans1.y, trans1.z, i_offset=-1, label='t1m1')
trans2 = Vector(1.0,math.tan(alpha))*251.16*mm
HistoryLocation(trans2.x, trans2.y, trans2.z, label='t2p0')
HistoryLocation(trans2.x, trans2.y, trans2.z, i_offset=1, label='t2p1')
HistoryLocation(trans2.x, trans2.y, trans2.z, i_offset=-1, label='t2m1')
trans3 = Vector(1.0,math.tan(alpha))*270.70*mm
HistoryLocation(trans3.x, trans3.y, trans3.z, label='t3p0')
HistoryLocation(trans3.x, trans3.y, trans3.z, i_offset=1, label='t3p1')
HistoryLocation(trans3.x, trans3.y, trans3.z, i_offset=-1, label='t3m1')
trans4 = Vector(301.33,60.0)*mm
HistoryLocation(trans4.x, trans4.y, trans4.z, label='t4p0')
HistoryLocation(trans4.x, trans4.y, trans4.z, i_offset=1, label='t4p1')
HistoryLocation(trans4.x, trans4.y, trans4.z, i_offset=-1, label='t4m1')
trans5 = Vector(321.06,57.0)*mm
HistoryLocation(trans5.x, trans5.y, trans5.z, label='t5p0')
HistoryLocation(trans5.x, trans5.y, trans5.z, i_offset=1, label='t5p1')
HistoryLocation(trans5.x, trans5.y, trans5.z, i_offset=-1, label='t5m1')

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.0
# TUSQ runs for up to 250ms.
gdata.max_time = 5.0e-3 # long enough to test unsteadiness
gdata.max_step = 800000 
gdata.dt = 1.0e-8
gdata.dt_plot = 0.25e-3
gdata.dt_history = 1.0e-6

sketch.xaxis(0.0, 0.50, 0.1, -0.010)
sketch.yaxis(0.0, 0.50, 0.1, -0.010)
sketch.window(0.0, 0.0, 0.50, 0.50, 0.05, 0.05, 0.25, 0.25)
