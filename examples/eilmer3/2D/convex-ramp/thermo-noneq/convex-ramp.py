# convex-ramp.py
# PJ, Dan and Rowan, 15-Aug-2013
# Model of Mohammadian's convex-ramp experiment with thermal nonequilibrium.
# Revised 26-Aug-2013 to take his polynomial at face value.

gdata.title = "Mohammadian convex ramp, 2T thermo."
print gdata.title

# Gas-model
species = select_gas_model(model='two temperature gas', species=['N2','O2'])
gm = get_gas_model_ptr()
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()
# Energy exchange model (only if there are nonequilibrium temperatures)
if ntm > 1:
    set_energy_exchange_update("N2-O2-TV.lua")
# Conditions to match those reported in JFM paper with a guess for Tvib.
p_inf = 66.43 # Pa
u_inf = 1589.8 # m/s
# Temperatures
T_inf = [ 0.0 ] * ntm 
T_inf[0] = 41.92 # gas static temperature: degree K
Tv_inf = 1000.0 # NOTE: freestream vibrational temperature closer to stagnation T
for itm in range(1,ntm):
    T_inf[itm] = Tv_inf # nonequilibrium temperature: degree K
# Mass-fractions
massf_inf = [ 0.0 ] * nsp
massf_inf[species.index("N2")] = 0.767 # standard air
massf_inf[species.index("O2")] = 0.233 # standard air
#
inflow  = FlowCondition(p=p_inf, u=u_inf, T=T_inf, massf=massf_inf)
initial = FlowCondition(p=p_inf/5, u=0, T=T_inf, massf=massf_inf)
#
T_wall = 296.0 # degree K --assumed cold-wall temperature

# Mohammadian used the inch as his length scale.
m_per_inch = 0.0254
mm = 1.0e-3 # metres per mm

def ramp(t):
    """
    Parametric definition of ramp profile in xy-plane.

    Here, we join the initial straight 18-degree ramp to the polynomial.
    """
    alpha = 18.0*math.pi/180.0 # angle of initial straight section
    sin18 = math.sin(alpha)
    cos18 = math.cos(alpha)
    tan18 = math.tan(alpha)
    x_join_inch = 3.0
    y_join_inch = x_join_inch * tan18
    L1 = x_join_inch/cos18 # length of initial straight section
    L2 = 4.14677 # length of fairing (computed via maxima)
    t2 = (L1+L2) * t
    if t2 < L1:
        x_inch = t2 * cos18
        y_inch = t2 * sin18
    else:
        s = (t2 - L1)/L2 * 4.0577 
        g = 0.0026 * s**4 - 0.0211 * s**3
        x_inch = x_join_inch + s * cos18 - g * sin18
        y_inch = y_join_inch + s * sin18 + g * cos18
    return (x_inch*m_per_inch, y_inch*m_per_inch, 0.0)

# leading edge of ramp
x,y,z = ramp(0.0); a0 = Vector(x,y,z); a1 = a0+Vector(0.0,5*mm) 
# downstream end of transition curve
x,y,z = ramp(1.0); b0 = Vector(x,y,z); b1 = b0+Vector(-10.0*mm,40*mm)
# For the final straight section, angle continues at final angle of transition.
x_length = 10*m_per_inch - b0.x
beta = -1.90*math.pi/180.0
# end of model
c0 = Vector(b0.x+x_length,b0.y+x_length*math.tan(beta))
c1 = Vector(c0.x,b1.y)

rcfx = RobertsClusterFunction(1,0,1.2)
rcfy = RobertsClusterFunction(1,0,1.1)
ni0 = 200; nj0 = 40 # We'll scale discretization off these values
factor = 2.0
ni0 = int(ni0*factor); nj0 = int(nj0*factor)

wedge = SuperBlock2D(make_patch(Line(a1,b1),Line(b0,b1),
                                PyFunctionPath(ramp),Line(a0,a1)),
                     nni=ni0, nnj=nj0, nbi=10, nbj=2,
                     bc_list=[SupInBC(inflow),None,
                              FixedTBC(T_wall),SupInBC(inflow)],
                     cf_list=[rcfx,rcfy,rcfx,rcfy],
                     fill_condition=inflow, label="wedge")
tail = SuperBlock2D(CoonsPatch(b0,c0,c1,b1),
                    nni=int(ni0/4), nnj=nj0, nbi=4, nbj=2,
                    bc_list=[SupInBC(inflow),FixedPOutBC(p_inf/5),
                             FixedTBC(T_wall),None],
                    cf_list=[None,rcfy,None,rcfy],
                    fill_condition=initial, label="tail")
identify_block_connections()

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.0
gdata.max_time = 1.0e-3 # long enough for several flow lengths
gdata.max_step = 2000000
gdata.dt = 1.0e-9
gdata.dt_plot = 0.1e-3

sketch.xaxis(0.0, 0.20, 0.05, -0.010)
sketch.yaxis(0.0, 0.20, 0.05, -0.010)
sketch.window(0.0, 0.0, 0.20, 0.20, 0.05, 0.05, 0.25, 0.25)
