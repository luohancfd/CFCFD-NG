# couette.py
# turbulent plane couette flow

gdata.dimensions = 2
select_gas_model(model='ideal gas', species=['air'])
gdata.title = "Telbany's experiment"
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.turbulence_model = "k_omega"
gdata.max_time = 5.0
gdata.max_step = 8000000000
gdata.dt = 1.0e-10
gdata.dt_plot = 0.5
gdata.dt_history = 1e-2
gdata.print_count = 20000

# -------------- Some useful functions -----------------------------
def sutherland_viscosity(T, gas_type):
    "Sutherland's viscosity law"
    # Available gas types : Air, N2, O2, H2, CO2, CO, Ar 
    # Reference : White FM (2006) Viscous Fluid Flow, pp.28
    if gas_type == "Air":
        mu_ref = 1.716e-5; T_ref = 273.0; S = 111.0
    elif gas_type == "N2":
        mu_ref = 1.663e-5; T_ref = 273.0; S = 107.0
    elif gas_type == "O2":
        mu_ref = 1.919e-5; T_ref = 273.0; S = 139.0
    elif gas_type == "H2":
        mu_ref = 8.411e-6; T_ref = 273.0; S = 97.0
    elif gas_type == "CO2":
        mu_ref = 1.370e-5; T_ref = 273.0; S = 222.0
    elif gas_type == "CO":
        mu_ref = 1.657e-5; T_ref = 273.0; S = 136.0
    elif gas_type == "Ar":
        mu_ref = 1.716e-5; T_ref = 273.0; S = 144.0
    else:
        print "ERROR! Gas type not recognised for this function!"
    return mu_ref * (T / T_ref)**1.5 * (T_ref + S)/(T + S);
# ------------------------------------------------------------------


p_inf = 101325.0  # Pa
T_inf = 293.0
rho_inf = p_inf / (287.0 * T_inf)  # kg/m**3
u_max = 12.84    # m/s
I_turb = 0.1
u_turb_lam = 100.0
tke_inf = 1.5*(I_turb*u_max)**2
mu_t_inf = u_turb_lam * sutherland_viscosity(T_inf, "Air")
omega_inf = rho_inf * tke_inf / mu_t_inf

print "tke_inf is", tke_inf
print "omega_inf is", omega_inf

h0 = 0.0
h1 = 66.0e-3;

l0 = 0.0
l1 = 2440.0e-3

a = Node(l0,h0)
b = Node(l0,h1)
c = Node(l1,h0)
d = Node(l1,h1)

ab = Line(a,b); cd = Line(c,d); 
ac = Line(a,c); bd = Line(b,d); 

HistoryLocation(l1/2.0,h1)
HistoryLocation(l1/1.5,h1/2.0)
HistoryLocation(l1,0.0)

initial = FlowCondition(p=p_inf, u=u_max, v=0.0, tke=tke_inf, omega=omega_inf, T=T_inf)

mv = MovingWallBC(v_trans=[u_max,0.0,0.0], Twall_flag=True, Twall=T_inf)
ad = FixedTBC(Twall=T_inf)
inletoutlet = InletOutletBC(Pout=p_inf, I_turb=I_turb, u_turb_lam=u_turb_lam)

bc_list0 = [ad,inletoutlet,mv,inletoutlet]
nx0 = 480; ny0 = 160;

blk0 = SuperBlock2D(make_patch(bd, cd, ac, ab), 
              nni=nx0, nnj=ny0,
              nbi=48, nbj=1,
              fill_condition=initial,
              bc_list=bc_list0,
              cf_list=[RobertsClusterFunction(1,1,1.2),
                       RobertsClusterFunction(1,1,1.05),
                       RobertsClusterFunction(1,1,1.2),
                       RobertsClusterFunction(1,1,1.05)])

identify_block_connections()

# The following scales provide a reasonable picture.
sketch.xaxis(l0, l1, 0.02, -0.005)
sketch.yaxis(h0, h1, 0.2, -0.004)
sketch.window(l0, h0, l1, h1, 0.05, 0.05, 0.15, 0.075)

