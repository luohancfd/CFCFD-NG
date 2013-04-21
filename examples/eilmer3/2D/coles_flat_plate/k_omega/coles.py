## \file coles.py
## \brief Turbulent flow over a flat plate -- Coles test case
## \author PJ 
## \version Jan Pieter Nap, 1-Feb-2007, and again PJ October 2007
## \version for Elmer3, 31-Aug-2008
##

job_title = "Coles Mach 3.7 flow over a flat plate (k-omega)"
print job_title
gdata.dimensions = 2
# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Coles' data set 53010801
p_inf = 1.358e3  # Pa
u_inf = 677.4    # m/s
T_inf = 83.34    # degrees K

# Estimate turbulence quantities for free stream
# by specifying the intensity as 0.05 and estimating the
# turbulence viscosity as 100 times the laminar viscosity.
tke_inf = 1.5 * (u_inf * 0.05)**2
rho_inf = p_inf / (287.0 * T_inf)
def mu_air(T):
    "Sutherland expression for air viscosity."
    from math import sqrt
    mu_ref = 17.89e-6; T_ref = 273.1; S = 110.4
    T_T0 = T / T_ref
    return mu_ref * (T_ref + S)/(T + S) * T_T0 * sqrt(T_T0);
mu_t_inf = 100.0 * mu_air(T_inf)
omega_inf = rho_inf * tke_inf / mu_t_inf
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,],
                       tke=tke_inf, omega=omega_inf )
if 1:
    # It seems that everyone has been starting the flow field with
    # inflow conditions throughout and then letting the boundary layer
    # grow out into the main stream.
    initial = inflow
else:
    # Here is a low-pressure initial state more like a shock tunnel.
    # Unfortunately, it seems to play havoc with the turbulence.
    initial = FlowCondition(p= 0.1*p_inf, u=0.0, v=0.0, T=296.0, massf=[1.0,],
                            tke=tke_inf/100.0, omega=omega_inf/10.0 )


# Geometry of plate and flow domain.
L = 0.60 # metres
H = 0.4 * L
NB = 4 	 # number of blocks

#         wall
#        c---------b
# flow=> |         |
#        d         |
#          -\-     |
#    flow=>    -\- |
#        0         a ----> x
# 
a = Node(L, 0.0); b = Node(L, H); c = Node(0.0, H); d = Node(0.0, 3.0*H/4.0)
north = Line(c,b); east = Line(a,b); south = Line(d,a); west = Line(d,c)

# Define the blocks, boundary conditions and set the discretisation.
blk = SuperBlock2D(make_patch(north, east, south, west), 
                   nni=100, nnj=90, nbi=NB, nbj=1, 
                   fill_condition=initial,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.01),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.05)],
                   bc_list=[AdiabaticBC(), ExtrapolateOutBC(),
                            SupInBC(inflow), SupInBC(inflow)])


# Do a little more setting of global data.
gdata.turbulence_flag = 1  # to activate our turbulence model
gdata.turbulence_model = "k_omega"  # no need to specify; it is the default
gdata.title = job_title
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
 	
gdata.max_time = 5.0e-3  # should allow a few flow lengths   
gdata.dt_plot =  1.0e-3
gdata.dt_history = 1.0e-5
gdata.max_step = 3000000 

gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.0	
gdata.cfl_count = 3
gdata.stringent_cfl = 0 # 1 is more robust
gdata.dt = 1.0e-9	# only an initial guess, the simulation will take this over

# The following scales provide a reasonable picture.
sketch.xaxis(0.0, 0.6, 0.1, -0.05)
sketch.yaxis(0.0, 0.5, 0.1, -0.04)
sketch.window(0.0, 0.0, 0.6, 0.6, 0.05, 0.05, 0.17, 0.17)

