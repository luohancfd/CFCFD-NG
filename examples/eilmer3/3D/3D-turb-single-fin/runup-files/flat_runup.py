## \file flat_runup.py
## \brief Turbulent flow over a flat plate 
##         - Built off turb-flat-plate.py by W. Chan & P. Jacobs
##         - k-omega turbulence model
##
## \author Samuel Stennett 6 July 2014


from cfpylib.gasdyn import sutherland

gdata.title = "Test Case 2 runup - flow over a flat plate (k-omega)"
print gdata.title
gdata.turbulence_model = "k_omega"
gdata.viscous_flag = 1
gdata.flux_calc = AUSMDV
        
gdata.max_time = 0.92e-3  # About 3 flow lengths (1 flow length ~ 0.31 ms)
gdata.dt_plot =  0.1e-3
gdata.dt_history = 1.0e-3
gdata.max_step = 30000000

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-9  # Only an initial guess - the simulation will take over

select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Kim's experimental setup
p_inf = 10.3e3  # Pa
u_inf = 707.0   # m/s (Mach 4 (3.98))
T_inf = 70.4   # degrees K
rho_inf = p_inf / (287.1 * T_inf)

# Estimate turbulence quantities for free stream by specifying turbulence
# intensity and the ratio of turbulent-to-laminar viscosity.
turbulence_intensity = 0.01
turb_lam_viscosity_ratio = 1.0
tke_inf = 1.5 * (turbulence_intensity * u_inf)**2
mu_t_inf = turb_lam_viscosity_ratio * sutherland.mu(T_inf, "Air")
omega_inf = rho_inf * tke_inf / mu_t_inf
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,],
                       tke=tke_inf, omega=omega_inf )
initial = inflow

# Geometry of plate and flow domain.
L = 0.240 # metres
H = 0.02769

#        d---------c
# flow=> |         |
#        |         |
#        |         |
#        |         |
#        a---------b ----> x
# 	    wall
#

a = Node(0.0, 0.0); b = Node(L, 0.0); c = Node(L, H); d = Node(0.0, H)
north = Line(d,c); east = Line(b,c); south = Line(a,b); west = Line(a,d)

# Define the blocks, boundary conditions and set the discretisation.
tc2_runup_blk = SuperBlock2D(make_patch(north, east, south, west), 
                   nni=256, nnj=111, nbi=8, nbj=8, 
                   fill_condition=initial,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(1,0, 1.001),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(1,0, 1.001)],
                   bc_list=[ExtrapolateOutBC(),ExtrapolateOutBC(),AdiabaticBC(), SupInBC(inflow)])

