## \file flat_runup.py
## \brief Turbulent flow over a flat plate 
##         - Built off turb-flat-plate.py by W. Chan & P. Jacobs
##         - k-omega turbulence model
##
## \author Samuel Stennett 25 September 2014


from cfpylib.gasdyn import sutherland

gdata.title = "Test Case 3 runup - flow over a flat plate (k-omega)"
print gdata.title
gdata.turbulence_model = "k_omega"
gdata.viscous_flag = 1
gdata.flux_calc = AUSMDV
        
gdata.max_time = 5.0e-3  # About 3 flow lengths (1 flow length ~ 1.78 ms)
gdata.dt_plot =  1.0e-4
gdata.dt_history = 5.0e-4
gdata.max_step = 3000000

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-9  # Only an initial guess - the simulation will take over

select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Kim's experimental setup
p_inf = 7.1e3  # Pa
u_inf = 671.92   # m/s (Mach 4 (3.98))
T_inf = 70.3   # degrees K
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
L = 1.60 # metres
H = 0.1524

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
tc3_runup_blk = SuperBlock2D(make_patch(north, east, south, west), 
                   nni=512, nnj=152, nbi=4, nbj=4, 
                   fill_condition=initial,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(1,0, 1.0013),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(1,0, 1.0013)],
                   bc_list=[ExtrapolateOutBC(),ExtrapolateOutBC(),AdiabaticBC(), SupInBC(inflow)])

