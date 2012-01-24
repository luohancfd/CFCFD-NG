## \file turb_flat_plate.py
## \brief Turbulent flow over a flat plate 
##         - Mabey test case (AGARDograph 223 - Test series 7402)
##           (Referenced from Fernholz & Finley (1977), 
##           AGARDograph No. 223, "A critical compilation of 
##           compressible turbulent boundary layer data.")
##         - k-omega turbulence model
##
## \author Wilson Chan 18 June 2010
##         PJ, clean up for barrine, November 2010

from cfpylib.gasdyn import sutherland

gdata.title = "Mabey's Mach 4.5 flow over a flat plate (k-omega)"
print gdata.title
gdata.turbulence_flag = 1 
gdata.turbulence_model = "k_omega"
gdata.viscous_flag = 1
gdata.flux_calc = AUSMDV
        
gdata.max_time = 1.6e-3  # About 3 flow lengths (1 flow length ~ 0.56 ms)  
gdata.dt_plot =  0.2e-3
gdata.dt_history = 1.0e-3
gdata.max_step = 3000000

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-9  # Only an initial guess - the simulation will take over

select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Mabey's data set 74021802
p_inf = 3.16e3  # Pa
u_inf = 712.9   # m/s
T_inf = 62.16   # degrees K
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
L = 0.4 # metres
H = 0.4 * L

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
                   nni=128, nnj=96, nbi=4, nbj=4, 
                   fill_condition=initial,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.0014),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.0074)],
                   bc_list=[AdiabaticBC(), ExtrapolateOutBC(),
                            SupInBC(inflow), SupInBC(inflow)])

