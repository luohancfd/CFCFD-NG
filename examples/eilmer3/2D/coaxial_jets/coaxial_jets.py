## \file modified and renamed coaxial_jets.py
## \brief Simulation of AD Cutler's coaxial mixing test case
## \author Wilson Chan
## \verson Wilson Chan May 2008
##         Wilson Chan 29 Apr 2009 Re-formatted to Elmer3 format
##         Wilson Chan 05 May 2009 Added nozzle and exterior parts
##                                 of the flow domain (The way the
##                                 blocks are specified is as per
##                                 logbook page 197.
##
##         Wilson Chan 11 Aug 2009 Refinement to fine level overdone
##                                 for the lip thickness region..
##                                 Will change number of cells in
##                                 that region back to medium level
##                                 and re-adjust grids to suit.
##	   
##	   Sam Stennett 3 Nov 2014 Edited case for upload to the
##				   Eilmer3 repository. Cleaned up
##				   custom functions.
##         
##         IMPORTANT: 	New version of Eilmer3 seems to have broken
##			seems to have broken the turbulence model used
##			in this case. Further investigations must be
##			done in order to determine the cause.
##				   

from cfpylib.gasdyn import sutherland

job_title = "Simulation of AD Cutler's coaxial mixing test case"
print job_title

# Set up global data
gdata.dimensions = 2
gdata.axisymmetric_flag = 1
gdata.viscous_flag = 1
#gdata.turbulence_model = "k_omega"
gdata.diffusion_flag = 1
gdata.diffusion_model = "FicksFirstLaw"
gdata.title = job_title
gdata.flux_calc = ADAPTIVE

gdata.max_time = 1.2e-3
gdata.dt_plot =  0.5e-5
gdata.dt_history = 1.0e-5
gdata.max_step = 30000000
gdata.print_count = 20 # Consider changing to 100 for smaller .log files

gdata.cfl = 0.1
gdata.cfl_count = 3
gdata.stringent_cfl = 0 # 1 is more robust
gdata.dt = 1.0e-14
gdata.max_mu_t_factor = 5000.0

# Some additional flags to control what to consider in simulations
include_external_domain = 1
include_nozzles = 1

# Select gas model .. accepting default specific heat ratios and gas constants
include_O2 = 1   # Either 0 (No) or 1 (Yes)
if include_O2 == 0:
    select_gas_model(model='ideal gas', species=['air', 'He'])
elif include_O2 == 1:
    select_gas_model(model='ideal gas', species=['air', 'He', 'O2'])

# ------------------------ Set flow data -------------------------------
# Coflow flow conditions
p_coflow = 101300   # Pa
u_coflow = 487.4    # m/s
T_coflow = 183.8    # degrees K
rho_coflow = p_coflow / (287.1 * T_coflow)  # Uses R from Air

# Centrejet flow conditions
p_cj = 101300   # Pa
u_cj = 1107.3   # m/s
T_cj = 149.05   # degrees K
rho_cj = p_cj / (1537.3 * T_cj)    # Uses R from He-O2

# Coflow flow conditions at throat
p_coflow_th = 306403.0  # Pa
u_coflow_th = 317.0     # m/s
T_coflow_th = 250.0     # degrees K
rho_coflow_th = p_coflow_th / (287.1 * T_coflow_th)  # Uses R from Air

# Centrejet flow conditions at throat
p_cj_th = 313262.0  # Pa
u_cj_th = 759.8     # m/s
T_cj_th = 236.5     # degrees K
rho_cj_th = p_cj_th / (1537.3 * T_cj_th)    # Uses R from He-O2

# Ambient conditions
p_amb = 101300  # Pa
u_amb = 0.0     # m/s
T_amb = 300.0   # degrees K
rho_amb = p_amb / (287.1 * T_amb)  # Uses R from Air

# Estimate turbulence quantities for free stream by specifying turbulence
# intensity and the ratio of turbulent-to-laminar viscosity

# for coflow jet
tke_coflow = 1.5 * (0.01 * u_coflow)**2
mu_t_coflow = 1.0 * sutherland.mu(T_coflow, "Air")
omega_coflow = rho_coflow * tke_coflow / mu_t_coflow

# for centrejet
tke_cj = 0.0  # 1.5 * (0.005 * u_cj)**2
mu_t_cj = 1.0 * sutherland.mu(T_cj, "Air")
omega_cj = rho_cj * tke_cj / mu_t_cj

# for coflow jet at throat
tke_coflow_th = 1.5 * (0.035 * u_coflow_th)**2
mu_t_coflow_th = 5000.0 * sutherland.mu(T_coflow_th, "Air")
omega_coflow_th = rho_coflow_th * tke_coflow_th / mu_t_coflow_th

# for centrejet at throat
tke_cj_th = 1.5 * (0.020 * u_cj_th)**2
mu_t_cj_th = 5000.0 * sutherland.mu(T_cj_th, "Air")
omega_cj_th = rho_cj_th * tke_cj_th / mu_t_cj_th

# for ambient conditions
tke_amb = 1.5 * (0.001 * u_amb)**2
mu_t_amb = 1.0 * sutherland.mu(T_amb, "Air")
omega_amb = rho_amb * tke_amb / mu_t_amb

print "Inflow turbulence (coflow) : tke=", tke_coflow, "omega=", omega_coflow
print "Inflow turbulence (centrejet) : tke=", tke_cj, "omega=", omega_cj
print "Inflow turbulence (coflow throat) : tke=", tke_coflow_th, "omega=", omega_coflow_th
print "Inflow turbulence (centrejet throat) : tke=", tke_cj_th, "omega=", omega_cj_th
print "Inflow turbulence (ambient) : tke=", tke_amb, "omega=", omega_amb

# This set of options is in place because considering a mixture of He-O2 in the
# centrejet causes the simulations to crash. (Need to get Rowan's help)
if include_O2 == 0:
    coflow_cond = FlowCondition(p=p_coflow, u=u_coflow, T=T_coflow,
                                massf=[1.0,0.0], tke=tke_coflow, omega=omega_coflow)
    cj_cond = FlowCondition(p=p_cj, u=u_cj, T=T_cj,
                            massf=[0.0,1.0], tke=tke_cj, omega=omega_cj)
    coflow_cond_th = FlowCondition(p=p_coflow_th, u=u_coflow_th, T=T_coflow_th,
                                   massf=[1.0,0.0], tke=tke_coflow_th, omega=omega_coflow_th)
    cj_cond_th = FlowCondition(p=p_cj_th, u=u_cj_th, T=T_cj_th,
                               massf=[0.0,1.0], tke=tke_cj_th, omega=omega_cj_th)
    amb_cond = FlowCondition(p=p_amb, u=u_amb, T=T_amb,
                             massf=[1.0,0.0], tke=tke_amb, omega=omega_amb)
elif include_O2 == 1:
    coflow_cond = FlowCondition(p=p_coflow, u=u_coflow, T=T_coflow,
                                massf=[1.0,0.0,0.0], tke=tke_coflow, omega=omega_coflow)
    cj_cond = FlowCondition(p=p_cj, u=u_cj, T=T_cj,
                            massf=[0.0,0.7039,0.2961], tke=tke_cj, omega=omega_cj)
    coflow_cond_th = FlowCondition(p=p_coflow_th, u=u_coflow_th, T=T_coflow_th,
                                   massf=[1.0,0.0,0.0], tke=tke_coflow_th, omega=omega_coflow_th)
    cj_cond_th = FlowCondition(p=p_cj_th, u=u_cj_th, T=T_cj_th,
                               massf=[0.0,0.7039,0.2961], tke=tke_cj_th, omega=omega_cj_th)
    amb_cond = FlowCondition(p=p_amb, u=u_amb, T=T_amb,
                             massf=[1.0,0.0,0.0], tke=tke_amb, omega=omega_amb)

# -------------------------- Define geometry ---------------------------------
# Set x-coordinates (for easier handling of geometry)
x0 = 0.0    # indicates x starts at injection plane
x1 = 0.0126556
x2 = 0.265  # Physical domain = 0.261 m, but extended to avoid 
            # having to examine data at the exit boundary
x_ext = -0.246

# Set y-coordinates (for easier handling of geometry)
y0 = 0.0
y1 = 0.005
y2 = 0.00525
y3 = 0.030233
y4 = 0.24425
y5 = 0.33025

# Start building lines with nodes
north0 = Line(Node(x0,y1), Node(x1,y1))
east0 = Line(Node(x1,y0), Node(x1,y1))
south0 = Line(Node(x0,y0), Node(x1,y0))
west0 = Line(Node(x0,y0), Node(x0,y1))

north1 = Line(Node(x0,y2), Node(x1,y2))
east1 = Line(Node(x1,y1), Node(x1,y2))
south1 = north0
west1 = Line(Node(x0,y1), Node(x0,y2))

north2 = Line(Node(x0,y3), Node(x1,y3))
east2 = Line(Node(x1,y2), Node(x1,y3))
south2 = north1
west2 = Line(Node(x0,y2), Node(x0,y3))

north3 = Line(Node(x1,y1), Node(x2,y1))
east3 = Line(Node(x2,y0), Node(x2,y1))
south3 = Line(Node(x1,y0), Node(x2,y0))
west3 = east0

north4 = Line(Node(x1,y2), Node(x2,y2))
east4 = Line(Node(x2,y1), Node(x2,y2))
south4 = north3
west4 = east1

north5 = Line(Node(x1,y3), Node(x2,y3))
east5 = Line(Node(x2,y2), Node(x2,y3))
south5 = north4
west5 = east2

north6 = Line(Node(x1,y5), Node(x2,y5))
east6 = Line(Node(x2,y3), Node(x2,y5))
south6 = north5
west6 = Line(Node(x1,y3), Node(x1,y5))

north7 = Line(Node(x_ext,y5), Node(x1,y5))
east7 = west6
south7 = Line(Node(x_ext,y4), Node(x1,y3))
west7 = Line(Node(x_ext,y4), Node(x_ext,y5))

north8 = Spline2("centrejet_nozzle.dat")
east8 = west0
south8 = Line(Node(north8.eval(0.0).x,y0),east8.eval(0.0))
west8 = Line(south8.eval(0.0),north8.eval(0.0))

north9 = Spline2("coflow_outer_nozzle.dat")
east9 = west2 
south9 = Spline2("coflow_inner_nozzle.dat")
west9 = Line(south9.eval(0.0),north9.eval(0.0))


# ------- Define blocks, boundary conditions and set discretisation ----------
#
#   NOTE: ASCII IMAGE IS NOT TO SCALE, ONLY USE AS GUIDE TO BLOCK DEFINITION
#
#		         _________________________
#		        |    |                    |
#		        | 7  |                    |
#		         \   |         6          |
#		          \  |                    |
#		           \ |                    |
#		   _________\|____________________|
#      	  INFLOW>>|____9__ |2|         5          |
#			 _\|=|====================| <---4 (1 is between 0 and 2)
#		INFLOW>>|_8|0|_________3__________| _ _ _ _ 
#		               SYMMETRY AXIS
#

if include_nozzles == 0:
    blk0 = SuperBlock2D(make_patch(north0, east0, south0, west0),
                        nni=64, nnj=52, nbi=1, nbj=1,
                        fill_condition=cj_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(0,1,1.030),
                                 RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(0,1,1.010)],
                        bc_list=[SlipWallBC(), SlipWallBC(),
                                 SlipWallBC(), SupInBC(cj_cond)])

    blk1 = SuperBlock2D(make_patch(north1, east1, south1, west1),
                        nni=64, nnj=16, nbi=1, nbj=1,
                        fill_condition=coflow_cond,              
                        cf_list=[RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(0,0,0),
                                 RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(1,1,1.09)],
                        bc_list=[SlipWallBC(), SlipWallBC(),
                                 SlipWallBC(), SlipWallBC()])

    blk2 = SuperBlock2D(make_patch(north2, east2, south2, west2),
                        nni=64, nnj=120, nbi=1, nbj=2,
                        fill_condition=coflow_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(1,1,1.020),
                                 RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(1,1,1.005)],
                        bc_list=[AdiabaticBC(), SlipWallBC(),
                                 SlipWallBC(), SupInBC(coflow_cond)])

elif include_nozzles == 1:
    blk0 = SuperBlock2D(make_patch(north0, east0, south0, west0),
                        nni=64, nnj=52, nbi=1, nbj=1,
                        fill_condition=cj_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(0,1,1.030),
                                 RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(0,1,1.010)],
                        bc_list=[SlipWallBC(), ExtrapolateOutBC(),
                                 SlipWallBC(), SlipWallBC()])

    blk1 = SuperBlock2D(make_patch(north1, east1, south1, west1),
                        nni=64, nnj=16, nbi=1, nbj=1,
                        fill_condition=coflow_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(0,0,0),
                                 RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(1,1,1.09)],
                        bc_list=[SlipWallBC(), ExtrapolateOutBC(),
                                 SlipWallBC(), SlipWallBC()])

    blk2 = SuperBlock2D(make_patch(north2, east2, south2, west2),
                        nni=64, nnj=120, nbi=1, nbj=2,
                        fill_condition=coflow_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(1,1,1.020),
                                 RobertsClusterFunction(1,0,1.07),
                                 RobertsClusterFunction(1,1,1.005)],
                        bc_list=[AdiabaticBC(), ExtrapolateOutBC(),
                                 SlipWallBC(), SlipWallBC()])

blk3 = SuperBlock2D(make_patch(north3, east3, south3, west3),
                    nni=240, nnj=52, nbi=6, nbj=1,
                    fill_condition=cj_cond,
                    cf_list=[RobertsClusterFunction(1,0,1.10),
                             RobertsClusterFunction(0,1,1.030),
                             RobertsClusterFunction(1,0,1.10),
                             RobertsClusterFunction(0,1,1.030)],
                    bc_list=[SlipWallBC(), ExtrapolateOutBC(),
                             SlipWallBC(), SlipWallBC()])

blk4 = SuperBlock2D(make_patch(north4, east4, south4, west4),
                    nni=240, nnj=16, nbi=6, nbj=1,
                    fill_condition=coflow_cond,              
                    cf_list=[RobertsClusterFunction(1,0,1.10),
                             RobertsClusterFunction(0,0,0),
                             RobertsClusterFunction(1,0,1.10),
                             RobertsClusterFunction(0,0,0)],
                    bc_list=[SlipWallBC(), ExtrapolateOutBC(),
                             SlipWallBC(), SlipWallBC()])

if include_external_domain == 0:
    blk5 = SuperBlock2D(make_patch(north5, east5, south5, west5),
                        nni=240, nnj=120, nbi=6, nbj=2,
                        fill_condition=coflow_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.10),
                                 RobertsClusterFunction(1,1,1.020),
                                 RobertsClusterFunction(1,0,1.10),
                                 RobertsClusterFunction(1,1,1.020)],
                        bc_list=[ExtrapolateOutBC(), ExtrapolateOutBC(),
                                 SlipWallBC(), SlipWallBC()])

elif include_external_domain == 1:
    blk5 = SuperBlock2D(make_patch(north5, east5, south5, west5),
                        nni=240, nnj=120, nbi=6, nbj=2,
                        fill_condition=coflow_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.10),
                                 RobertsClusterFunction(1,1,1.020),
                                 RobertsClusterFunction(1,0,1.10),
                                 RobertsClusterFunction(1,1,1.020)],
                        bc_list=[SlipWallBC(), ExtrapolateOutBC(),
                                 SlipWallBC(), SlipWallBC()])

    blk6 = SuperBlock2D(make_patch(north6, east6, south6, west6),
                        nni=240, nnj=64, nbi=6, nbj=1,
                        fill_condition=amb_cond,
                        cf_list=[RobertsClusterFunction(1,0,1.10),
                                 RobertsClusterFunction(1,0,1.0006),
                                 RobertsClusterFunction(1,0,1.10),
                                 RobertsClusterFunction(1,0,1.0006)],
                        bc_list=[SlipWallBC(), ExtrapolateOutBC(),
                                 SlipWallBC(), SlipWallBC()])

    blk7 = SuperBlock2D(make_patch(north7, east7, south7, west7),
                        nni=48, nnj=64, nbi=1, nbj=1,
                        fill_condition=amb_cond,
                        cf_list=[RobertsClusterFunction(0,1,1.01),
                                 RobertsClusterFunction(1,0,1.0006),
                                 RobertsClusterFunction(0,1,1.01),
                                 RobertsClusterFunction(1,0,1.0006)],
                        bc_list=[SlipWallBC(), SlipWallBC(),
                                 SlipWallBC(), SlipWallBC()])

if include_nozzles == 1:
    blk8 = SuperBlock2D(make_patch(north8, east8, south8, west8),
                        nni=52, nnj=52, nbi=1, nbj=1,
                        fill_condition=cj_cond,
                        cf_list=[RobertsClusterFunction(1,1,1.05),
                                 RobertsClusterFunction(0,1,1.01),
                                 RobertsClusterFunction(1,1,1.05),
                                 RobertsClusterFunction(0,1,1.01)],
                        bc_list=[AdiabaticBC(), SlipWallBC(),
                                 SlipWallBC(), SupInBC(cj_cond_th)])
    
    blk9 = SuperBlock2D(make_patch(north9, east9, south9, west9),
                        nni=100, nnj=120, nbi=2, nbj=2,
                        fill_condition=coflow_cond,
                        cf_list=[RobertsClusterFunction(1,1,1.02),
                                 RobertsClusterFunction(1,1,1.005),
                                 RobertsClusterFunction(1,1,1.02),
                                 RobertsClusterFunction(1,1,1.005)],
                        bc_list=[AdiabaticBC(), SlipWallBC(),
                                 AdiabaticBC(), SupInBC(coflow_cond_th)])

identify_block_connections()
