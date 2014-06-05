## \file backstep.py
## \brief Turbulent supersonic flow over a backward facing step
## \author Wilson Chan 
## \version 30-May-2008 Original mbcns2 runs
##          29-Dec-2008 Ran only for nozzle part (in a separate
##                      script) and used the exit profile as a
##                      udf inflow to the actual backward-facing
##                      step configuration.
##          05-Jan-2009 Extract the exit profile from the nozzle 
##                      calculations and use user-defined-function.
##                      capability to start off actual calculations.
##          24-Apr-2009 Worked out best turbulence inflow parameters
##                      for the nozzle to match the experimental
##                      results, then re-ran the whole domain again
##                      (with the nozzle). This is because we are
##                      having issues with the udf-supersonic-in.lua
##                      not giving the correct inflow to the 
##                      simulations.
##

from math import *

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


job_title = "Eklund's turbulent supersonic flow over a backward facing step"
print job_title

# -------- Do a little more setting of global data -----------------------------
gdata.turbulence_flag = 1  # Activate turbulence model
gdata.turbulence_model = "k_omega"
gdata.title = job_title
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE

gdata.max_time = 1.2e-3    # should allow a few flow lengths (about 5 here)
gdata.dt_plot =  0.2e-3
gdata.max_step = 30000000

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.stringent_cfl = 0    # 1 is more robust
gdata.dt = 1.0e-14         # only an initial guess, the simulation will take over

gdata.dimensions = 2


# --------------- Define gas model and flow conditions -----------------
# ------------ Accept defaults giving perfect air (gamma=1.4) ---------
select_gas_model(model='ideal gas', species=['air'])

# Define inflow conditions
p_inf = 35.0e3    # Pa
u_inf = 518.0     # m/s
T_inf = 167.0     # degrees K
rho_inf = p_inf / (287.0 * T_inf)  # kg/m**3

# Estimate turbulence quantities for free stream by specifying turbulence
# intensity and the ratio of turbulent-to-laminar viscosity
turbulence_intensity = 0.01
turb_lam_viscosity_ratio = 100.0
tke_inf = 1.5 * (turbulence_intensity * u_inf)**2
mu_t_inf = turb_lam_viscosity_ratio * sutherland_viscosity(T_inf, "Air")
omega_inf = rho_inf * tke_inf / mu_t_inf
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,],
                       tke=tke_inf, omega=omega_inf )


# -------------------------- Define geometry ---------------------------------
x1_noz = -0.00326 
x0 = 0.0
x1 = 0.0345   # Actual experimental domain ... 0.03429 m

y0 = 0.0
y1 = 0.003175  # Step height
y2 = 0.0213    # previous y2 ... 0.0196
y0_noz = -x1_noz * tan(radians(0.5)) + y1  

x0_noz = -0.08  # Start of 'nozzle'
y1_noz = -x0_noz * tan(radians(0.5)) + y1

north0 = Line(Node(x0,y1), Node(x1,y1))
east0 = Line(Node(x1,y0), Node(x1,y1))                          
south0 = Line(Node(x0,y0), Node(x1,y0))
west0 = Line(Node(x0,y0), Node(x0,y1))

north1 = Line(Node(x0,y2), Node(x1,y2))
east1 = Line(Node(x1,y1), Node(x1,y2))                          
south1 = north0
west1 = Line(Node(x0,y1), Node(x0,y2))

north2 = Line(Node(x1_noz,y2), Node(x0,y2))
east2 = west1                          
south2 = Line(Node(x1_noz,y0_noz), Node(x0,y1))
west2 = Line(Node(x1_noz,y0_noz), Node(x1_noz,y2))

north_noz = Line(Node(x0_noz,y2), Node(x1_noz,y2))
east_noz = Line(Node(x1_noz,y0_noz), Node(x1_noz,y2))
south_noz = Line(Node(x0_noz,y1_noz), Node(x1_noz,y0_noz))
west_noz = Line(Node(x0_noz,y1_noz), Node(x0_noz,y2))


# ------- Define blocks, boundary conditions and set discretisation ---------
nB_I = 5  # Number of i-blocks for SuperBlock2D
nB_J = 1  # Number of j-blocks for SuperBlock2D
blk0 = SuperBlock2D(make_patch(north0, east0, south0, west0),
                    nni=140, nnj=60, nbi=nB_I, nbj=nB_J,
                    fill_condition=inflow,
                    cf_list=[RobertsClusterFunction(1,0,1.08),
                             RobertsClusterFunction(1,1,1.01),
                             RobertsClusterFunction(1,0,1.08),
                             RobertsClusterFunction(1,1,1.01)],
                    bc_list=[SlipWallBC(), ExtrapolateOutBC(),
                             AdiabaticBC(), AdiabaticBC()]) 
                             # WEST can be SlipWallBC()

nB_I = 5  # Number of i-blocks for SuperBlock2D
nB_J = 2  # Number of j-blocks for SuperBlock2D
blk1 = SuperBlock2D(make_patch(north1, east1, south1, west1),
                    nni=140, nnj=120, nbi=nB_I, nbj=nB_J,
                    fill_condition=inflow,              
                    cf_list=[RobertsClusterFunction(1,0,1.08),
                             RobertsClusterFunction(1,1,1.0025),
                             RobertsClusterFunction(1,0,1.08),
                             RobertsClusterFunction(1,1,1.0025)],
                    bc_list=[AdiabaticBC(), ExtrapolateOutBC(),
                             SlipWallBC(), SlipWallBC()])

nB_I = 1  # Number of i-blocks for SuperBlock2D
nB_J = 2  # Number of j-blocks for SuperBlock2D
blk2 = SuperBlock2D(make_patch(north2, east2, south2, west2),
                    nni=24, nnj=120, nbi=nB_I, nbj=nB_J,
                    fill_condition=inflow,
                    cf_list=[RobertsClusterFunction(0,1,1.18),  
                             RobertsClusterFunction(1,1,1.0025),  
                             RobertsClusterFunction(0,1,1.18),  
                             RobertsClusterFunction(1,1,1.0025)],
                    bc_list=[AdiabaticBC(), SlipWallBC(), 
                             AdiabaticBC(), SlipWallBC()]) 

nB_I = 5   # Number of i-blocks for SuperBlock2D
nB_J = 2   # Number of j-blocks for SuperBlock2D
blk_noz = SuperBlock2D(make_patch(north_noz, east_noz, 
                                  south_noz, west_noz),
                       nni=100, nnj=120, nbi=nB_I, nbj=nB_J,
                       fill_condition=inflow,
                       cf_list=[RobertsClusterFunction(1,1,1.06),
                                RobertsClusterFunction(1,1,1.0025),
                                RobertsClusterFunction(1,1,1.06),
                                RobertsClusterFunction(1,1,1.0025)],
                       bc_list=[AdiabaticBC(), SlipWallBC(),
                                AdiabaticBC(), SupInBC(inflow)])

identify_block_connections()
