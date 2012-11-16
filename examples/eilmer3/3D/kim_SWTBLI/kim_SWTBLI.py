## \file kim_SWTBLI.py     
## \brief 3D computations of the single-fin shock-wave-turbulent-
##        boundary-layer-interaction experiments by Kim et. al.
##        (AIAA/J 1991 v29 n10)
##
## \author Wilson Chan 21 June 2010
## \version1 Full domain with extra blocks used instead of
##           using ExtrapolateOutBC
## \version2 Just use the ExtrapolateOutBC, so only using 
##           2 blocks now 
## \version3 Will shorten the first block and just use the udf.lua
##           technique to take in the inflow - this should save us
##           some compute time.
## \version4 Will move away from the simplistic 2-block grid and
##           use a radiating grid like how the experts do it. This
##           should hopefully help us when doing the postprocessing
##           as well, when we try to extract compare results at beta


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

# Set up global data
job_title = "Single-fin shock-wave-turbulent-boundary-layer-interaction M=4.0 D=20deg"
print job_title
gdata.dimensions = 3
gdata.turbulence_flag = 1 
gdata.turbulence_model = "k_omega"
gdata.title = job_title
gdata.viscous_flag = 1
gdata.flux_calc = AUSMDV  # ADAPTIVE
gdata.print_count = 100
        
gdata.max_time = 1.0e-4 #8.0e-4  # 1 flow length ~ 2e-4 secs ; 4 flow lengths ~ 8e-4 secs  
gdata.dt_plot =  0.5e-4 
gdata.dt_history = 1.0e-3
gdata.max_step = 30000000

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.stringent_cfl = 0  # 1 is more robust
gdata.dt = 1.0e-10  # Only an initial guess - the simulation will take over

# Select gas model - accepting defaults for air (R=287.1, gamma=1.4)
select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions this to match the stagnation conditions 
# provided by Kim et. al.
p_inf = 10.3e3  # Pa
u_inf = 707.0   # m/s
T_inf = 70.4    # degrees K
rho_inf = p_inf / (287.1 * T_inf)
# From pitot survey conducted by Lu, F K with a condition quite
# similar to Kim's
#p_inf = 10.483e3  # Pa
#u_inf = 647.34    # m/s
#T_inf = 69.4      # degrees K
#rho_inf = 0.523   # kg/m^3 #p_inf / (287.1 * T_inf)

# Estimate turbulence quantities for free stream by specifying turbulence
# intensity and the ratio of turbulent-to-laminar viscosity
turbulence_intensity = 0.01
turb_lam_viscosity_ratio = 1.0
tke_inf = 1.5 * (turbulence_intensity * u_inf)**2
mu_t_inf = turb_lam_viscosity_ratio * sutherland_viscosity(T_inf, "Air")
omega_inf = rho_inf * tke_inf / mu_t_inf
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,],
                       tke=tke_inf, omega=omega_inf)

# It seems that everyone has been starting the flow field with inflow 
# conditions throughout and then letting the boundary layer grow out 
# into the main stream, instead of having an initial state that is
# more representative of the actual tunnel conditions (which
# unfortunately seems to wreck havoc with our turbulence model).
initial_flag = "inflow"   # "inflow" or "actual"
if initial_flag == "inflow":
    initial = inflow
elif initial_flag == "actual":
    initial = FlowCondition(p= 0.1*p_inf, u=0.0, v=0.0, T=296.0, massf=[1.0,],
                            tke=tke_inf/100.0, omega=omega_inf/10.0)


# Fin dimensions
finL = 0.14  # Length of fin (extended from 0.127m to 0.14m)
finH = 0.076 # Height of fin


# Build nodes
inflowBlk1_p0 = Node(-0.01, 0.0, -0.01*sin(radians(45.0)))
inflowBlk1_p1 = Node(-0.01+0.01*cos(radians(45.0)), 0.0, -0.01*sin(radians(45.0)))
inflowBlk1_p2 = Node(-0.01+0.01*cos(radians(45.0)), finH, -0.01*sin(radians(45.0))) 
inflowBlk1_p3 = Node(-0.01, finH, -0.01*sin(radians(45.0))) 
inflowBlk1_p4 = Node(-0.01, 0.0, 0.0)
inflowBlk1_p5 = Node(0.0, 0.0, 0.0)
inflowBlk1_p6 = Node(0.0, finH, 0.0)
inflowBlk1_p7 = Node(-0.01, finH, 0.0)

inflowBlk2_p0 = Node(-0.01, 0.0, -1*finL*sin(radians(60.0)))
inflowBlk2_p1 = Node(finL*cos(radians(60.0)), 0.0, -1*finL*sin(radians(60.0)))
inflowBlk2_p2 = Node(finL*cos(radians(60.0)), finH, -1*finL*sin(radians(60.0)))
inflowBlk2_p3 = Node(-0.01, finH, -1*finL*sin(radians(60.0)))
inflowBlk2_p4 = inflowBlk1_p0
inflowBlk2_p5 = inflowBlk1_p1
inflowBlk2_p6 = inflowBlk1_p2
inflowBlk2_p7 = inflowBlk1_p3

mainBlk_p0 = inflowBlk2_p1      
mainBlk_p1 = Node(finL*cos(radians(20.0)), 0.0, -1*finL*sin(radians(20.0)))
mainBlk_p2 = Node(finL*cos(radians(20.0)), finH, -1*finL*sin(radians(20.0)))
mainBlk_p3 = inflowBlk2_p2      
mainBlk_p4 = inflowBlk1_p1      
mainBlk_p5 = inflowBlk1_p5      
mainBlk_p6 = inflowBlk1_p6      
mainBlk_p7 = inflowBlk1_p2      


# Build paths
inflowBlk1_p01 = Line(inflowBlk1_p0, inflowBlk1_p1)
inflowBlk1_p12 = Line(inflowBlk1_p1, inflowBlk1_p2)
inflowBlk1_p32 = Line(inflowBlk1_p3, inflowBlk1_p2)
inflowBlk1_p03 = Line(inflowBlk1_p0, inflowBlk1_p3)
inflowBlk1_p45 = Line(inflowBlk1_p4, inflowBlk1_p5)
inflowBlk1_p56 = Line(inflowBlk1_p5, inflowBlk1_p6)
inflowBlk1_p76 = Line(inflowBlk1_p7, inflowBlk1_p6)
inflowBlk1_p47 = Line(inflowBlk1_p4, inflowBlk1_p7)
inflowBlk1_p04 = Line(inflowBlk1_p0, inflowBlk1_p4)
inflowBlk1_p15 = Arc(inflowBlk1_p1, inflowBlk1_p5, inflowBlk1_p4)
inflowBlk1_p26 = Arc(inflowBlk1_p2, inflowBlk1_p6, inflowBlk1_p7)
inflowBlk1_p37 = Line(inflowBlk1_p3, inflowBlk1_p7)

inflowBlk2_p01 = Line(inflowBlk2_p0, inflowBlk2_p1)
inflowBlk2_p12 = Line(inflowBlk2_p1, inflowBlk2_p2)
inflowBlk2_p32 = Line(inflowBlk2_p3, inflowBlk2_p2)
inflowBlk2_p03 = Line(inflowBlk2_p0, inflowBlk2_p3)
inflowBlk2_p45 = Line(inflowBlk2_p4, inflowBlk2_p5)
inflowBlk2_p56 = Line(inflowBlk2_p5, inflowBlk2_p6)
inflowBlk2_p76 = Line(inflowBlk2_p7, inflowBlk2_p6)
inflowBlk2_p47 = Line(inflowBlk2_p4, inflowBlk2_p7)
inflowBlk2_p04 = Line(inflowBlk2_p0, inflowBlk2_p4)
inflowBlk2_p15 = Line(inflowBlk2_p1, inflowBlk2_p5)
inflowBlk2_p26 = Line(inflowBlk2_p2, inflowBlk2_p6)
inflowBlk2_p37 = Line(inflowBlk2_p3, inflowBlk2_p7)

mainBlk_p01 = Arc(mainBlk_p0, mainBlk_p1, mainBlk_p5)
mainBlk_p12 = Line(mainBlk_p1, mainBlk_p2)
mainBlk_p32 = Arc(mainBlk_p3, mainBlk_p2, mainBlk_p6)
mainBlk_p03 = Line(mainBlk_p0, mainBlk_p3)
mainBlk_p45 = inflowBlk1_p15 
mainBlk_p56 = Line(mainBlk_p5, mainBlk_p6)
mainBlk_p76 = inflowBlk1_p26
mainBlk_p47 = Line(mainBlk_p4, mainBlk_p7)
mainBlk_p04 = Line(mainBlk_p0, mainBlk_p4)
mainBlk_p15 = Line(mainBlk_p1, mainBlk_p5)
mainBlk_p26 = Line(mainBlk_p2, mainBlk_p6)
mainBlk_p37 = Line(mainBlk_p3, mainBlk_p7)


# Build parametric volumes
inflowBlk1_pvol = WireFrameVolume(inflowBlk1_p01, inflowBlk1_p12,
                                  inflowBlk1_p32, inflowBlk1_p03,
                                  inflowBlk1_p45, inflowBlk1_p56, 
                                  inflowBlk1_p76, inflowBlk1_p47, 
                                  inflowBlk1_p04, inflowBlk1_p15, 
                                  inflowBlk1_p26, inflowBlk1_p37)

inflowBlk2_pvol = WireFrameVolume(inflowBlk2_p01, inflowBlk2_p12,
                                  inflowBlk2_p32, inflowBlk2_p03,
                                  inflowBlk2_p45, inflowBlk2_p56, 
                                  inflowBlk2_p76, inflowBlk2_p47, 
                                  inflowBlk2_p04, inflowBlk2_p15, 
                                  inflowBlk2_p26, inflowBlk2_p37)

mainBlk_pvol = WireFrameVolume(mainBlk_p01, mainBlk_p12,
                               mainBlk_p32, mainBlk_p03,
                               mainBlk_p45, mainBlk_p56,
                               mainBlk_p76, mainBlk_p47,
                               mainBlk_p04, mainBlk_p15,
                               mainBlk_p26, mainBlk_p37)

# Value given to Robert's cluster function for all walls
wallClustVal = 1.050

# Define the blocks, boundary conditions and set the discretisation.
inflowBlk1 = SuperBlock3D(parametric_volume=inflowBlk1_pvol,
                          nni=15, nnj=54, nnk=54,
                          nbi=1, nbj=2, nbk=2,
                          fill_condition=initial, label="inflowBlk1",
                          cf_list=[RobertsClusterFunction(0,0, 1.00),    # p01
                                   RobertsClusterFunction(1,0, wallClustVal),   # p12
                                   RobertsClusterFunction(0,0, 1.00),    # p32
                                   RobertsClusterFunction(1,0, wallClustVal),   # p03
                                   RobertsClusterFunction(0,0, 1.00),    # p45   
                                   RobertsClusterFunction(1,0, wallClustVal),   # p56
                                   RobertsClusterFunction(0,0, 1.00),    # p76
                                   RobertsClusterFunction(1,0, wallClustVal),   # p47
                                   RobertsClusterFunction(0,1, wallClustVal),   # p04
                                   RobertsClusterFunction(0,1, wallClustVal),   # p15
                                   RobertsClusterFunction(0,1, wallClustVal),   # p26
                                   RobertsClusterFunction(0,1, wallClustVal)],  # p37
                          bc_list=[SlipWallBC(), AdjacentBC(), AdiabaticBC(), # NESWTB
                                   UserDefinedBC("udf-supersonic-in.lua"), # SupInBC(inflow), 
                                   SlipWallBC(), AdjacentBC()])
      
inflowBlk2 = SuperBlock3D(parametric_volume=inflowBlk2_pvol,
                          nni=15, nnj=54, nnk=54,
                          nbi=1, nbj=2, nbk=2,
                          fill_condition=initial, label="inflowBlk2",
                          cf_list=[RobertsClusterFunction(0,0, 1.00),    # p01
                                   RobertsClusterFunction(1,0, wallClustVal),   # p12
                                   RobertsClusterFunction(0,0, 1.00),    # p32
                                   RobertsClusterFunction(1,0, wallClustVal),   # p03
                                   RobertsClusterFunction(0,0, 1.00),    # p45   
                                   RobertsClusterFunction(1,0, wallClustVal),   # p56
                                   RobertsClusterFunction(0,0, 1.00),    # p76
                                   RobertsClusterFunction(1,0, wallClustVal),   # p47
                                   RobertsClusterFunction(0,1, 1.050),   # p04
                                   RobertsClusterFunction(0,1, 1.050),   # p15
                                   RobertsClusterFunction(0,1, 1.050),   # p26
                                   RobertsClusterFunction(0,1, 1.050)],  # p37
                          bc_list=[SlipWallBC(), AdjacentBC(), AdiabaticBC(), # NESWTB
                                   UserDefinedBC("udf-supersonic-in.lua"), # SupInBC(inflow), 
                                   AdjacentBC(), SlipWallBC()])

mainBlk = SuperBlock3D(parametric_volume=mainBlk_pvol,
                       nni=54, nnj=54, nnk=54,
                       nbi=2, nbj=2, nbk=2,
                       fill_condition=initial, label="mainBlk",
                       cf_list=[RobertsClusterFunction(0,1, wallClustVal-0.025),   # p01
                                RobertsClusterFunction(1,0, wallClustVal),   # p12
                                RobertsClusterFunction(0,1, wallClustVal-0.025),   # p32
                                RobertsClusterFunction(1,0, wallClustVal),   # p03
                                RobertsClusterFunction(0,1, wallClustVal),   # p45   
                                RobertsClusterFunction(1,0, wallClustVal),   # p56
                                RobertsClusterFunction(0,1, wallClustVal),   # p76
                                RobertsClusterFunction(1,0, wallClustVal),   # p47
                                RobertsClusterFunction(0,1, 1.050),   # p04
                                RobertsClusterFunction(0,1, 1.050),   # p15
                                RobertsClusterFunction(0,1, 1.050),   # p26
                                RobertsClusterFunction(0,1, 1.050)],  # p37
                       bc_list=[SlipWallBC(), AdiabaticBC(), AdiabaticBC(), # NESWTB
                                AdjacentBC(), AdjacentBC(), ExtrapolateOutBC()])

identify_block_connections()
