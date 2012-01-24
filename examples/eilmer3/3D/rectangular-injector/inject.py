# inject.py -- 16:1 rectagular fuel injector.
# PJ, VW
# Eilmer3 port: 06-Feb-2010

# -------------------------------------------------------------------
# Some handy definitions for later.
import math
# from cfpylib.geom.box3d import makeSimpleBox

# ---------------- First, set the global data ----------------------
gdata.title = "16:1 rectagular fuel injector."
gdata.dimensions = 3
gdata.dt = 4.0e-10
gdata.t_order = 1
gdata.max_time = 4.0e-5
gdata.max_step = 80000
gdata.reacting_flag = 0
gdata.dt_plot = 0.25e-5
gdata.dt_history = 1.0e-7

gdata.turbulence_flag = 1  
gdata.turbulence_model = "k_omega"  
gdata.viscous_flag = 1
gdata.diffusion_model = "FicksFirstLaw"
gdata.diffusion_flag = 1
gdata.cfl = 0.4	
gdata.cfl_count = 3


# ------------ Second, set up flow conditions -------------------
# These will be used for fill and boundary conditions.
species_list = select_gas_model(model='ideal gas', 
                                species=['H2', 'air'])

# Define free stream conditions
p_inf = 13.861e3 # Pa
u_inf = 2780.0   # m/s
T_inf = 770.0    # degrees K

# Estimate turbulence quantities for free stream
# by specifying the intensity as 0.05 and estimating the
# turbulence viscosity as 0.1 times the laminar viscosity
# to match Viti et al.
tke_inf = 1.5 * (u_inf * 0.05)**2
rho_inf = p_inf / (287.0 * T_inf)
def mu_air(T):
    "Sutherland expression for air viscosity."
    from math import sqrt
    mu_ref = 17.89e-6; T_ref = 273.1; S = 110.4
    T_T0 = T / T_ref
    return mu_ref * (T_ref + S)/(T + S) * T_T0 * sqrt(T_T0);
mu_t_inf = 0.1 * mu_air(T_inf)
omega_inf = rho_inf * tke_inf / mu_t_inf
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

# Estimate turbulence quantities for H2 inflow turbulence
# by specifying the intensity and length scale for pipe flow
L2 = sqrt(2.0)*4.0e-3           # streamwise length of injector
pH2 = 528.28e3			# Pa
uH2 = 842.42			# m/s
wH2 = 842.42			# m/s
TH2 = 245.833			# K
rhoH2 = pH2/(4121*TH2)		# kg/m**3
mu_H2 = 1.4741e-5 		# Pa.s, NIST chemical webbook at 1atm, 625K
ReH2 = rhoH2*sqrt(2)*uH2*L2/mu_H2 # Reynolds number based on jet width
lt = 0.07*L2			# Turbulence length scale, pipe flow approx
It  = 0.16*ReH2**(-0.125) 	# Turbulence intensity
tke_H2 = 1.5*(sqrt(2)*uH2*It)**2
omega_H2 = 0.09**(-0.25)*sqrt(tke_H2)/lt
print "Fuel inflow turbulence: tke=", tke_H2, "omega=", omega_H2

inflowCond = FlowCondition(p=p_inf, u=u_inf, v=0.0, w=0.0, T=T_inf, massf={'air':1.0},tke=tke_inf, omega=omega_inf)
injectCond = FlowCondition(p=pH2, u=uH2, v=0.0, w=wH2, T=TH2, massf={'H2':1.0}, tke=tke_H2, omega=omega_H2)

# ------------ Third, set up the blocks ---------------------
# Parameters defining the duct...
Lmix = 33.0e-3    	# distance from jet centroid to domain end
L1 = 10.0e-3     	# distance from leading edge to injector
Whalf0 = 12.0e-3 	# half-width of duct
Whalf1 = 1.25e-4 	# half-width of injector
H = 16.0e-3      	# height of duct
L0 = L1+Lmix    	# length of duct in flow direction

# Plan of blocks
#                 NORTH BNDRY
#         +--------+---+---------------+
#         |........|.:.|....:....:.....|
#         |...01...|.11|....:.21.:.....|
# inflow> |........|.:.|....:....:.....| outflow>
# (WEST)  +--------+---+---------------+ (EAST)
#         |   00   | 10|    : 20 :     |
#         +--------+---+---------------+
#                 SOUTH BNDRY
#
#                    ^
#                 injector
cluster_k  = RobertsClusterFunction(1, 0, 1.025) # cluster down, toward the bottom surface
cluster_kend = RobertsClusterFunction(1, 0, 1.15) # cluster down, toward the bottom surface
cluster_i0 = RobertsClusterFunction(0, 1, 1.1) # cluster streamwise toward injector
cluster_i2 = RobertsClusterFunction(1, 0, 1.125)
cluster_j1 = RobertsClusterFunction(1, 0, 1.015) # cluster cross-stream toward injector

ni0 = 40; ni1 = 100; ni2 = 150;
nj0 = 10; nj1 = 72;
nk0 = 100;

# upstream pair of blocks
pv = makeSimpleBox(xPos=0.0, yPos=0.0, xSize=L1, ySize=Whalf1, zSize=H)
cflist = [cluster_i0,None,cluster_i0,None]*2 + [cluster_k,]*4;
# list of boundary conditions in N, E, S, W, T, B order:
bclist = [AdjacentBC(), AdjacentBC(), SlipWallBC(),
          UserDefinedBC("udf-supersonic-in.lua"),
	  SlipWallBC(), AdiabaticBC()]
# 12 edges is a full complement; see elmer_prep.py for the order of edges
blk00 = Block3D(parametric_volume=pv, cf_list=cflist,
                fill_condition=inflowCond,
                nni=ni0, nnj=nj0, nnk=nk0, 
		bc_list=bclist)
	       
pv = makeSimpleBox(xPos=0.0,yPos=Whalf1, xSize=L1,ySize=Whalf0-Whalf1, zSize=H)
cflist = [cluster_i0,cluster_j1,cluster_i0,cluster_j1]*2 + [cluster_k,]*4;
bclist = [SlipWallBC(),AdjacentBC(), AdjacentBC(),
          UserDefinedBC("udf-supersonic-in.lua"),
	  SlipWallBC(), AdiabaticBC()]
blk01 = SuperBlock3D(parametric_volume=pv, cf_list=cflist,
                     fill_condition=inflowCond,
                     nni=ni0, nnj=nj1, nnk=nk0,
		     nbi=1, nbj=3, nbk=1,
		     bc_list=bclist)

# injector and part of plate beside it
pv = makeSimpleBox(xPos=L1,yPos=0.0, xSize=L2, ySize=Whalf1, zSize=H)
cflist = [None,None,None,None]*2 + [cluster_k,]*4;
bclist = [AdjacentBC(), AdjacentBC(), SlipWallBC(), AdjacentBC(),
	  SlipWallBC(), SupInBC(injectCond)]
blk10 = SuperBlock3D(parametric_volume=pv, cf_list=cflist,
                     fill_condition=inflowCond,
                     nni=ni1, nnj=nj0, nnk=nk0,
		     nbi=2, nbj=1, nbk=1,
		     bc_list=bclist)

pv = makeSimpleBox(xPos=L1,yPos=Whalf1, xSize=L2, ySize=Whalf0-Whalf1, zSize=H)
cflist = [None,cluster_j1,None,cluster_j1]*2 + [cluster_k,]*4;
bclist = [SlipWallBC(),AdjacentBC(), AdjacentBC(), AdjacentBC(),
	  SlipWallBC(), AdiabaticBC()]
blk11 = SuperBlock3D(parametric_volume=pv, cf_list=cflist,
                     fill_condition=inflowCond,
                     nni=ni1, nnj=nj1, nnk=nk0,
		     nbi=2, nbj=3, nbk=1,
		     bc_list=bclist)

# blocks downstream of injector
pv = makeSimpleBox(xPos=L1+L2,yPos=0.0, xSize=L0-(L1+L2), ySize=Whalf1, zSize=H)
cflist = [cluster_i2,None,cluster_i2,None]*2 + [cluster_k,cluster_kend,cluster_kend,cluster_k];
bclist = [AdjacentBC(), ExtrapolateOutBC(sponge_flag=0), SlipWallBC(), AdjacentBC(),
	  SlipWallBC(), AdiabaticBC()]
blk20 = SuperBlock3D(parametric_volume=pv, cf_list=cflist,
                     fill_condition=inflowCond,
                     nni=ni2, nnj=nj0, nnk=nk0,
		     nbi=3, nbj=1, nbk=1,
		     bc_list=bclist, hcell_list=[(ni2-1,0,0.5*nk0)])

pv = makeSimpleBox(xPos=L1+L2,yPos=Whalf1, xSize=L0-(L1+L2), ySize=Whalf0-Whalf1, zSize=H)
cflist = [cluster_i2,cluster_j1,cluster_i2,cluster_j1]*2 + [cluster_k,cluster_kend,cluster_kend,cluster_k];
bclist = [SlipWallBC(), ExtrapolateOutBC(sponge_flag=0), AdjacentBC(), AdjacentBC(),
	  SlipWallBC(), AdiabaticBC()]
blk21 = SuperBlock3D(parametric_volume=pv, cf_list=cflist,
                     fill_condition=inflowCond,
                     nni=ni2, nnj=nj1, nnk=nk0,
		     nbi=3, nbj=3, nbk=1,
		     bc_list=bclist)

identify_block_connections()
