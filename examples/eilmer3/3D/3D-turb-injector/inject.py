# injector.py -- single discrete-hole injection.
# PJ
# Elmer2 original: Nov-2006
# Eilmer3 port: 06-Feb-2010
# Modified 24/9/14 by Samuel Stennett for Undergraduate Thesis Work

# -------------------------------------------------------------------
# Some handy definitions for later.
import math
import numpy as np
from cfpylib.gasdyn import sutherland

# ---------------- First, set the global data ----------------------
gdata.title = "Test Case 3 - 3D Flat Plate Injector."
gdata.dimensions = 3
gdata.viscous_flag = 1
gdata.turbulence_model = "k_omega"
gdata.flux_calc = AUSMDV

gdata.max_time = 1.2e-3 # About 3 flow lengths (1 flow length ~ 0.41 ms) 
gdata.max_step = 30000000
gdata.dt_plot = 2.5e-5			#CONSIDER CHANGING TO GET MORE DATA ENTRIES
gdata.dt_history = 1.0e-5

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-9  # Only an initial guess - the simulation will take over

# Estimate turbulence quantities for free stream by specifying turbulence
# intensity and the ratio of turbulent-to-laminar viscosity.
u_inf = 671.92
T_inf = 70.3
rho_inf = 7.1e3 / (287.1 * 70.3)
turbulence_intensity = 0.05 #assume 5%, from paper
turb_lam_viscosity_ratio = 10.0 #assume turb is 1/10th lam, from paper
tke_inf = 1.5 * (turbulence_intensity * u_inf)**2
mu_t_inf = turb_lam_viscosity_ratio * sutherland.mu(T_inf, "Air")
omega_inf = rho_inf * tke_inf / mu_t_inf
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

# ------------ Second, set up flow conditions -------------------
# These will be used for fill and boundary conditions.
species_list = select_gas_model(model='ideal gas', species=['air'])
initialCond = FlowCondition(p=7.1e3, u=671.92, T=70.3, massf=[1.0,], tke=tke_inf, omega=omega_inf)
injectCond = FlowCondition(p=2006.0e3, w=323.67, T=261.0, massf=[1.0,], tke=tke_inf, omega=omega_inf)
 
# ------------ Third, set up the blocks ---------------------
# Parameters defining the geometry (SEE PICTURE IN THESIS BOOK FOR LABELS)
tlen = 0.2769
twidth = 0.1143
theight = 0.1524
L0 = 0.0762 				#Leading distance to front of injector
rad1 = 0.00412/2.0 			#Injector radius
L1 = tlen - (L0 + 2.0*rad1)		#Trailing distance from back of injector
L2 = L0 + rad1 - rad1*np.cos(np.pi/4.0)	#Leading distance, hitting side of injector 		MAY NEED TO BE ROUNDED round(L2,4)
L3 = 2.0*rad1*np.cos(np.pi/4.0)		#Distance between side intersections of injector 	AS ABOVE
L4 = L1 + rad1 - rad1*np.cos(np.pi/4.0)	#Trailing distance, from the injector intersection 	AS ABOVE
W0 = rad1*np.sin(np.pi/4.0)		#Width to intersection of injector			AS ABOVE
W1 = twidth - W0			#Width from intersection to boundary			AS ABOVE
#Hall = theight/3.0			#Height between blocks to height boundary		
	
# Plan of blocks
#                 NORTH BNDRY
#         +--------+---+---------------+
#         |        |   |               |
#         |   01   | 11|      21       |
# inflow> |        |   |               | outflow>
# (WEST)  +--------+/-\+---------------+ (EAST)
#         |   00   / 10\      20       |
#         +-------+--+--+--------------+
#                 SOUTH BNDRY
#                    ^
#                 injector
#

#ARC CENTRE POINTS
jetcentre0 = Node(L0+rad1,0,0)
jetcentre1 = Node(L0+rad1,0,theight)

#BOTTOM SURFACE POINTS
p000 = Node(0,0,0)
p001 = Node(L0,0,0)
p002 = Node(L0+2.0*rad1,0,0)
p003 = Node(tlen,0,0)
p010 = Node(0,W0,0)
p011 = Node(L2,W0,0)		#MAY NEED TO BE ROUNDED
p012 = Node(L2+L3,W0,0)		#MAY NEED TO BE ROUNDED
p013 = Node(tlen,W0,0)		#MAY NEED TO BE ROUNDED
p020 = Node(0,twidth,0)
p021 = Node(L2,twidth,0)
p022 = Node(L2+L3,twidth,0)
p023 = Node(tlen,twidth,0)	

#TOP SURFACE POINTS
p100 = Node(0,0,theight)
p101 = Node(L0,0,theight)
p102 = Node(L0+2.0*rad1,0,theight)
p103 = Node(tlen,0,theight)
p110 = Node(0,W0,theight)
p111 = Node(L2,W0,theight)		#MAY NEED TO BE ROUNDED
p112 = Node(L2+L3,W0,theight)	#MAY NEED TO BE ROUNDED
p113 = Node(tlen,W0,theight)	#MAY NEED TO BE ROUNDED
p120 = Node(0,twidth,theight)
p121 = Node(L2,twidth,theight)
p122 = Node(L2+L3,twidth,theight)
p123 = Node(tlen,twidth,theight)

#BOTTOM SURFACE LINES
Line000 = Line(p000,p001)
Line001 = Line(p001,p002)
Line002 = Line(p002,p003)
Line003 = Line(p000,p010)
Line004 = Arc(p001,p011,jetcentre0)
Line005 = Arc(p002,p012,jetcentre0)
Line006 = Line(p003,p013)
Line007 = Line(p010,p011)
Line008 = Arc(p011,p012,jetcentre0)
Line009 = Line(p012,p013)
Line010 = Line(p010,p020)
Line011 = Line(p011,p021)
Line012 = Line(p012,p022)
Line013 = Line(p013,p023)
Line014 = Line(p020,p021)
Line015 = Line(p021,p022)
Line016 = Line(p022,p023)

#TOP SURFACE LINES
Line100 = Line(p100,p101)
Line101 = Line(p101,p102)
Line102 = Line(p102,p103)
Line103 = Line(p100,p110)
Line104 = Arc(p101,p111,jetcentre1)
Line105 = Arc(p102,p112,jetcentre1)
Line106 = Line(p103,p113)
Line107 = Line(p110,p111)
Line108 = Arc(p111,p112,jetcentre1)
Line109 = Line(p112,p113)
Line110 = Line(p110,p120)
Line111 = Line(p111,p121)
Line112 = Line(p112,p122)
Line113 = Line(p113,p123)
Line114 = Line(p120,p121)
Line115 = Line(p121,p122)
Line116 = Line(p122,p123)

#VERTICAL LINES
Vert000 = Line(p000,p100)
Vert001 = Line(p001,p101)
Vert002 = Line(p002,p102)
Vert003 = Line(p003,p103)
Vert010 = Line(p010,p110)
Vert011 = Line(p011,p111)
Vert012 = Line(p012,p112)
Vert013 = Line(p013,p113)
Vert020 = Line(p020,p120)
Vert021 = Line(p021,p121)
Vert022 = Line(p022,p122)
Vert023 = Line(p023,p123)


#NEED TO FIX NUMBER OF CELLS AND CLUSTERING


# upstream pair of blocks
ni1 = 25 #76
ni2 = 5 #5
ni3 = 50 #196
nj1 = 5 #2
nj2 = 70 #112
nk1 = 80

#BLOCK 00
cluster_k  = RobertsClusterFunction(1, 0, 1.0008)
cluster_i = RobertsClusterFunction(0, 1, 1.1)
#cluster_j = RobertsClusterFunction(1, 0, 1.0)

cflist = [cluster_i,None,cluster_i,None,]*2 + [cluster_k,]*4;
bclist = [AdjacentBC(),AdjacentBC(),SlipWallBC(),StaticProfBC(filename="REPLACE.dat",n_profile=2),ExtrapolateOutBC(),AdiabaticBC()];
blk00 = SuperBlock3D(nni=ni1, nnj=nj1, nnk=nk1, nbi=1, nbj=1, nbk=8, parametric_volume=WireFrameVolume(Line000,Line004,Line007,Line003,Line100,Line104,Line107,Line103,Vert000,Vert001,Vert011,Vert010),
                cf_list=cflist, bc_list=bclist, fill_condition=initialCond)

#BLOCK 01
cluster_k  = RobertsClusterFunction(1, 0, 1.0008)
cluster_i = RobertsClusterFunction(0, 1, 1.1)
cluster_j = RobertsClusterFunction(1, 0, 1.1)

cflist = [cluster_i,cluster_j,cluster_i,cluster_j,]*2 + [cluster_k,]*4;
bclist = [ExtrapolateOutBC(),AdjacentBC(),AdjacentBC(),StaticProfBC(filename="REPLACE.dat",n_profile=2),ExtrapolateOutBC(),AdiabaticBC()];
blk01 = SuperBlock3D(nni=ni1, nnj=nj2, nnk=nk1, nbi=1, nbj=2, nbk=8, parametric_volume=WireFrameVolume(Line007,Line011,Line014,Line010,Line107,Line111,Line114,Line110,Vert010,Vert011,Vert021,Vert020),
                cf_list=cflist, bc_list=bclist, fill_condition=initialCond)

# injector and part of plate beside it

#BLOCK 10
cluster_k  = RobertsClusterFunction(1, 0, 1.0008)
#cluster_i = RobertsClusterFunction(0, 1, 1.1)
#cluster_j = RobertsClusterFunction(1, 0, 1.1)

cflist = [None,None,None,None,]*2 + [cluster_k,]*4;
bclist = [AdjacentBC(),AdjacentBC(),SlipWallBC(),AdjacentBC(),ExtrapolateOutBC(),SupInBC(inflow_condition=injectCond)];
blk10 = SuperBlock3D(nni=ni2, nnj=nj1, nnk=nk1, nbi=1, nbj=1, nbk=8, parametric_volume=WireFrameVolume(Line001,Line005,Line008,Line004,Line101,Line105,Line108,Line104,Vert001,Vert002,Vert012,Vert011),
                cf_list=cflist, bc_list=bclist, fill_condition=initialCond)

#BLOCK 11
cluster_k  = RobertsClusterFunction(1, 0, 1.0008)
#cluster_i = RobertsClusterFunction(0, 1, 1.1)
cluster_j = RobertsClusterFunction(1, 0, 1.1)

cflist = [None,cluster_j,None,cluster_j,]*2 + [cluster_k,]*4;
bclist = [ExtrapolateOutBC(),AdjacentBC(),AdjacentBC(),AdjacentBC(),ExtrapolateOutBC(),AdiabaticBC()];
blk11 = SuperBlock3D(nni=ni2, nnj=nj2, nnk=nk1, nbi=1, nbj=2, nbk=8, parametric_volume=WireFrameVolume(Line008,Line012,Line015,Line011,Line108,Line112,Line115,Line111,Vert011,Vert012,Vert022,Vert021),
                cf_list=cflist, bc_list=bclist, fill_condition=initialCond)

# blocks downstream of injector

#BLOCK 20
cluster_k  = RobertsClusterFunction(1, 0, 1.0008)
cluster_i = RobertsClusterFunction(1, 0, 1.1)
#cluster_j = RobertsClusterFunction(1, 0, 1.1)

cflist = [cluster_i,None,cluster_i,None,]*2 + [cluster_k,]*4;
bclist = [AdjacentBC(),ExtrapolateOutBC(),SlipWallBC(),AdjacentBC(),ExtrapolateOutBC(),AdiabaticBC()];
blk20 = SuperBlock3D(nni=ni3, nnj=nj1, nnk=nk1, nbi=2, nbj=1, nbk=8, parametric_volume=WireFrameVolume(Line002,Line006,Line009,Line005,Line102,Line106,Line109,Line105,Vert002,Vert003,Vert013,Vert012),
                cf_list=cflist, bc_list=bclist, fill_condition=initialCond)

#BLOCK 21
cluster_k  = RobertsClusterFunction(1, 0, 1.0008)
cluster_i = RobertsClusterFunction(1, 0, 1.1)
cluster_j = RobertsClusterFunction(1, 0, 1.1)

cflist = [cluster_i,cluster_j,cluster_i,cluster_j,]*2 + [cluster_k,]*4;
bclist = [ExtrapolateOutBC(),ExtrapolateOutBC(),AdjacentBC(),AdjacentBC(),ExtrapolateOutBC(),AdiabaticBC()];
blk21 = SuperBlock3D(nni=ni3, nnj=nj2, nnk=nk1, nbi=2, nbj=2, nbk=8, parametric_volume=WireFrameVolume(Line009,Line013,Line016,Line012,Line109,Line113,Line116,Line112,Vert012,Vert013,Vert023,Vert022),
                cf_list=cflist, bc_list=bclist, fill_condition=initialCond)

identify_block_connections()

