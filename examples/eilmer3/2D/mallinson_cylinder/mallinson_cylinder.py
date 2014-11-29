## \file mallinson_cylinder.py
## \brief Mach 9 turbulent flow over a hollow cylinder
## \author Wilson Chan
## \version Wilson Chan, 24-Oct-2008                           
##
#  Reference: - Mallinson SG et al. (2000), Shock Waves, v.10, pp.313-322
#             - Boyce RR & Hillier R (2000), AIAA 2000-2226

# Basically, we are increasing axial resolution of the flow domain
# in attempt to remove the instabilities in the surface heat flux
# results. For these changes, just look for the # INCREASE AXIAL RESOLUTION
# comment. We will attempt this in 2 ways ..
#   1. Shorten domain using the same amount of grids (1st try), or
#   2. Increase number of cells in the axial direction. 

# Now on 22 June 2009, we will try to see what happens when we allow the
# k-omega model transition naturally.. or will it even transition??

# Now on 16 July 2009, we will try to extend the flow domain back to the
# original length (while trying to keep the axial grid resolution similar
# to the resolution in the shorter flow domain).

# Now on 18 July 2009, after we have reverted to r502 of the cfcfd2 package
# (where major changes have been made to how viscous terms were computed),
# we have to now re-adjust the x_tr location downstream. A good result that
# came out from this new revision is that the noise we had before in the
# heat transfer results is all gone now!! Must inform Peter about this..

# Modified by Samuel Stennett, Nov 2014
# Removed redundant code, revised outdated sections, combined separate cases
# (different cell counts, clustering, blocking) into this version.

job_title = "Mach 9 flow over a hollow cylinder (k-omega)"
print job_title

# Set global data
gdata.title = job_title
gdata.dimensions = 2        # 2D calculations
gdata.axisymmetric_flag = 1 # Axisymmetric calculations
gdata.viscous_flag = 1      # Viscous calculation
gdata.turbulence_model = "k_omega"
gdata.flux_calc = ADAPTIVE

gdata.max_time = 5.0e-3  
gdata.dt_plot =  0.1e-3
gdata.max_step = 30000000

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.stringent_cfl = 0  # 1 is more robust
gdata.dt = 1.0e-14  # only an initial guess, the simulation will take over

# Select gas model - accept defaults for N2 giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['N2'])

# Define inflow conditions
p_inf = 3.3e3    # Pa
u_inf = 1498.0   # m/s
T_inf = 69.7     # degrees K
rho_inf = 0.159  # kg/m**3

# Set inflow turbulence quantities to almost-laminar values
tke_inf = 1.0e-12   
omega_inf = 1.0   
print "Inflow turbulence: tke=", tke_inf, "omega=", omega_inf

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,],
                       tke=tke_inf, omega=omega_inf )

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
                            tke=tke_inf/100.0, omega=omega_inf/10.0 )

# Define geometry (all dimensions in metres)
x0 = 0.0
x1 = 0.85
y0 = 0.0
y1 = 0.0375
y2 = 0.057
y3 = 0.2 

# Build lines with nodes
north = Line(Node(x0,y2), Node(x1,y3))
east = Line(Node(x1,y1), Node(x1,y3))
south = Line(Node(x0,y1), Node(x1,y1))
west = Line(Node(x0,y1), Node(x0,y2))

# Point to define TurbulenceZone
# However, since we have set the inflow turbulence quantities to 
# almost-laminar values and allowed the k-omega model to transition
# the boundary layer (transitions at about 0.031 m) and we want the
# boundary layer to transition at 0.085 m, we'll set the TurbulenceZone
# to start at 0.071 - 0.031 = 0.040 metres.
x_tr = 0.047    # Previously 0.04 m.. Added 0.007 m for r502 cfcfd2 package.

# Define turbulent zone
# First vector defines lower left corner of turbulent zone
# Second vector defines upper right corner of turbulent zone
TurbulenceZone(Vector(x_tr,y1,0.0), Vector(x1,y3,0.0))

casenum = 0 # SELECT A CASE, FROM 0 to 8
# CASES DEFINED AS FOLLOWS:	
#	0:	100x33cells-yplus3.0-ar215
#	1:	200x65cells-yplus1.5-ar227
#	2:	200x65cells-yplus0.4-ar922
#	3:	400x130cells-yplus0.8-ar235
#	4:	400x130cells-yplus1.7-ar114
#	5:	400x130cells-yplus0.2-ar961
#	6:	400x130cells-yplus0.2-ar864
#	7:	600x130cells-yplus0.2-ar577
#	8:	520x169cells-yplus0.6-ar237

# 100x33cells-yplus3.0-ar215
if casenum == 0:
	nC_I = 100
	nC_J = 33
	nB_I = 4
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.001), 
                  RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.01)]
# 200x65cells-yplus1.5-ar227
elif casenum == 1:
	nC_I = 200
	nC_J = 65
	nB_I = 8
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.001), 
                  RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.01)]
# 200x65cells-yplus0.4-ar922
elif casenum == 2:
	nC_I = 200
	nC_J = 65
	nB_I = 8
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.0002), 
                  RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.002)]
# 400x130cells-yplus0.8-ar235
elif casenum == 3:
	nC_I = 400
	nC_J = 130
	nB_I = 16
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.001), 
                  RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.01)]
# 400x130cells-yplus1.7-ar114
elif casenum == 4:
	nC_I = 400
	nC_J = 130
	nB_I = 16
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.002), 
                  RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.02)]
# 400x130cells-yplus0.2-ar961
elif casenum == 5:
	nC_I = 400
	nC_J = 130
	nB_I = 16
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.0002), 
                  RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.002)]
# 400x130cells-yplus0.2-ar864
elif casenum == 6:
	nC_I = 400
	nC_J = 130
	nB_I = 16
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.05), 
                  RobertsClusterFunction(1,0, 1.0002), 
                  RobertsClusterFunction(1,0, 1.05), 
                  RobertsClusterFunction(1,0, 1.002)]
# 600x130cells-yplus0.2-ar577
elif casenum == 7:
	nC_I = 600
	nC_J = 130
	nB_I = 8
	nB_J = 2
	cflist = [RobertsClusterFunction(1,0, 1.05), 
                  RobertsClusterFunction(1,0, 1.0002), 
                  RobertsClusterFunction(1,0, 1.05), 
                  RobertsClusterFunction(1,0, 1.002)]
# 520x169cells-yplus0.6-ar237
elif casenum == 8:
	nC_I = 520
	nC_J = 169
	nB_I = 16
	nB_J = 1
	cflist = [RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.001), 
                  RobertsClusterFunction(1,0, 1.03), 
                  RobertsClusterFunction(1,0, 1.01)]
else:
	print("Not a valid case, prepare to crash...")

# Define the blocks, boundary conditions and set the discretisation
blk = SuperBlock2D(make_patch(north, east, south, west),
                   nni=nC_I, nnj=nC_J, nbi=nB_I, nbj=nB_J,
                   fill_condition=initial,
                   cf_list=cflist,
                   bc_list=[UserDefinedBC("udf-supersonic-in.lua"), ExtrapolateOutBC(),
                            FixedTBC(295.0), UserDefinedBC("udf-supersonic-in.lua")])
