## \file tc2_3D.py
## \brief 3D Segment for Test Case 2, using tc2runup as input profile
##		-Developed as part of Samuel Stennett's Undergraduate Thesis, 2014
##		-Utilising modified code written by W. Chan, P. Jacobs, R. Gollan.
##
## \author Samuel Stennett 4 August 2014
##

from cfpylib.gasdyn import sutherland

gdata.title = "Kim's Mach 4 deflected flow (k-omega)"
print gdata.title
gdata.turbulence_model = "k_omega"
gdata.dimensions = 3
gdata.viscous_flag = 1
gdata.flux_calc = AUSMDV
        
gdata.max_time = 1.6e-3  # About 3 flow lengths (1 flow length ~ 0.56 ms)  #####
gdata.max_step = 30000000
gdata.dt_plot = 0.5e-6

gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-10  # Only an initial guess - the simulation will take over

select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Mabey's data set 74021802
p_inf = 10.3e3  # Pa
u_inf = 707.0  # m/s
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

initial = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,],
                       tke=tke_inf, omega=omega_inf )

def toRadians(degrees):
    import math
    return degrees * math.pi / 180.0

def simpleBoxCorners(xPos=0.0, yPos=0.0, zPos=0.0, xSize=1.0, ySize=1.0, zSize=1.0):
    """\brief Creates a corner coordinate list for a simple box."""
    p0 = Node(xPos,       yPos,       zPos)
    p1 = Node(xPos+xSize, yPos,       zPos)
    p2 = Node(xPos+xSize, yPos+ySize, zPos)
    p3 = Node(xPos,       yPos+ySize, zPos)
    p4 = Node(xPos,       yPos,       zPos+zSize)
    p5 = Node(xPos+xSize, yPos,       zPos+zSize)
    p6 = Node(xPos+xSize, yPos+ySize, zPos+zSize)
    p7 = Node(xPos,       yPos+ySize, zPos+zSize)
    return [p0, p1, p2, p3, p4, p5, p6, p7]

def makeSimpleBox(p):
    return SimpleBoxVolume(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7])


#block definiton

##INPUT ANGLE##
alpha1 = 16.0 #degrees, as well as 20 degrees

block0Corners = simpleBoxCorners(xSize=0.01,ySize=0.075+0.09*math.sin(toRadians(alpha1)),zSize=0.02)
cluster_i = RobertsClusterFunction(0,1,1.05)
cluster_j = RobertsClusterFunction(0,1,1.0015)
cluster_k = RobertsClusterFunction(1,0,1.0013)
cflist = [cluster_i,cluster_j,]*4+[cluster_k,]*4;
bclist = [SlipWallBC(),AdjacentBC(),ExtrapolateOutBC(),StaticProfBC(filename="REPLACE.dat",n_profile=2),ExtrapolateOutBC(),AdiabaticBC()]

blk0=SuperBlock3D(nni=15,nnj=100,nnk=100,nbi=1,nbj=4,nbk=4,parametric_volume=makeSimpleBox(block0Corners),cf_list=cflist,fill_condition=initial,bc_list=bclist)

block1Corners = simpleBoxCorners(xSize=0.01+0.09*math.cos(toRadians(alpha1)),ySize=0.075+0.09*math.sin(toRadians(alpha1)),zSize=0.02)
cluster_i = RobertsClusterFunction(1,0,1.03)
cluster_j1 = RobertsClusterFunction(0,1,1.0015) #small edge
cluster_j2 = RobertsClusterFunction(0,1,1.0015) #large edge
cluster_k = RobertsClusterFunction(1,0,1.0013)
cflist = [cluster_i,cluster_j1,cluster_i,cluster_j2,]*2+[cluster_k,]*4;
bclist = [AdiabaticBC(),ExtrapolateOutBC(),ExtrapolateOutBC(),AdjacentBC(),ExtrapolateOutBC(),AdiabaticBC()]

#Shift corners, to create deflection 
block1Corners[2].y=0.075;
block1Corners[6].y=0.075;
block1Corners[0].x=0.01;
block1Corners[3].x=0.01;
block1Corners[4].x=0.01;
block1Corners[7].x=0.01;
blk1=SuperBlock3D(nni=85,nnj=100,nnk=100,nbi=3,nbj=4,nbk=4,parametric_volume=makeSimpleBox(block1Corners),cf_list=cflist,fill_condition=initial,bc_list=bclist)

#Connecting Blocks
identify_block_connections()

