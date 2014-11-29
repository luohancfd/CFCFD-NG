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
##	   Samuel Stennett April 2014 - Converted to 3D

#Initialise sutherland for gas properties
from cfpylib.gasdyn import sutherland

#Initialise simulation attributes
gdata.title = "Mabey's Mach 4.5 flow over a flat plate (k-omega)"
print gdata.title
gdata.turbulence_model = "k_omega"
gdata.dimensions = 3
gdata.viscous_flag = 1
gdata.flux_calc = AUSMDV    
gdata.max_time = 1.6e-3  # About 3 flow lengths (1 flow length ~ 0.56 ms)  
gdata.dt_history = 1.0e-3
gdata.max_step = 3000000
gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-9  # Only an initial guess - the simulation will take over

#Assume ideal gas, use properties for air
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

block1Corners = simpleBoxCorners(xSize=0.4,ySize=0.01,zSize=0.16)
cluster_i = RobertsClusterFunction(1,0,1.05)
cluster_k1 = RobertsClusterFunction(0,1,1.0014)
cluster_k2 = RobertsClusterFunction(0,1,1.0074)
cflist = [cluster_i,None,]*4+[cluster_k2,cluster_k1,cluster_k1,cluster_k2,];
bclist = [SlipWallBC(),ExtrapolateOutBC(x_order=0,sponge_flag=0),SlipWallBC(),SupInBC(inflow_condition=inflow),AdiabaticBC(),SupInBC(inflow_condition=inflow)]

block1Corners[0].z=0.12;
block1Corners[3].z=0.12;

blk0=SuperBlock3D(label="turb3D",nni=128,nnj=10,nnk=96,nbi=4,nbj=2,nbk=2,parametric_volume=makeSimpleBox(block1Corners),cf_list=cflist,fill_condition=initial,bc_list=bclist)

