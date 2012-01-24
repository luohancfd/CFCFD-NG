## \file ablating_sphere.py
## \brief Test simulation for AblatingBC class
## \author DFP, ported from cyl50 example, 21-Oct-2009
##

global x_from_y
from cfpylib.gasdyn.billig import x_from_y, y_from_x
global tan, pi
from math import cos, sin, tan, sqrt, pi

job_title = "Mach 7.32 flow over an 8.9cm radius sphere."
print job_title

# We can set individual attributes of the global data object.
gdata.title = job_title
gdata.axisymmetric_flag = 1

# to start from a saved partA::9999 solution
transfer_solution = False
# to turn on ablation
with_ablation = False 
mdot_ratio = 0.004

select_gas_model(model='ideal gas', species=['air'])
gm = get_gas_model_ptr()
Q = Gas_data(gm)

# Define flow conditions
T_wall = 300.0
p_inf = 552.0
T_inf = 300.0
global M_inf
M_inf = 7.32
massf_inf = [ 1.0 ]
# do some calculations to get velocity and mass flux
for itm in range(gm.get_number_of_modes()):
    Q.T[itm] = T_inf
Q.p = p_inf
for isp,massf in enumerate(massf_inf):
    Q.massf[isp] = massf
gm.eval_thermo_state_pT(Q)
u_inf = M_inf * Q.a
rho_inf = Q.rho
mdot_inf = u_inf * rho_inf
print "Freestream velocity corresponding to Mach %0.2f is %0.1f" % ( M_inf, u_inf )
inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=massf_inf)
if transfer_solution:
    initial = ExistingSolution(rootName="sphere", solutionWorkDir="../partA/", nblock=12, tindx=9999)
else:
    initial = FlowCondition(p=p_inf/10.0, u=0.0, v=0.0, T=T_inf, massf=massf_inf)

# setup a quarter sphere with the specified radius and an inflow boundary
# described by the Billig correlation
global Rn
Rn = 8.9e-2
global bx_scale, by_scale
bx_scale = 1.04
by_scale = 1.04
# origin
o = Node(0.0,0.0)
# sphere surface nodes
a = Node(-Rn,0.0)
b = Node(0.0,Rn)
# inflow boundary nodes
inflow_nodes = []
np = 16
y_top = by_scale * y_from_x(b.x/bx_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
dy = y_top / ( np - 1 )
for iy in range(np):
    y = dy * iy
    x = - bx_scale * x_from_y(y/by_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
    inflow_nodes.append( Node(x,y) )
    print inflow_nodes[-1].str()
# curves
west = Spline(inflow_nodes)
south = Line(inflow_nodes[0],a)
north = Line(inflow_nodes[-1],b)
east = Arc(a,b,o)

# Define the blocks, boundary conditions and set the discretisation.
nnx = 60; nny = 60
nbx = 4; nby = 4;

betaNS = 1.0              # clustering
betaSN = 1.0
betaWE = 1.01
betaEW = 1.01

cluster_functions = [RobertsClusterFunction(0, 1, betaEW),
                     RobertsClusterFunction(1, 0, betaNS),
                     RobertsClusterFunction(0, 1, betaWE),
                     RobertsClusterFunction(1, 0, betaSN)]

if with_ablation:
    mdot = -mdot_inf*mdot_ratio
    print "Ablating wall boundary: mdot_inf = %e kg/m2-s, mdot_wall = %e kg/m2-s" % ( mdot_inf, mdot )
    eastBC = AblatingBC(T_wall,[mdot])
else:
    eastBC = FixedTBC(T_wall)

blk_0 = SuperBlock2D(psurf=make_patch(north, east, south, west),
		     fill_condition=initial,
		     nni=nnx, nnj=nny,
		     nbi=nbx, nbj=nby,
		     cf_list=cluster_functions,
		     bc_list=[ExtrapolateOutBC(), eastBC, SlipWallBC(), SupInBC( inflow )],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-0", hcell_list=[(nnx,1)])
             
identify_block_connections()

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = Rn * 20 / u_inf    # 20 body lengths
gdata.max_step = 2300000
gdata.dt = 1.0e-10
gdata.stringent_cfl = 0
gdata.dt_plot = Rn * 1 / u_inf    # 20 solutions
gdata.cfl = 0.5

sketch.scales(0.3/Rn, 0.3/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-1.5*Rn, 0.0, 0.25*Rn, -0.1*Rn)
sketch.yaxis(0.0, 1.5*Rn, 0.25*Rn, 0.0)
