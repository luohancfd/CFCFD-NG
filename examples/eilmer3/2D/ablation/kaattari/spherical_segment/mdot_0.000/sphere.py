## \file ablating_sphere.py
## \brief Test simulation for AblatingBC class
## \author DFP, ported from cyl50 example, 21-Oct-2009
##

global x_from_y
from cfpylib.gasdyn.billig import x_from_y, y_from_x
global tan, pi
from math import cos, sin, tan, sqrt, pi, asin

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
global bx_scale, by_scale
bx_scale = 1.01
by_scale = 1.01
global Rn
Rn = 21.3e-2
Rb = 8.9e-2
global theta
theta = asin( Rb / Rn )
Rs = 1.0e-2
Ls = 1.0e-2
# origin
o = Node(0.0,0.0)
# sphere surface nodes
a = Node(-Rn,0.0)
b = Node(-Rn*cos(theta),Rn*sin(theta))
s0 = Node(b.x + Rs * cos( theta ), b.y - Rs * sin( theta ),label='s0')
c = Node(s0.x, s0.y + Rs,label='c')
d = Node(c.x+Ls,c.y,label='d') 
# inflow boundary nodes
y_top = by_scale * y_from_x(-d.x/bx_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
inflow_nodes = []
np = 12
dy = y_top / ( np - 1 )
for iy in range(np):
    y = dy * iy
    x = - bx_scale * x_from_y(y/by_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
    inflow_nodes.append( Node(x,y) )
# scale the top points to make more round
inflow_nodes[-1].y *= 0.90
inflow_nodes[-2].y *= 0.92
inflow_nodes[-3].y *= 0.94
inflow_nodes[-4].y *= 0.96
inflow_nodes[-5].y *= 0.98
# intersection inflow node
def zero_func( y ):
    # from cfpylib.gasdyn.billig import x_from_y
    return - bx_scale * x_from_y(y/by_scale, M_inf, theta=0.0, axi=1, R_nose=Rn) - y/tan(pi-theta)
# bisection method
by = 0.0
uy = Rn
while abs(uy-by) > 1e-9:
    mp = ( by + uy ) / 2
    if zero_func(by) * zero_func(mp) > 0: by = mp
    else: uy = mp
y_int = 0.5*(uy+by)
x_int = y_int/tan(pi-theta)
# split inflow_nodes list in two parts
for i in range(len(inflow_nodes)):
    if inflow_nodes[i].y > y_int: break
ef_nodes = inflow_nodes[:i]; ef_nodes.append(Node(x_int,y_int))
fg_nodes = inflow_nodes[i:]; fg_nodes.insert(0,Node(x_int,y_int))
# curves
# block0
west0 = Spline(ef_nodes)
south0 = Line(ef_nodes[0],a)
north0 = Line(ef_nodes[-1],b)
east0 = Arc(a,b,o)
# block1
west1 = Spline(fg_nodes)
south1 = Line(ef_nodes[-1],b)
north1 = Line(fg_nodes[-1],d)
east1 = Polyline([Arc(b,c,s0),Line(c,d)])


# Define the blocks, boundary conditions and set the discretisation.
nnx = 40; nny0 = 30; nny1 = 10
nbx = 3; nby0 = 3; nby1 = 1 

betaNS = 1.0              # clustering
betaSN = 1.0
betaWE = 1.02
betaEW = 1.02

cluster_functions = [RobertsClusterFunction(0, 1, betaEW),
                     RobertsClusterFunction(1, 0, betaNS),
                     RobertsClusterFunction(0, 1, betaWE),
                     RobertsClusterFunction(1, 0, betaSN)]

if with_ablation:
    mdot = -mdot_inf*mdot_ratio
    print "Ablating wall boundary: mdot_inf = %e kg/m2-s, mdot_wall = %e kg/m2-s" % ( mdot_inf, mdot )
    blk0_eastBC = AblatingBC(T_wall,[mdot])
else:
    blk0_eastBC = FixedTBC(T_wall)

blk_0 = SuperBlock2D(psurf=make_patch(north0, east0, south0, west0),
		     fill_condition=initial,
		     nni=nnx, nnj=nny0,
		     nbi=nbx, nbj=nby0,
		     cf_list=cluster_functions,
		     bc_list=[AdjacentBC(), blk0_eastBC, SlipWallBC(), SupInBC( inflow )],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-0", hcell_list=[(nnx,1)])
             
blk_1 = SuperBlock2D(psurf=make_patch(north1, east1, south1, west1),
		     fill_condition=initial,
		     nni=nnx, nnj=nny1,
		     nbi=nbx, nbj=nby1,
		     cf_list=cluster_functions,
		     bc_list=[ExtrapolateOutBC(), FixedTBC(T_wall), SlipWallBC(), SupInBC( inflow )],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-0", hcell_list=[(nnx,1)])
             
identify_block_connections()

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = Rb * 2 * 10 / u_inf    # 10 body diameters
gdata.max_step = 230000
gdata.dt = 1.0e-10
gdata.stringent_cfl = 0
gdata.dt_plot = Rb * 2 / u_inf    # 10 solutions
gdata.cfl = 0.5
gdata.cfl_count = 1

sketch.scales(0.3/Rn, 0.3/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-1.5*Rn, 0.0, 0.25*Rn, -0.1*Rn)
sketch.yaxis(0.0, 1.5*Rn, 0.25*Rn, 0.0)
