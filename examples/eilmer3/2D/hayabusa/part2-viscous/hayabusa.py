## \file hayabusa.py
## \brief Simulating the JAXA Hayabusa sample return capsule
## \author DFP, 30-Oct-2012
## 
## Part2: Viscous solution on a finer grid

from cfpylib.gasdyn.billig import x_from_y, y_from_x
from cfpylib.nm.zero_solvers import bisection
from math import *
from cfpylib.flow.shock_layer_surface import ShockLayerSurface

# for shock fitting
import matplotlib.pyplot as plt
from scipy.optimize import *
from numpy import *
from cfpylib.util.YvX import *

def BRH(a,b,beta,dx0,dx1,gamma):
    RCF = RobertsClusterFunction(a,b,beta0)
    HCF = HypertanClusterFunction(dx0,dx1)
    return DiscontinuousUnivariateFunction(gamma,RCF,HCF)

job_title = "JAXA Hayabusa sample return capsule."
print job_title

gdata.title = job_title
gdata.axisymmetric_flag = 1

fit2shock = True

#
# 1. Setup the gas model
#
species = select_gas_model(model='two temperature gas', species=['N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus'])
set_reaction_update("Park93-s03-AIC-EIIC.lua")
set_energy_exchange_update("air-TV-TE.lua")
gm = get_gas_model_ptr()
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()

#
# 2. Define flow conditions
#
rho_inf = 1.645e-4
T_inf = 233.25
u_inf = 11.6e3
massf_inf = [ 0.0 ] * gm.get_number_of_species()
massf_inf[species.index('N2')] = 0.767
massf_inf[species.index('O2')] = 0.233

# do some calculations to get pressure, Mach number and total mass-flux
Q = Gas_data(gm)
for itm in range(ntm):
    Q.T[itm] = T_inf
Q.rho = rho_inf
for isp,massf in enumerate(massf_inf):
    Q.massf[isp] = massf
gm.eval_thermo_state_rhoT(Q)
M_inf = u_inf / Q.a
p_inf = Q.p
print "p_inf = %0.1f, M_inf = %0.1f" % ( p_inf, M_inf )
inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=[T_inf]*ntm, massf=massf_inf)

# use the inviscid solution as the initial condition
initial = ExistingSolution(rootName="hayabusa", solutionWorkDir="../part1-inviscid/", nblock=4, tindx=9999)

# find the shock points from the previous solution
shock_points = []
shock_y_coords = []
shock_x_coords = []
for flow in initial.flow:
    for i in range(flow.ni):
        for j in range(flow.nj):
            for k in range(flow.nk):
                if flow.data["S"][i][j][k]==1:
                    print flow.data.keys()
                    shock_points.append( Vector3( flow.data["pos.x"][i][j][k], flow.data["pos.y"][i][j][k], flow.data["pos.z"][i][j][k] ) ) 
                    shock_x_coords.append( flow.data["pos.x"][i][j][k] )
                    shock_y_coords.append( flow.data["pos.y"][i][j][k] )
print "Found %d shock points" % len(shock_points)
shock_x_coords, shock_y_coords = zip(*sorted(zip(shock_x_coords, shock_y_coords)))

#
# 3. Define the geometry
#
Rn = 200.0e-3
D = 400.0e-3
theta = 45.0 * pi / 180.0
L = 200.0e-3
global bx_scale, by_scale
bx_scale = 0.95
by_scale = 0.95

# origin
o = Node(0.0,0.0)

# surface nodes
global c,d
a = Node(-Rn,0.0, label='a')
b = Node(-Rn*cos(theta),Rn*sin(theta), label='b')
c = Node( b.x + ( D/2.0-b.y)/tan(theta), D/2.0, label='c')
d = Node( 0.0, c.y - abs(c.x), label='d')

body = Polyline( [Arc(a,b,o),Line(b,c)] )

# inflow boundary nodes
x_limit = c.x
inflow_nodes = []
np = 32
if fit2shock:
    # make a best fit to tthe shock location point cloud
    def shock_x_from_y( p, y ):
        print p
        _bx_scale, _by_scale, _M_inf, _Rn = p
        x = []
        for _y in y:
            _x = - _bx_scale * x_from_y(_y/_by_scale, _M_inf, theta=0.0, axi=1, R_nose=_Rn)
            x.append( _x )
        return array(x)

    def shock_y_from_x( p, x ):
        print p
        _bx_scale, _by_scale, _M_inf, _Rn = p
        y = []
        for _x in x:
            _y = _by_scale * y_from_x(-_x*_bx_scale, _M_inf, theta=0.0, axi=1, R_nose=_Rn)
            y.append( _y )
        return array(y)

    def residuals( p, x, y ):
        x_dash = shock_x_from_y(p,y)
        return sqrt(sum((array(x) - x_dash)**2)/len(x))
    
    p0 = [ bx_scale, by_scale, M_inf, Rn ]
    plsq = fmin_slsqp(residuals, p0, args=(shock_x_coords, shock_y_coords), bounds=[(1.0e-1,1e1),(1.0e-1,1e1),(1.0e0,1e2),(1.0e-4,1e0)], fprime=None)
    # p = p0
    p = plsq
    fit_x = []
    fit_y = copy(shock_y_coords)
    fit_y.sort()
    insert(fit_y,0,0.0)
    fit_x = shock_x_from_y(p,fit_y)
    
    plt.plot(fit_x,fit_y,"b-",label="fit")
    plt.plot(shock_x_coords,shock_y_coords,"g.",label="points")
    plt.grid()
    plt.legend()
    plt.show()  

    y_top = p[1] * y_from_x(-x_limit/p[0], p[2], theta=0.0, axi=1, R_nose=p[3])
    dy = y_top / ( np - 1 )
    tmp_shock_nodes = []
    for iy in range(np):
        y = dy * iy
        x = - p[0] * x_from_y(y/p[1], p[2], theta=0.0, axi=1, R_nose=p[3])
        tmp_shock_nodes.append( Node(x,y) )
    shock_spline = Spline( tmp_shock_nodes )

    inflow_nodes = []
    for _is in range(np):
        s = _is / float( np - 1 )
        # point on the body
        bp = body.eval(s)
	# eval normal angle on the body
	dpdt = body.dpdt(s)
	if dpdt.x == 0.0:
            theta = 0.0
        else:
  	    theta = - atan( dpdt.y / dpdt.x ) + pi / 2.	
	#  find the shock point
	if s==0.0: t = 0
        elif s==1.0: t = 1
	elif theta != 0.0:
            def f( t ):
	        wp = shock_spline.eval(t)
	        L = ( wp.y - bp.y ) / sin( theta ) 
	        lp = Vector3( bp.x - L * cos(theta), bp.y + L * sin(theta) ) 
	        return wp.x - lp.x
	    t = bisection( f, 0.0, 1.0 )
	else:
            def f( t ):
	        wp = shock_spline.eval(t)
	        L = ( bp.x - wp.x ) / cos( theta ) 
	        lp = Vector3( bp.x - L * cos(theta), bp.y + L * sin(theta) ) 
	        return wp.y - lp.y
	    t = bisection( f, 0.0, 1.0 )
	sp = shock_spline.eval(t)
        # calculate the inflow_node location by extending by some factor
        f_extend = 1.2
        inflow_nodes.append( Node( bp.x + ( sp.x - bp.x ) * f_extend, bp.y + ( sp.y - bp.y ) * f_extend ) )
else:
    y_top = by_scale * y_from_x(-x_limit/bx_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
    dy = y_top / ( np - 1 )
    for iy in range(np):
        y = dy * iy
        x = - bx_scale * x_from_y(y/by_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
        inflow_nodes.append( Node(x,y) )

# find intersection of surface normal with inflow boundary at the top most point
global inflow_spline
inflow_spline = Spline(inflow_nodes)
def zero_func(y):
    dc_line = Line( d, Node( d.x+(c.x-d.x)*2.0, d.y+(c.y-d.y)*2.0 ) )
    t = (y-d.y)/((c.y-d.y)*2.0)
    return inflow_spline.eval_from_y(y).x - dc_line.eval(t).x
y_int = bisection( zero_func, by=0.0, uy=D )
x_int = inflow_spline.eval_from_y(y_int).x

# split the inflow spline
west_nodes = []
for node in inflow_nodes:
    if node.y < y_int: west_nodes.append(node)
west_nodes.append( Node( x_int, y_int ) )

# curves - block0
west0 = Spline(west_nodes)
east0 = Polyline( [Arc(a,b,o),Line(b,c)] )

#
# 4. Define the blocks, boundary conditions and set the discretisation
#
nnx = 60; nny=60
nbx = 2; nby = 2

# clustering at shock and boundary layer [ outflow, surface, axis, inflow ]
beta0 = 1.1; dx0 = 5.0e-1; dx1 = 2.0e-2; gamma = 0.2
beta1 = 1.2
cf_list = [BRH(0,1,beta0,dx0,dx1,gamma),
           RobertsClusterFunction(0,1,beta1),
           BRH(0,1,beta0,dx0,dx1,gamma),
           RobertsClusterFunction(0,1,beta1)]

# boundary conditions [ outflow, surface, axis, inflow ]
bc_list=[ ExtrapolateOutBC(),
          FixedTBC(3000.0),
          SlipWallBC(),
          SupInBC(inflow)]

# catalycity [ outflow, surface, axis, inflow ]
wc_bc_list = [NonCatalyticWBC(),
              SuperCatalyticWBC(massf_inf),
              NonCatalyticWBC(),
              NonCatalyticWBC()]

def RCF(a,b,beta):
    return RobertsClusterFunction(a, b, beta)

blk_0 = SuperBlock2D(psurf=ShockLayerSurface(east0, west0),
		     fill_condition=initial,
		     nni=nnx, nnj=nny,
		     nbi=nbx, nbj=nby,
		     cf_list=cf_list,
		     bc_list=bc_list,
		     wc_bc_list=wc_bc_list,
		     label="BLOCK-0", hcell_list=[(nnx,1)])
 
#
# 5. Simulation control parameters 
#
gdata.viscous_flag = 1
gdata.viscous_delay = 0.1 * Rn / u_inf
gdata.viscous_factor_increment = 1.0e-4
gdata.diffusion_flag = 1
gdata.diffusion_model = "ConstantLewisNumber"
gdata.flux_calc = ADAPTIVE
gdata.max_time = Rn * 5 / u_inf    # 5 body lengths
gdata.max_step = 230000
gdata.dt = 1.0e-9
gdata.reaction_time_start = 0
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1 / u_inf    # 5 solutions
gdata.cfl = 0.5
gdata.cfl_count = 1
gdata.print_count = 1

#
# 6. svg sketch parameters
#
sketch.scales(0.3/Rn, 0.3/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-1.5*Rn, 0.0, 0.25*Rn, -0.1*Rn)
sketch.yaxis(0.0, 1.5*Rn, 0.25*Rn, 0.0)
