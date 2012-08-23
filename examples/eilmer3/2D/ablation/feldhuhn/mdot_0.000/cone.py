## \file cone.py
## \brief Job-specification file for 5 degree sphere-cone in
## \brief Mach 5 flow (White Oak Lab. Hypersonic Tunnel)
## \author DFP et al, 02-November-2009
##
## Ref: Marvin, J.G. and Akin, C.M.
##      Combined Effects of Mass Addition and Nose Bluntness 
##      on Boundary-Layer Transition
##      AIAA Journal Vol. 8 No. 5 pp 857-863 May 1970 

import sys
global cos, pi, sin
from math import cos, pi, sin
from gaspy import *
from cfpylib.gasdyn.billig import x_from_y, y_from_x

job_title = "5 degree sphere-cone in Mach 5 Flow."
print job_title

# We can set individual attributes of the global data object.
gdata.title = job_title
gdata.axisymmetric_flag = 1

gas_model = "IG"
transfer_solution = False
with_ablation = False
global with_blk3; with_blk3 = False
injected_gas = 'N2' 
mdot_ratio = 0.0
cell_scale = 1.0
viscous_delay_BLs = 1
max_BLs = 10
plot_BLs = 1

# use ideal gas-model
if gas_model == "IG":
    select_gas_model(model='ideal gas', species=['air',injected_gas])
    massf_inf = [ 1.0, 0.0 ]
else:
   print "Only ideal gas model is currently available."
   sys.exit()
gm =  get_gas_model_ptr()

# Define flow condition
M_inf = 5
p_total_atm = 10.0     # total pressure in atmospheres
T_total_R = 945.0      # total temperature in degrees Rankine
gamma_air = 1.4
R_air = 287.0
T_inf = ( T_total_R * 0.555555556 ) / ( 1.0 + ( gamma_air - 1.0 ) / 2.0 * M_inf**2 )
u_inf = M_inf * sqrt( gamma_air * R_air * T_inf )
p_inf = ( p_total_atm * 101.325e3 ) / ( 1.0 + ( gamma_air - 1.0 ) / 2.0 * M_inf**2 )**(gamma_air/(gamma_air-1.0))
Twall_R = 520.0         # wall temperature in degrees Rankine
Twall = Twall_R * 0.555555556
Q = Gas_data(gm)
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()
for itm in range(ntm): Q.T[itm] = T_inf
Q.p = p_inf 
for isp in range(nsp): Q.massf[isp] = massf_inf[isp]
gm.eval_thermo_state_pT(Q)
gm.eval_transport_coefficients(Q)
T_total = ( gm.mixture_enthalpy(Q) + 0.5 * u_inf**2 ) / gm.Cp(Q) 
print "Given T_total = %f K, calculated T_total = %f K" % ( T_total_R * 0.555555556, T_total )
ReL = Q.rho * u_inf / Q.mu
Q.print_values()
print "Reynolds number = %0.3e per foot" % (ReL/3.28084)

inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=massf_inf)

if transfer_solution:
    initial = ExistingSolution(rootName="cone", solutionWorkDir=".", nblock=18, tindx=9999)
else:
    initial = FlowCondition(p=p_inf/10.0, u=0.0, v=0.0, T=T_inf, massf=massf_inf )
    
# Set up the geometry
global theta
Rn = 2.0 * 0.0254           # 2 inch nose radius
Lx = 15.0 * 0.0254          # 5 inch total x length
theta = 5.0 * pi / 180.0    # 5 degree cone
s_Rn_list = [ 0.52, 1.48, 3.24 ]

# firstly the points on the surface
global o, b, d, e
o = Node(0.0, 0.0, label="o")
a = Node(-Rn, 0.0, label="a")
b = Node( -Rn*cos(s_Rn_list[0]), Rn*sin(s_Rn_list[0]), label="b" )
c = Node( -Rn*cos(pi/2-theta), Rn*sin(pi/2-theta), label="c" )
L_flank_c02 = s_Rn_list[1] * Rn - ( pi/2-theta ) * Rn
d = Node( c.x + L_flank_c02 * cos( theta ), c.y + L_flank_c02 * sin( theta ), label="d" )
L_flank_c23 = s_Rn_list[2] * Rn - s_Rn_list[1] * Rn
e = Node( d.x + L_flank_c23 * cos( theta ), d.y + L_flank_c23 * sin( theta ), label="e" )
f = Node( a.x + Lx, ( e.y - c.y ) / ( e.x - c.x ) * ( a.x + Lx - c.x ) + c.y, label="f" )

# now the points on the curved inflow boundary
by_scale = 1.04; bx_scale = 1.04
inflow_points = []
np = 30
y_top = by_scale * y_from_x(-f.x/bx_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
dy = ( y_top ) / ( np - 1 )
for iy in range( np ):
    y = dy * iy
    x = - bx_scale * x_from_y(y/by_scale, M_inf, theta=0.0, axi=1, R_nose=Rn)
    inflow_points.append( Node(x,y) )
global inflow_spline
inflow_spline = Spline( inflow_points )

def bisection_method( f, by, uy ):
	while abs(uy-by) > 1e-6:
	    mp = ( by + uy ) / 2
	    if f(by) * f(mp) > 0:
		by = mp
	    else:
		uy = mp
	return 0.5*(uy+by)
    
def restrict_points_list( plist, p0, p1 ):
   new_plist = [ p0 ]
   for pt in plist:
       if pt.y <= p0.y*1.01 or pt.y >= p1.y*0.99: continue
       else: new_plist.append( pt )
   new_plist.append(p1)
   return new_plist
    
# calc end of chamber 1 inflow coords
g = Node( inflow_points[0].x, inflow_points[0].y, label="g" )
def zero_func( y ):
    oH_line = Line( o, Node( b.x*10.0, b.y*10.0 ) )
    t = y/(b.y*10.0)
    return inflow_spline.eval_from_y(y).x - oH_line.eval(t).x
y_int = bisection_method( zero_func, by=0.0, uy=y_top )
h = Node( inflow_spline.eval_from_y(y_int).x, y_int, label="h" )
gh_points = restrict_points_list( inflow_points, g, h )
# gh_points.pop(0)
# calc end of chamber 2 inflow coords
def zero_func( y ):
    dI_line = Line( d, Node( d.x-10.0*cos(pi/2.0-theta), d.y+10.0*sin(pi/2.0-theta) ) )
    t = y/(10.0*sin(pi/2.0-theta))
    return inflow_spline.eval_from_y(y).x - dI_line.eval(t).x
y_int = bisection_method( zero_func, by=0.0, uy=y_top*2 )
i = Node( inflow_spline.eval_from_y(y_int).x, y_int, label="i" )
hi_points = restrict_points_list( inflow_points, h, i )
# calc end of chamber 3 inflow coords
def zero_func( y ):
    if with_blk3: angle = theta
    else: angle = 0.0
    eJ_line = Line( e, Node( e.x-10.0*cos(pi/2.0-angle), e.y+10.0*sin(pi/2.0-angle) ) )
    t = y/(10.0*sin(pi/2.0-angle))
    return inflow_spline.eval_from_y(y).x - eJ_line.eval(t).x
y_int = bisection_method( zero_func, by=0.0, uy=y_top*2 )
j = Node( inflow_spline.eval_from_y(y_int).x, y_int, label="j" )
ij_points = restrict_points_list( inflow_points, i, j )
# ij_points.pop(-1)
k = Node( inflow_points[-1].x, inflow_points[-1].y, label="k" )
jk_points = restrict_points_list( inflow_points, j, k )
# gh_points.pop(-1)

# block0 boundaries
north0 = Line( h, b )
east0 = Arc( a, b, o )
south0 = Line( g, a )
west0 = Spline( gh_points )

# block1 boundaries
north1 = Line( i, d )
east1 = Polyline( [ Arc(b,c,o), Line(c,d) ] )
south1 = north0
west1 = Spline( hi_points )

# block2 boundaries
north2 = Line( j, e )
east2 = Line( d, e )
south2 = north1
west2 = Spline( ij_points )

# block3 boundaries
north3 = Line( k, f )
east3 = Line( e, f )
south3 = north2
west3 = Spline( jk_points )

# clustering
betaEW = 1.03; betaWE = 1.03; betaSN = 0.0; betaNS = 0.0
cf = [RobertsClusterFunction(0, 1, betaEW),
      RobertsClusterFunction(1, 0, betaNS),
      RobertsClusterFunction(0, 1, betaWE),
      RobertsClusterFunction(1, 0, betaSN)]
       
# cell discretization
inflow_spline = Spline( inflow_points )
nnx = 120
nny = 400
nny0 = int ( cell_scale * west0.length() / inflow_spline.length() * nny )
nny1 = int ( cell_scale * west1.length() / inflow_spline.length() * nny )
nny2 = int ( cell_scale * west2.length() / inflow_spline.length() * nny )
nny3 = int ( cell_scale * west3.length() / inflow_spline.length() * nny )
print "nny0 = %d, nny1 = %d, nny2 = %d, nny3 = %d" % ( nny0, nny1, nny2, nny3 )

# super-blocking
NBX = 3
NBY = 16
NBY0 = int( west0.length() / inflow_spline.length() * NBY )
NBY1 = int( west1.length() / inflow_spline.length() * NBY )
NBY2 = int( west2.length() / inflow_spline.length() * NBY )
NBY3 = int( west3.length() / inflow_spline.length() * NBY )
print "NBY0 = %d, NBY1 = %d, NBY2 = %d, NBY3 = %d" % ( NBY0, NBY1, NBY2, NBY3 )

if with_ablation:
    b0_eastBC = AblatingBC( Twall=Twall, mdot=[0.0, mdot[0] ] )
else:
    b0_eastBC = FixedTBC( Twall )

blk0 = SuperBlock2D(psurf=make_patch(north0, east0, south0, west0),
		     fill_condition=initial,
		     nni=nnx, nnj=nny0,
		     nbi=NBX, nbj=NBY0,
		     cf_list=cf,
		     bc_list=[AdjacentBC(), b0_eastBC, SlipWallBC(), SupInBC( inflow ) ],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-0", hcell_list=[(nnx,1)])
          
if with_ablation:
    b1_eastBC = AblatingBC( Twall=Twall, mdot=[0.0, mdot[1] ] )
else:
    b1_eastBC = FixedTBC( Twall )
             
blk1 = SuperBlock2D(psurf=make_patch(north1, east1, south1, west1),
		     fill_condition=initial,
		     nni=nnx, nnj=nny1,
		     nbi=NBX, nbj=NBY1,
		     cf_list=cf,
		     bc_list=[AdjacentBC(), b1_eastBC, AdjacentBC(), SupInBC( inflow ) ],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-1", hcell_list=[(nnx,1)])
             
if with_ablation:
    b2_eastBC = AblatingBC( Twall=Twall, mdot=[0.0, mdot[2] ] )
else:
    b2_eastBC = FixedTBC( Twall )
    
if with_blk3:
    b2_northBC = AdjacentBC()
else:
    b2_northBC = ExtrapolateOutBC()
             
blk2 = SuperBlock2D(psurf=make_patch(north2, east2, south2, west2),
		     fill_condition=initial,
		     nni=nnx, nnj=nny2,
		     nbi=NBX, nbj=NBY2,
		     cf_list=cf,
		     bc_list=[b2_northBC, b2_eastBC, AdjacentBC(), SupInBC( inflow ) ],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-2", hcell_list=[(nnx,1)])
             
if with_blk3:
    blk3 = SuperBlock2D(psurf=make_patch(north3, east3, south3, west3),
                 fill_condition=initial,
                 nni=nnx, nnj=nny3,
                 nbi=NBX, nbj=NBY3,
                 cf_list=cf,
                 bc_list=[ExtrapolateOutBC(), FixedTBC( Twall ), AdjacentBC(), SupInBC( inflow ) ],
                 wc_bc_list=[NonCatalyticWBC()]*4,
                 label="BLOCK-3", hcell_list=[(nnx,1)])

identify_block_connections()

if with_blk3:
   L = f.x - a.x
   h = k.y
else:
   L = e.x - a.x
   h = j.y

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = L / u_inf * max_BLs
gdata.max_step = 1000000
gdata.dt = 1.0e-10
gdata.stringent_cfl = 1
gdata.cfl = 0.5
gdata.cfl_count = 10
gdata.print_count = 20
# gdata.fixed_time_step = True
gdata.dt_plot = L / u_inf * plot_BLs
gdata.dt_history = 1.0e-6

sketch.scales(0.3/L, 0.3/L)
sketch.origin(0.0, 0.0)
sketch.xaxis(-0.1*L,L, 0.25*L, -0.1*L)
sketch.yaxis(0.0, 1.5*h, 0.25*h, 0.0)
