## \file blunted_cone.py
## \brief Job-specification file for 5 degree cone in
## \brief Mach 7.4 flow (Ames 3.5 foot blowdown tunnel)
## \author DFP et al, 30-October-2009
##
## Ref: Marvin, J.G. and Akin, C.M.
##      Combined Effects of Mass Addition and Nose Bluntness 
##      on Boundary-Layer Transition
##      AIAA Journal Vol. 8 No. 5 pp 857-863 May 1970 

import sys
from math import *
from gaspy import *

job_title = "5 degree cone in Mach 7.4 Flow."
print job_title

# We can set individual attributes of the global data object.
gdata.title = job_title
gdata.axisymmetric_flag = 1

gas_model = "IG"
transfer_solution = False
with_ablation = False
with_porous_section = True
injected_gas = 'air' 
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
M_inf = 7.4
p_total_psi = 600.0     # total pressure in psi
T_total_R = 1500.0      # total temperature in degrees Rankine
gamma_air = 1.4
R_air = 287.0
T_inf = ( T_total_R * 0.555555556 ) / ( 1.0 + ( gamma_air - 1.0 ) / 2.0 * M_inf**2 )
u_inf = M_inf * sqrt( gamma_air * R_air * T_inf )
p_inf = ( p_total_psi * 6894.75729 ) / ( 1.0 + ( gamma_air - 1.0 ) / 2.0 * M_inf**2 )**(gamma_air/(gamma_air-1.0))
Twall = 0.37 * ( T_total_R * 0.555555556 )
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
Re = Q.rho * u_inf * ( 3.75 * 0.0254 + 16.25 * 0.0254 ) / Q.mu
Q.print_values()
print "Reynolds number = ", Re

inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=massf_inf)

if transfer_solution:
    initial = ExistingSolution(rootName="cone", solutionWorkDir="../mdot_0.000", nblock=14, tindx=9999)
else:
    initial = FlowCondition(p=p_inf/10.0, u=0.0, v=0.0, T=T_inf, massf=massf_inf )

# Set up the geometry
L_impermiable = 3.75 * 0.0254
L_porous = 16.25 * 0.0254
theta = 5.0 * pi / 180.0	# 5 degree cone in radians
x_total = ( L_impermiable + L_porous ) * cos( theta )
y_total = ( L_impermiable + L_porous ) * sin( theta )
f_height_start = 0.05
f_height_end = 2.4
alpha = atan( ( f_height_end - f_height_start ) * y_total / x_total )

# inflow block points
o = Node( 0.0, 0.0, label="o" )
a = Node( - 0.01 * x_total, 0.0, label="a" )
b = Node( a.x, y_total * f_height_start, label="b" )
c = Node( 0.0, b.y, label="c" )
# impermiable cone block points
d = Node ( L_impermiable, L_impermiable*sin(alpha) + f_height_start*y_total, label="d" )
g = Node( d.x, d.x*sin(theta), label="g" )
# porous cone block points
e = Node ( x_total, x_total*sin(alpha), label="e" )
f = Node( x_total, y_total, label="f" )

# inflow block edges
north0 = Line(b,c)
east0  = Line(o,c)
south0 = Line(a,o)
west0  = Line(a,b)
# impermiable cone block edges
north1 = Line(c,d)
east1 = Line(g,d)
south1 = Line(o,g)
west1 = east0
# porous cone block edges
north2 = Line(d,e)
east2 = Line(f,e)
south2 = Line(g,f)
west2 = east1

# Define the blocks, boundary conditions and set the discretisation.
nnx_cone = int(cell_scale*300)
nnx0 = int(4 * cell_scale)
nnx1 = int(L_impermiable/x_total*nnx_cone)
nnx2 = nnx_cone - nnx1
nny = int(40*cell_scale)
NBX0 = 1		# super-blocking
if with_porous_section: NBX1 = 1
else: NBX1 = 4
NBX2 = 4
NBY = 3

print "block set 0: cells-per-sub-block = ", int(nnx0*nny/float(NBX0*NBY))
print "block set 1: cells-per-sub-block = ", int(nnx1*nny/float(NBX1*NBY))
print "block set 2: cells-per-sub-block = ", int(nnx2*nny/float(NBX2*NBY))

betaNS = 1.1               # clustering
betaSN = 1.1
betaWE = 0.0
betaEW = 0.0

cluster_functions = [RobertsClusterFunction(0, 1, betaEW),
                     RobertsClusterFunction(1, 0, betaNS),
                     RobertsClusterFunction(0, 1, betaWE),
                     RobertsClusterFunction(1, 0, betaSN)]

blk_0 = SuperBlock2D(psurf=make_patch(north0, east0, south0, west0),
		     fill_condition=initial,
		     nni=nnx0, nnj=nny,
		     nbi=NBX0, nbj=NBY,
		     cf_list=cluster_functions,
		     bc_list=[SupInBC( inflow ), AdjacentBC(), SlipWallBC(), SupInBC( inflow )],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-0", hcell_list=[(nnx0,1)])
             
if with_porous_section:
   blk1_eastBC = AdjacentBC()
else:
   blk1_eastBC = ExtrapolateOutBC()
             
blk_1 = SuperBlock2D(psurf=make_patch(north1, east1, south1, west1),
		     fill_condition=initial,
		     nni=nnx1, nnj=nny,
		     nbi=NBX1, nbj=NBY,
		     cf_list=cluster_functions,
		     bc_list=[SupInBC( inflow ), blk1_eastBC, FixedTBC(Twall), AdjacentBC()],
		     wc_bc_list=[NonCatalyticWBC()]*4,
		     label="BLOCK-1", hcell_list=[(nnx1,1)])
             
if with_porous_section:
    if with_ablation:
       mdot = [ 0.0 ] * nsp 
       mdot[-1] = - rho_inf * u_inf * mdot_ratio
       blk2_southBC = AblatingBC( Twall=Twall, mdot=mdot )
    else:
       blk2_southBC = FixedTBC(Twall)
       
    blk_2 = SuperBlock2D(psurf=make_patch(north2, east2, south2, west2),
                 fill_condition=initial,
                 nni=nnx2, nnj=nny,
                 nbi=NBX2, nbj=NBY,
                 cf_list=cluster_functions,
                 bc_list=[SupInBC( inflow ), ExtrapolateOutBC(), blk2_southBC, AdjacentBC()],
                 wc_bc_list=[NonCatalyticWBC()]*4,
                 label="BLOCK-2", hcell_list=[(nnx2,1)])

identify_block_connections()

if with_porous_section:
   L = x_total
else:
   L = L_impermiable

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

sketch.scales(0.3/x_total, 0.3/x_total)
sketch.origin(0.0, 0.0)
sketch.xaxis(-0.1*x_total,x_total, 0.25*x_total, -0.1*x_total)
sketch.yaxis(0.0, 1.5*y_total, 0.25*y_total, 0.0)
