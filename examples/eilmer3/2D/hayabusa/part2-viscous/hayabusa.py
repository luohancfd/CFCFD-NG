## \file hayabusa.py
## \brief Simulating the JAXA Hayabusa sample return capsule
## \author DFP, 30-Oct-2012
## 
## Part2: Viscous solution on a finer grid
##             Condition from Dan's thesis

from math import *
from cfpylib.grid.shock_layer_surface import *

job_title = "JAXA Hayabusa sample return capsule."
print job_title

gdata.title = job_title
gdata.axisymmetric_flag = 1

fit2shock = False

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
rho_inf = 1.73e-4
T_inf = 230.0
u_inf = 9.679e3
massf_inf = [ 0.0 ] * gm.get_number_of_species()
massf_inf[species.index('N2')] = 0.767
massf_inf[species.index('O2')] = 0.233
fixed_T_wall = 3000.0

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

# make the computational domain
x_limit = c.x; gamma = 0.2
if fit2shock:
    shock, nodes = fit_billig2shock( initial, gdata.axisymmetric_flag, M_inf, Rn, body )
    psurf = make_parametric_surface( M_inf=M_inf, R=Rn, axi=gdata.axisymmetric_flag, east=body, shock=shock, f_s=1.0/(1.0-gamma) )
else:  
    bx_scale = 0.95; by_scale = 0.95
    psurf, west = make_parametric_surface( bx_scale, by_scale, M_inf, Rn, axi=gdata.axisymmetric_flag )

#
# 4. Define the blocks, boundary conditions and set the discretisation
#
nnx = 60; nny=40
nbx = 4; nby = 4

# clustering at shock and boundary layer [ outflow, surface, axis, inflow ]
beta0 = 1.1; dx0 = 5.0e-1; dx1 = 2.0e-2 # gamma is previously defined
beta1 = 1.2
cf_list = [BHRCF(beta0,dx0,dx1,gamma),
           RCF(1,0,beta1),
           BHRCF(beta0,dx0,dx1,gamma),
           RCF(1,0,beta1)]

# boundary conditions [ outflow, surface, axis, inflow ]
if fixed_T_wall:
    eastBC = FixedTBC(fixed_T_wall)
else:
    eastBC = SurfaceEnergyBalanceBC(0.9)
bc_list=[ ExtrapolateOutBC(),
          eastBC,
          SlipWallBC(),
          SupInBC(inflow)]

# catalycity [ outflow, surface, axis, inflow ]
wc_bc_list = [NonCatalyticWBC(),
              SuperCatalyticWBC(massf_inf),
              NonCatalyticWBC(),
              NonCatalyticWBC()]

blk_0 = SuperBlock2D(psurf=psurf,
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
gdata.viscous_factor_increment = 1.0e-5
gdata.diffusion_flag = 1
gdata.diffusion_model = "ConstantLewisNumber"
gdata.electric_field_work_flag = 0   # should only be applied with ambipolar diffusion
gdata.flux_calc = ADAPTIVE
gdata.max_time = Rn * 5 / u_inf      # 5 body lengths
gdata.max_step = 230000
gdata.dt = 1.0e-10
gdata.reaction_time_start = 0
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1 / u_inf       # 5 solutions
gdata.cfl = 0.25
gdata.cfl_count = 1
gdata.print_count = 1

#
# 6. svg sketch parameters
#
sketch.scales(0.3/Rn, 0.3/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-1.5*Rn, 0.0, 0.25*Rn, -0.1*Rn)
sketch.yaxis(0.0, 1.5*Rn, 0.25*Rn, 0.0)
