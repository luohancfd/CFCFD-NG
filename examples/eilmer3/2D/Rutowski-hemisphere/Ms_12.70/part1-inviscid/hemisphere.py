## \file hemisphere.py
## \brief Mach 12.7 condition from Rutowski and Bershader (1964)
## \author DFP, 28-May-2014
##

from cfpylib.grid.shock_layer_surface import *

gdata.title  = "Shock heated argon flow over a 1/2 inch hemisphere"
gdata.title += "- part 1: inviscid"
print gdata.title

# axisymmetry
gdata.axisymmetric_flag = 1

# gas model
species = select_gas_model(  model = "two temperature gas", \
                           species = [ "Ar", "Ar_plus", "e_minus" ] )
gm = get_gas_model_ptr()
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()

# kinetics
set_reaction_update("../../kinetic-models/Ar-2T-chemical-reactions.lua")
set_energy_exchange_update("../../kinetic-models/Ar-2T-energy-exchange.lua")
    
# flow conditions - shock heated argon, initially at 10 Torr and 300K
T_wall = 300.0
Ms = 12.7
from cfpylib.gasdyn.cea2_gas import *
reactants = { "Ar" : 1.0, "Ar+" : 0.0, "e-" : 0.0 }
cea = Gas( reactants, onlyList=reactants.keys(), with_ions=True, \
           trace=1.0e-20 )
cea.set_pT(p=1333.3,T=300.0)
Us = cea.a * Ms
print "Us = ", Us
cea.shock_process(Us)
rho_inf = cea.rho
T_inf = [ cea.T ]*ntm
massf_inf = []
for isp,sp in enumerate(species):
    cea_sp = sp.replace("_plus","+").replace("_minus","-")
    massf_inf.append( cea.species[cea_sp] )
u_inf = Us - cea.u2

# do some calculations to get pressure and Mach number
Q = Gas_data(gm)
Q.rho = rho_inf
for itm in range(ntm):
    Q.T[itm] = T_inf[itm]
mf_sum = 0.0
for isp in range(nsp):
    Q.massf[isp] = massf_inf[isp]
    mf_sum += Q.massf[isp]
massf_inf = []
for isp in range(nsp):
    Q.massf[isp] /= mf_sum
    massf_inf.append( Q.massf[isp] )
gm.eval_thermo_state_rhoT(Q)
Q.print_values(False)
M_inf = u_inf / Q.a
p_inf = Q.p
print "M_inf = %0.2f" % ( M_inf )

# inflow and initial conditions
inflow  = FlowCondition(p=p_inf,      u=u_inf, v=0.0, T=T_inf, \
                        massf=massf_inf)
initial = FlowCondition(p=p_inf/10.0, u=0.0,   v=0.0, T=T_inf, \
                        massf=massf_inf)

# geometry
Rn = 1.27e-2
psurf, west = make_parametric_surface( 1.0, 1.0, M_inf, Rn, \
                                       axi=gdata.axisymmetric_flag )
                                       
# mesh clustering
cf_list=[ None ] * 4
                                       
# boundary conditions
bc_list=[ExtrapolateOutBC(),			# outflow
         FixedTBC(T_wall),				# surface
         SlipWallBC(),					# symmetry
         SupInBC(inflow)]				# inflow
         
# catalytic boundary conditions
wc_bc_list=[NonCatalyticWBC()]*4

blk_0 = SuperBlock2D(psurf=psurf,
             fill_condition=initial,
             nni=40, nnj=30,
             nbi=2, nbj=2,
             cf_list=cf_list,
             bc_list=bc_list,
             wc_bc_list=wc_bc_list,
             label="BLOCK-0")

identify_block_connections()

# global simulation parameters
gdata.viscous_flag = 0
gdata.viscous_delay = 0.1 * Rn / u_inf
gdata.viscous_factor_increment = 1.0e-3
gdata.diffusion_flag = 0
gdata.diffusion_model = "ConstantLewisNumber"
gdata.electric_field_work_flag = 0
gdata.reaction_time_start = 0 * Rn / u_inf
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.max_time = Rn * 10 / u_inf    # 10 body lengths
gdata.reaction_time_start = Rn * 1 / u_inf
gdata.max_step = 230000
gdata.dt = 1.0e-10
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1 / u_inf    # 10 solutions
gdata.cfl = 1.0
gdata.cfl_count = 10
gdata.print_count = 20

sketch.scales(0.03/Rn, 0.03/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-2.0e-2, 0.0, 0.5e-2, -0.3e-2)
sketch.yaxis(0.0, 4.0e-2, 1.0e-2, -0.3e-2)
