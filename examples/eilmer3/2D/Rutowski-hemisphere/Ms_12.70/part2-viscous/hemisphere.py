## \file hemisphere.py
## \brief Mach 12.7 condition from Rutowski and Bershader (1964)
## \author DFP, 28-May-2014
##

from cfpylib.grid.shock_layer_surface import *

gdata.title  = "Shock heated argon flow over a 1/2 inch hemisphere"
gdata.title += "- part 2: viscous"
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

# inflow and initial conditions (continuation from part 1)
inflow  = FlowCondition(p=p_inf,      u=u_inf, v=0.0, T=T_inf, \
                        massf=massf_inf)
initial = ExistingSolution(rootName="hemisphere", \
                           solutionWorkDir="../part1-inviscid/", \
                           nblock=4, tindx=10)

# geometry
Rn = 1.27e-2
gamma = 0.2
print "WARNING: the shock fitting procedure takes a long time as the"
print "         Billig function is difficult to solve at this Mach"
print "         number."
shock, nodes = fit_billig2shock( initial, gdata.axisymmetric_flag, \
                                 M_inf, Rn, None, show_plot=False )
psurf, west = make_parametric_surface( M_inf=M_inf, R=Rn, \
                                       axi=gdata.axisymmetric_flag, \
                                       east=None, shock=shock, \
                                       f_s=1.0/(1.0-gamma) )

# boundary conditions
bc_list=[ExtrapolateOutBC(),              # outflow
         FixedTBC(T_wall),                # surface
         SlipWallBC(),                    # symmetry
         SupInBC(inflow)]                 # inflow
          
# catalycity boundary conditions
wc_bc_list=[NonCatalyticWBC(),                # outflow
           SuperCatalyticWBC([1.0,0.0,0.0]),  # surface
           NonCatalyticWBC(),                 # symmetry
           NonCatalyticWBC()]                 # inflow
                                       
# mesh clustering
beta0 = 1.1; dx0 = 5.0e-1; dx1 = 5.0e-2
beta1 = 1.0
cf_list = [BHRCF(beta0,dx0,dx1,gamma),    # outflow
           RCF(0,1,beta1),                # surface
           BHRCF(beta0,dx0,dx1,gamma),    # symmetry
           RCF(0,1,beta1)]                # inflow         

# computation domain
blk_0 = SuperBlock2D(psurf=psurf,
             fill_condition=initial,
             nni=60, nnj=45,
             nbi=2, nbj=2,
             cf_list=cf_list,
             bc_list=bc_list,
             wc_bc_list=wc_bc_list,
             label="BLOCK-0")

identify_block_connections()

# global simulation parameters
gdata.viscous_flag = 1
gdata.viscous_delay = 0.001 * Rn / u_inf
gdata.viscous_factor_increment = 1.0e-4
# NOTE: diffusion is currently turned off
gdata.diffusion_flag = 0
gdata.diffusion_delay = 0.001 * Rn / u_inf
gdata.diffusion_factor_increment = 1.0e-4
gdata.diffusion_model = "Ramshaw-Chang"
# NOTE: if an ambipolar diffusion model is being used, the electric 
#       field work term should be included
gdata.electric_field_work_flag = gdata.diffusion_flag
gdata.reaction_time_start = 0 * Rn / u_inf
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.max_time = Rn * 1 / u_inf    # 1 body length
gdata.max_step = 2300000
gdata.dt = 1.0e-10
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1 / u_inf    # 1 solution
# NOTE: the CFL number can be increased to 0.5 after the viscous terms 
#       have been added
gdata.cfl = 1.0e-1
gdata.cfl_count = 1
gdata.print_count = 10

sketch.scales(0.03/Rn, 0.03/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-2.0e-2, 0.0, 0.5e-2, -0.3e-2)
sketch.yaxis(0.0, 4.0e-2, 1.0e-2, -0.3e-2)
