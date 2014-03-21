## \file hemisphere.py
## \brief 5.6km/s nitrogen flow over a 45mm cylinder
## \author DFP, 23-Nov-2009
##

from cfpylib.grid.shock_layer_surface import *

job_title = "Shock heated Argon flow over a 12.7mm diameter hemisphere (Rutowski and Bershader, 1964)"
print job_title

# define simulation control parameters
transfer_solution = True
fit2shock = True
inflow_source = "cea2"
gmodel = "two temperature gas"
species = [ "Ar", "Ar_plus", "e_minus" ]
rmodel = "Hoffert67-Ar-2T.lua"
eemodel = "Ar-TE.lua"

# We can set individual attributes of the global data object.
gdata.title = job_title
gdata.axisymmetric_flag = 1

# gas model
species = select_gas_model(model=gmodel, species=species)
gm = get_gas_model_ptr()
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()

# kinetics
if nsp>1:
    set_reaction_update(rmodel)
if ntm>1:
    set_energy_exchange_update(eemodel)

# radiation model
# select_radiation_model("rad-model.lua",-1)

# Define flow conditions - shock heated argon, initially at 10 Torr and 300K
T_wall = 300.0
Ms = 12.7
if inflow_source=="cea2":
    from cfpylib.gasdyn.cea2_gas import *
    reactants = { "Ar" : 1.0, "Ar+" : 0.0, "e-" : 0.0 }
    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=True, trace=1.0e-20 )
    cea.set_pT(p=1333.3,T=300.0)
    Us = cea.a * Ms
    print "Us = ", Us
    cea.shock_process(Us)
    rho_inf = cea.rho
    T_inf = [ cea.T ]*ntm
    massf_inf = []
    for isp,sp in enumerate(species):
        massf_inf.append( cea.species[sp.replace("_plus","+").replace("_minus","-")] )
    u_inf = Us - cea.u2
elif inflow_source=="poshax3":
    raise RuntimeError, "poshax3 inflow not implemented yet."
    T_inf = []*ntm
    rho_inf = 0.0 
    massf_inf = [ 0.0 ] * 3
    u_inf = Ms*sqrt(1.6667*PC_R_u/gm.molecular_weight(0)*300.0) - 0.0

# do some calculations to get pressure, Mach number and total mass-flux
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
print "p_inf = %0.2f" % ( p_inf )
print "rho_inf = %0.2f" % ( rho_inf )
print "T_inf = %0.2f" % ( T_inf[0] )
print "u_inf = %0.2f" % ( u_inf )
print "M_inf = %0.2f" % ( M_inf )

# define the inflow and initial conditions
inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=massf_inf)
if transfer_solution:
    initial = ExistingSolution(rootName="hemisphere", solutionWorkDir="../part1-inviscid/", nblock=4, tindx=10)
else:
    initial = FlowCondition(p=p_inf/10.0, u=0.0, v=0.0, T=T_inf, massf=massf_inf)

# make the geometry
Rn = 1.27e-2
beta0 = 1.1; dx0 = 5.0e-1; dx1 = 5.0e-2; gamma = 0.2
beta1 = 1.0
cf_list = [BHRCF(beta0,dx0,dx1,gamma),
           RCF(0,1,beta1),
           BHRCF(beta0,dx0,dx1,gamma),
           RCF(0,1,beta1)]
nnx = 60; nny = 45
if fit2shock:
    print "WARNING: the shock fitting procedure takes a long time as the Billig function"
    print "         is difficult to solve at this Mach number."
    shock, nodes = fit_billig2shock( initial, gdata.axisymmetric_flag, M_inf, Rn, None, show_plot=False )
    psurf, west = make_parametric_surface( M_inf=M_inf, R=Rn, axi=gdata.axisymmetric_flag, east=None, shock=shock, f_s=1.0/(1.0-gamma) )
else:  
    bx_scale = 1.0; by_scale = 1.0
    psurf, west = make_parametric_surface( bx_scale, by_scale, M_inf, Rn, axi=gdata.axisymmetric_flag )
    

blk_0 = SuperBlock2D(psurf=psurf,
		     fill_condition=initial,
		     nni=nnx, nnj=nny,
		     nbi=2, nbj=2,
		     cf_list=cf_list,
		     bc_list=[ExtrapolateOutBC(), FixedTBC(T_wall), SlipWallBC(), SupInBC(inflow)],
                     wc_bc_list=[NonCatalyticWBC(),SuperCatalyticWBC([1.0,0.0,0.0]),NonCatalyticWBC(),NonCatalyticWBC()],
		     label="BLOCK-0", hcell_list=[(nnx,1)])
                 
identify_block_connections()

# Do a little more setting of global data
gdata.viscous_flag = 1
gdata.viscous_delay = 0.001 * Rn / u_inf
gdata.viscous_factor_increment = 1.0e-4
# NOTE: diffusion is currently not working and is turned off
gdata.diffusion_flag = 0
gdata.diffusion_delay = 0.001 * Rn / u_inf
gdata.diffusion_factor_increment = 1.0e-4
gdata.diffusion_model = "Ramshaw-Chang"
# NOTE: if an ambipolar diffusion model is being used, the electric field work term should be included
gdata.electric_field_work_flag = gdata.diffusion_flag
gdata.reaction_time_start = 0 * Rn / u_inf
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.max_time = Rn * 1 / u_inf    # 1 body length
gdata.max_step = 2300000
gdata.dt = 1.0e-10
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1 / u_inf    # 1 solution
# NOTE: the CFL number can be increased to 0.5 after the viscous terms have been added
gdata.cfl = 1.0e-1
gdata.cfl_count = 1
gdata.print_count = 10

# flux calc options
# gdata.compression_tolerance = -0.5

sketch.scales(0.03/Rn, 0.03/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-2.0e-2, 0.0, 0.5e-2, -0.3e-2)
sketch.yaxis(0.0, 4.0e-2, 1.0e-2, -0.3e-2)
