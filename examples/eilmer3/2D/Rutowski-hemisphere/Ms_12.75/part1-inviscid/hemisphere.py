## \file hemisphere.py
## \brief 5.6km/s nitrogen flow over a 45mm cylinder
## \author DFP, 23-Nov-2009
##

from cfpylib.grid.shock_layer_surface import *

job_title = "Shock heated Argon flow over a 12.7mm diameter hemisphere (Rutowski and Bershader, 1964)"
print job_title

# define simulation control parameters
transfer_solution = False
fit2shock = False
inflow_source = "poshax3"
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

# Define flow conditions - shock heated argon, initially at 10 Torr and 300K
T_wall = 300.0
if inflow_source=="cea2":
    from cfpylib.gasdyn.cea2_gas import *
    M_inf = 12.75
    reactants = { "Ar" : 1.0, "Ar+" : 0.0, "e-" : 0.0 }
    cea = Gas( reactants, onlyList=reactants.keys(), with_ions=True, trace=1.0e-20 )
    cea.set_pT(p=1333.3,T=300.0)
    Us = cea.a * M_inf
    print "Us = ", Us
    cea.shock_process(Us)
    rho_inf = cea.rho
    T_inf = [ cea.T ]*ntm
    massf_inf = []
    for isp,sp in enumerate(species):
        massf_inf.append( cea.species[sp.replace("_plus","+").replace("_minus","-")] )
    # FIXME: get velocity from CEA2
    u_inf = 3360.44
elif inflow_source=="poshax3":
    # see 'freestream' directory for the source of these values
    T_inf = [ 11650.6025787 ]*ntm
    rho_inf = 0.116289607429
    massf_inf = [ 0.949450404198, 0.0505489016373, 6.94164723533e-07 ]
    u_inf = 12.75*sqrt(1.6667*PC_R_u/gm.molecular_weight(0)*300.0) - 755.248893203

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
    initial = ExistingSolution(rootName="hemisphere", solutionWorkDir="../inviscid-2T/", nblock=1, tindx=9999)
else:
    initial = FlowCondition(p=p_inf/10.0, u=0.0, v=0.0, T=T_inf, massf=massf_inf)

# make the geometry
Rn = 1.27e-2
beta0 = 1.1; dx0 = 5.0e-1; dx1 = 5.0e-2; gamma = 0.2
beta1 = 1.0
cf_list = [ None ] * 4
nnx = 40; nny = 30
if fit2shock:
    print "WARNING: the shock fitting procedure takes a long time as the Billig function"
    print "         is difficult to solve at this Mach number."
    shock = fit_billig2shock( initial, gdata.axisymmetric_flag, M_inf, Rn, None )
    psurf = make_parametric_surface( M_inf=M_inf, R=Rn, axi=gdata.axisymmetric_flag, east=None, shock=shock, f_s=1.0/(1.0-gamma) )
else:  
    bx_scale = 1.0; by_scale = 1.0
    psurf = make_parametric_surface( bx_scale, by_scale, M_inf, Rn, axi=gdata.axisymmetric_flag )
    

blk_0 = SuperBlock2D(psurf=psurf,
		     fill_condition=initial,
		     nni=nnx, nnj=nny,
		     nbi=2, nbj=2,
		     cf_list=cf_list,
		     bc_list=[ExtrapolateOutBC(), FixedTBC(T_wall), SlipWallBC(), SupInBC(inflow)],
                     # wc_bc_list=[NonCatalyticWBC(),SuperCatalyticWBC([1.0,0.0,0.0]),NonCatalyticWBC(),NonCatalyticWBC()],
                     wc_bc_list=[NonCatalyticWBC(),NonCatalyticWBC(),NonCatalyticWBC(),NonCatalyticWBC()],
		     label="BLOCK-0", hcell_list=[(nnx,1)])
                 
identify_block_connections()

# Do a little more setting of global data
gdata.viscous_flag = 0
gdata.viscous_delay = 0.1 * Rn / u_inf
gdata.viscous_factor_increment = 1.0e-3
gdata.diffusion_flag = 0
gdata.diffusion_model = "ConstantLewisNumber"
gdata.reaction_time_start = 0 * Rn / u_inf
gdata.flux_calc = ADAPTIVE
gdata.max_time = Rn * 10 / u_inf    # 10 body lengths
gdata.reaction_time_start = Rn * 1 / u_inf
gdata.max_step = 230000
gdata.dt = 1.0e-10
gdata.stringent_cfl = 1
gdata.dt_plot = Rn * 1 / u_inf    # 10 solutions
gdata.cfl = 0.5
gdata.cfl_count = 10
gdata.print_count = 20

# flux calc options
# gdata.compression_tolerance = -0.5

sketch.scales(0.3/Rn, 0.3/Rn)
sketch.origin(0.0, 0.0)
sketch.xaxis(-1.5*Rn, 0.0, 0.25*Rn, -0.1*Rn)
sketch.yaxis(0.0, 1.5*Rn, 0.25*Rn, 0.0)
