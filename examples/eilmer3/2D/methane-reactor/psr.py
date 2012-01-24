# psr.py
# Methane combustion in a perfectly-stirred reactor of fixed volume.
#
# This example demonstrates the use of Brendan's implementation of
# the GRI-Mech3.0 reaction scheme.
#
# PJ, 05-Sep-2011

gdata.title = "Methane combustion"
# We have built the gas model with the pieces that brendan provided,
# outside this script.
species_list = select_gas_model(fname="thermally-perfect-grimech30.lua")
print "species_list=", species_list
set_reaction_scheme("grimech30.lua",reacting_flag=1)

# Initially, the gas is a stoichiometric mix at a combustible condition.
gmodel = get_gas_model_ptr()
molef = {'CH4':1.0, 'O2':2.0, 'N2':7.52}
massf = gmodel.to_massf(molef)
initial = FlowCondition(p=1.01325e5, u=0.0, v=0.0, T=2000.0, massf=massf)
print "initial=", initial

# Geometry is a simple square and one cell is nominated as the point to
# record the history data.
p0 = Node(0.0,  0.0);  p1 = Node(0.01, 0.0)
p2 = Node(0.01, 0.01); p3 = Node(0.0,  0.01)
domain = make_patch(Line(p3,p2), Line(p1,p2), Line(p0,p1), Line(p0,p3))
blk = Block2D(domain, nni=2, nnj=2, fill_condition=initial, hcell_list=[(0,0)])

# Simulation parameters
gdata.flux_calc = ADAPTIVE
gdata.max_time = 400.0e-6
gdata.max_step = 4000
gdata.dt = 1.0e-8
gdata.cfl = 0.1  # force small steps so that we get the historical record
gdata.dt_plot = 20.0e-6
gdata.dt_history = 1.0e-6

