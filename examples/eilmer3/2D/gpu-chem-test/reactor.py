# Author: Rowan J. Gollan and Kyle Damm
# Date: 29-Jan-2015

gdata.title = "H2/O2 combustion reactor."
gdata.dimensions = 2

# Set gas model and chemistry scheme
select_gas_model(model='thermally perfect gas',
                 species=['O2', 'H2', 'H2O', 'O', 'H', 'OH'])
set_reaction_scheme("Evans_Schexnayder.lua")
gmodel = get_gas_model_ptr()

# Set the intial condition
init_moles = {'H2':2.0, 'O2':1.0}
tot_moles = sum(init_moles.values())
init_molef = { k:v/tot_moles for k,v in init_moles.iteritems() }
init_massf = gmodel.to_massf(init_molef)
init_p = 101.325e3
init_T = 1000.0
initial = FlowCondition(p=init_p, T=[init_T], massf=init_massf)

# Create a box, 0.1 x 0.1 m
a = Node(0.0, 0.0)
b = Node(0.0, 0.1)
c = Node(0.1, 0.0)
d = Node(0.1, 0.1)
blk0 = Block2D(make_patch(Line(b, d), Line(c, d), Line(a, c), Line(a, b)),
               nni=2, nnj=2,
               fill_condition=initial,
               hcell_list=[(0,0)])

# Set the simulation parameters.
gdata.fixed_time_step = True
gdata.dt = 1.0e-7
gdata.max_time = 5.0e-4
gdata.max_step = 2*(gdata.max_time/gdata.dt)
gdata.dt_history = 1.0e-7
gdata.gasdynamic_update_scheme = 'euler'





