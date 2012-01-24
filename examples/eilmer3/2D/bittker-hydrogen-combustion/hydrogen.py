## \file hydrogen.py
## \brief Eilmer3 job-spec file of a stream-tube with hydrogen combustion.
## \author Fabian Zander, 26th Oct 2010
# Edited by PJ to make it simpler for User Guide, 21-Jun-2011
#
# This example is trying to match test case 3 from the Bittker-Scullin code.
# See the report NASA TN D-6586.

gdata.title = "Hydrogen-Air Combustion Validation"
print gdata.title
gdata.dimensions = 2
gdata.axisymmetric_flag = 1
gdata.stringent_cfl = 1

# Hydrogen-air combustion model.
select_gas_model(model='thermally perfect gas',
                 species=['O', 'O2', 'N2', 'H', 'H2', 'H2O', 'HO2', 'OH', 'H2O2'])
set_reaction_scheme("Bittker_Scullin.lua", reacting_flag=1)
gmodel = get_gas_model_ptr()

# Initial mixture taken from Bittker-Scullin example on page 85 in mole-fraction.
molef = {'O2':0.1480, 'N2':0.5562, 'H2':0.2958}
inflow = FlowCondition(p=96.87e3, u=4551.73, v=0.0, T=1559.00, massf=gmodel.to_massf(molef))

# Setting up the geometry as a rectangle in the x,y-plane.
#   ^ y
#   | 
#   |
#   d----------c
#   |          |
#   a----------b --> x
d = Node(0.0, 0.01, label="d"); c = Node(0.1, 0.01, label="c")
a = Node(0.0, 0.00, label="a"); b = Node(0.1, 0.00, label="b")
blk = SuperBlock2D(make_patch(Line(d,c), Line(b,c), Line(a,b), Line(a,d)), 
                   nni=2000, nnj=2, nbi=200, nbj=1,
                   bc_list=[SlipWallBC(),ExtrapolateOutBC(),SlipWallBC(),SupInBC(inflow)],
                   fill_condition=inflow, label="Block")

gdata.sequence_blocks = 1
gdata.viscous_flag = 0
gdata.flux_calc = ADAPTIVE
gdata.t_order = 2
gdata.max_time = 5.0e-5
gdata.max_step = 1000000
gdata.dt = 1.0e-10
gdata.dt_plot = 5.0e-6
gdata.dt_history = 10.0e-5

sketch.xaxis(0.0, 0.10, 0.02, -0.005)
sketch.yaxis(0.0, 0.02, 0.02, -0.005)
sketch.window(0.0, 0.0, 0.1, 0.1, 0.05, 0.05, 0.17, 0.17)
