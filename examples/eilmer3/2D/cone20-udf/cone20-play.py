## \file cone20.py
## \brief Test job-specification file for e3prep.py
## \author PJ, 08-Feb-2005
##
## 15-Jan-2008 -- demonstrate user-defined boundary conditions
## 28-Apr-2009 -- finished implementation of AdjacentPlusUDF BC
##                so Dan can have a go at x-tube simulation with
##                slowly-opening diaphragms.
##
## This example contains lots of trial bits that are commented out
## because it is a test example for PJ to play, not a teaching example.

job_title = "Mach 1.5 flow over a 20 degree cone."
print job_title

# We can set individual attributes of the global data object.
gdata.dimensions = 2
gdata.title = job_title
gdata.case_id = 5
gdata.axisymmetric_flag = 1
gdata.stringent_cfl = 1  # to match the old mb_cns behaviour

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)
inflow  = FlowCondition(p=95.84e3, u=1000.0, v=0.0, T=1103.0)

# Set up two quadrilaterals in the (x,y)-plane be first defining
# the corner nodes, then the lines between those corners and then
# the boundary elements for the blocks.
# The labelling is not significant; it is just to make the MetaPost
# picture look the same as that produced by the Tcl scriptit program.
a = Node(0.0, 0.0, label="A")
b = Node(0.2, 0.0, label="B")
c = Node(1.0, 0.29118, label="C")
d = Node(1.0, 1.0, label="D")
e = Node(0.2, 1.0, label="E")
f = Node(0.0, 1.0, label="F")

ab = Line(a, b); bc = Line(b, c) # lower boundary including cone surface
fe = Line(f, e); ed = Line(e, d) # upper boundary
af = Line(a, f); be = Line(b, e); cd = Line(c, d) # vertical lines

# Define the blocks, boundary conditions and set the discretisation.
nx0 = 10; nx1 = 30; ny = 40
# help()
# help(make_patch)
# help(Block2D)
blk_0 = Block2D(make_patch(fe, be, ab, af), nni=nx0, nnj=ny,
                fill_condition=initial, label="BLOCK-1")
blk_1 = Block2D(make_patch(ed, cd, bc, be, "AO"), nni=nx1, nnj=ny,
                fill_condition=initial, label="BLOCK-1",
                hcell_list=[(9,0)], xforce_list=[0,0,1,0])
# --------------------------------------------------------------
# identify_block_connections()
# blk_0.bc_list[WEST] = SupInBC(inflow) # one way to set a BC
# blk_0.set_BC(WEST, TRANSIENT_UNI, filename="my_transient.dat")
# blk_1.set_BC(EAST, EXTRAPOLATE_OUT)   # another way
# --------------------------------------------------------------
connect_blocks_2D(blk_0, EAST, blk_1, WEST, with_udf=1, 
                  filename="diaphragm-test.lua", is_wall=0, use_udf_flux=0)
blk_0.set_BC(WEST, USER_DEFINED, filename="udf-supersonic-in.lua", use_udf_flux=1)
blk_0.set_BC(SOUTH, USER_DEFINED, filename="udf-slip-wall.lua", is_wall=1)
blk_1.set_BC(EAST, USER_DEFINED, filename="udf-extrapolate-out.lua")
blk_1.set_BC(SOUTH, USER_DEFINED, filename="udf-slip-wall.lua", is_wall=1)
# --------------------------------------------------------------

# Do a little more setting of global data.
gdata.udf_file = "udf-process.lua"
gdata.udf_source_vector_flag = 0
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
# gdata.max_time = 5.0e-3  # seconds
gdata.max_time = 2.0e-3  # seconds
gdata.max_step = 3000
gdata.dt = 1.0e-6
# gdata.print_count = 50
# gdata.cfl_count = 5
# gdata.dt = 4.0e-6
# gdata.fixed_time_step = True
gdata.dt_plot = 1.5e-3
gdata.dt_history = 10.0e-5
# -------------------------------------------------
# gdata.heat_time_start = 0.5e-3
# gdata.heat_time_stop = 2.5e-3
# gdata.heat_factor_increment = 0.01
# HeatZone(5.0e9, Vector(0.4,0.5), Vector(0.2,0.3))
# -------------------------------------------------
# ReactionZone(Vector(0.5,0.5), Vector(0.7,0.7))
# ReactionZone(Vector(0.6,0.6), Vector(0.8,0.8))
# gdata.reaction_time_start = 1.0e-3
# -------------------------------------------------
HistoryLocation(1.0, 2.0, i_offset=-2, j_offset=1, label="here")

sketch.xaxis(0.0, 1.0, 0.2, -0.05)
sketch.yaxis(0.0, 1.0, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)
