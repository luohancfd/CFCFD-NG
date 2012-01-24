## \file cone20.py
## \brief Simple job-specification file for e3prep.py
## \author PJ, 08-Feb-2005
##
## 15-Sep-2008 -- simplified version for Elmer3
##
## We have set this file up very much like the cone20.sit file
## so that users may more-easily see the correspondence between
## the Tcl and Python elements.

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

# Define flow conditions.
# We will get the initial flow condition from an old solution.
initial = ExistingSolution('cone20_old', ".", 2, 1, gdata.dimensions)
inflow  = FlowCondition(p=95.84e3, u=1000.0, v=0.0, T=1103.0)

# Set up two quadrilaterals in the (x,y)-plane be first defining
# the corner nodes, then the lines between those corners and then
# the boundary elements for the blocks.
# The labelling is not significant; it is just to make the SVG picture
# look the same as that produced by the Tcl scriptit program.
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
blk_0 = Block2D(make_patch(fe, be, ab, af), nni=nx0, nnj=ny,
                fill_condition=initial, label="BLOCK-0")
blk_1 = Block2D(make_patch(ed, cd, bc, be, "AO"), nni=nx1, nnj=ny,
                fill_condition=initial, label="BLOCK-1",
                hcell_list=[(9,0)], xforce_list=[0,0,1,0])
identify_block_connections()
blk_0.bc_list[WEST] = SupInBC(inflow) # one way to set a BC
blk_1.set_BC(EAST, EXTRAPOLATE_OUT)   # another way

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 3.5e-3  # seconds
gdata.max_step = 3000
gdata.dt = 1.0e-6
gdata.dt_plot = 1.5e-3
gdata.dt_history = 10.0e-5

sketch.xaxis(0.0, 1.0, 0.2, -0.05)
sketch.yaxis(0.0, 1.0, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)
