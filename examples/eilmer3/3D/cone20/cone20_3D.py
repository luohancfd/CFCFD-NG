# cone20_3D.py
# 20-degree half-angle cone in 3D
# PJ, 08-Feb-2005, Sep-2008 2D
# RJG Sep-2009 3D    

job_title = "Mach 1.5 flow over a 20 degree cone."
print job_title

# We can set individual attributes of the global data object.
gdata.dimensions = 3
gdata.title = job_title
gdata.stringent_cfl = 1  # to match the old mb_cns behaviour

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)
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

patch0 = make_patch(fe, be, ab, af)
patch1 = make_patch(ed, cd, bc, be)

def revolved_vol0(r, s, t):
    global patch0
    P = patch0.eval(r, s)
    phi = t*math.pi
    y = P.y*math.cos(phi)
    z = P.y*math.sin(phi)
    return (P.x, y, z)

def revolved_vol1(r, s, t):
    global patch1
    P = patch1.eval(r, s)
    phi = t*math.pi
    y = P.y*math.cos(phi)
    z = P.y*math.sin(phi)
    return (P.x, y, z)

# Define the blocks, boundary conditions and set the discretisation.
nx0 = 5; nx1 = 30; ny = 30; nz = 30
clustery = RobertsClusterFunction(1, 0, 1.2)

blk_0 = Block3D(PyFunctionVolume(revolved_vol0),
                nni=nx0, nnj=ny, nnk=nz,
                fill_condition=initial, label="blk-0")
# FIX-ME PJ Sep-2012: blk_0 needs compatible clustering. 
blk_1 = Block3D(PyFunctionVolume(revolved_vol1),
                nni=nx1, nnj=ny, nnk=nz,
                cf_list=[None, clustery, None, clustery,
                         None, clustery, None, clustery,
                         None, None, None, None],
                fill_condition=initial, label="blk-1")
identify_block_connections()
blk_1.bc_list[WEST] = SupInBC(inflow)
blk_1.bc_list[EAST] = ExtrapolateOutBC()

# Do a little more setting of global data.
gdata.viscous_flag = 0
gdata.flux_calc = ADAPTIVE
gdata.max_time = 5.0e-3  # seconds
gdata.max_step = 3000
gdata.dt = 1.0e-7
gdata.dt_plot = 1.5e-3
gdata.dt_history = 10.0e-5
