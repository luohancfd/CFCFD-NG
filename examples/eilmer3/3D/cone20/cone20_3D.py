# cone20_3D.py
# 20-degree half-angle cone in 3D
# PJ, 08-Feb-2005, Sep-2008 2D
# RJG Sep-2009 3D
# PJ revisited Sep-2012 (seems to be a Spring thing...)  

job_title = "Mach 1.5 flow over a 20 degree cone."
print job_title

# We can set individual attributes of the global data object.
gdata.dimensions = 3
gdata.title = job_title
gdata.stringent_cfl = 0  # to match the old mb_cns behaviour

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)
inflow  = FlowCondition(p=95.84e3, u=1000.0, v=0.0, T=1103.0)

# Set up two quadrilaterals in the (x,y)-plane be first defining
# the corner nodes, then the lines between those corners and then
# the boundary elements for the blocks.
a = Node(0.0, 0.0, label="A")
b = Node(0.2, 0.0, label="B")
c = Node(1.0, 0.29118, label="C")
d = Node(1.0, 1.0, label="D")
e = Node(0.2, 1.0, label="E")
f = Node(0.0, 1.0, label="F")

ab = Line(a, b); bc = Line(b, c) # inner edge including cone surface
fe = Line(f, e); ed = Line(e, d) # outer edge/boundary
af = Line(a, f); be = Line(b, e); cd = Line(c, d) # vertical lines

patch0 = make_patch(fe, be, ab, af)
patch1 = make_patch(ed, cd, bc, be)

# Rotate each patch about the x-axis to get volumes for the gas flow.
# For vol0, the entire SOUTH boundary is collapsed to the x-axis.
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
# We're going to slice the two regions into three pieces each in the
# circumferential (k-index) direction.  That way, we can pack them
# into a 4-proc MPI job that happens to fit PJ's workstation nicely :) 
nx0 = 5; nx1 = 30; ny = 30; nz = 30; nbk = 3
clustery = RobertsClusterFunction(1, 0, 1.2)
cf_list=[None, clustery, None, clustery,
         None, clustery, None, clustery,
         None, None, None, None]
upstr = SuperBlock3D(PyFunctionVolume(revolved_vol0), nbk=nbk,
                     nni=nx0, nnj=ny, nnk=nz, cf_list=cf_list,
                     fill_condition=initial, label="upstr")
cone = SuperBlock3D(PyFunctionVolume(revolved_vol1), nbk=nbk,
                    nni=nx1, nnj=ny, nnk=nz, cf_list=cf_list,
                    fill_condition=initial, label="cone")
# The coincidental vertices confuse the automated search for
# block connections so we need to make the connections manually.
# identify_block_connections()
for k in range(nbk):
    connect_blocks_3D(upstr.blks[0][0][k], cone.blks[0][0][k],
                      ((2,3),(6,7),(1,0),(5,4))) # EAST --> WEST
    upstr.blks[0][0][k].bc_list[WEST] = SupInBC(inflow)
    cone.blks[0][0][k].bc_list[EAST] = ExtrapolateOutBC()

# Do a little more setting of global data.
gdata.viscous_flag = 0
gdata.flux_calc = ADAPTIVE
gdata.max_time = 5.0e-3  # seconds
gdata.max_step = 10000
gdata.dt = 1.0e-7
gdata.dt_plot = 1.5e-3
gdata.dt_history = 10.0e-5
