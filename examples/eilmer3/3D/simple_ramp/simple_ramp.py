# A sample job description file ... is actually Python code.
# This is a fudged version of the cone20 case from mb_cns in 3D.
# It is now a ramp at 10 degrees rather than a conical surface.
# PJ, August 2004, Jan 2006, Jul 2006 (new thermochemistry module)
#     July 2008 Eilmer3 port by adding gdata.dimensions=3
# -------------------------------------------------------------------

# ---------------- First, set the global data ----------------------
# To see what parameters one can set, look up the class definition
# in the file e3prep.py.

gdata.title = "Ramp at 20 degrees."
gdata.dimensions = 3

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
gdata.viscous_flag = 0
gdata.max_time = 5.0e-3
gdata.max_step = 1000
gdata.reacting_flag = 0
# Set some of the other properties separately, just for fun.
gdata.t_order = 1
gdata.x_order = 2
# gdata.stringent_cfl = 1
gdata.dt_plot = 1.0e-3
gdata.dt_history = 1.0e-5

# ------------ Second, set up flow conditions -------------------
# These will be used for fill and boundary conditions.
initialCond = FlowCondition(p=5.955e3, u=0.0, T=304.0, massf=[1.0,])
inflowCond = FlowCondition(p=95.84e3, u=1000.0, T=1103.0, massf=[1.0,])

# ------------ Third, set up the blocks ---------------------
# These may explicitly reference previously defined flow conditions
# but, even if they don't, their setup implicitly references the
# first flow condition.

# Note that we can use the Python language to do some of our
# calculations.  Here are some handy definitions for later.

def toRadians(degrees):
    import math
    return degrees * math.pi / 180.0

def simpleBoxCorners(xPos=0.0, yPos=0.0, zPos=0.0, xSize=1.0, ySize=1.0, zSize=1.0):
    """\brief Creates a corner coordinate list for a simple box."""
    p0 = Node(xPos,       yPos,       zPos)
    p1 = Node(xPos+xSize, yPos,       zPos)
    p2 = Node(xPos+xSize, yPos+ySize, zPos)
    p3 = Node(xPos,       yPos+ySize, zPos)
    p4 = Node(xPos,       yPos,       zPos+zSize)
    p5 = Node(xPos+xSize, yPos,       zPos+zSize)
    p6 = Node(xPos+xSize, yPos+ySize, zPos+zSize)
    p7 = Node(xPos,       yPos+ySize, zPos+zSize)
    return [p0, p1, p2, p3, p4, p5, p6, p7]

def makeSimpleBox(p):
    return SimpleBoxVolume(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7])

# -------------------------------------------------------------------
# First block is the region in front of the ramp. 10x40(x4)
pvolume = makeSimpleBox(simpleBoxCorners(xSize=0.2,ySize=0.1))
cluster_k = RobertsClusterFunction(1, 0, 1.2) # cluster down, toward the wedge surface
cflist = [None,]*8 + [cluster_k,]*4; # 12 edges is a full complement
blk0 = Block3D(label="first-block", nni=10, nnk=40,
               parametric_volume=pvolume,
               cf_list=cflist,
               fill_condition=initialCond)
blk0.set_BC("WEST", "SUP_IN", inflow_condition=inflowCond)

# For the grid over the ramp, start with a regular box... 30x40(x4)
blk1Corners = simpleBoxCorners(xPos=0.2,xSize=0.8,ySize=0.1)
# Now, raise the end of the ramp.
blk1Corners[1].z = 0.8 * math.tan(toRadians(10.0))
blk1Corners[2].z = blk1Corners[1].z
blk1 = Block3D(label="second-block", nni=30, nnk=40,
               parametric_volume=makeSimpleBox(blk1Corners),
               cf_list=cflist,
               fill_condition=initialCond,
               hcell_list=[(1,1,2),(20,1,1)])
blk1.set_BC("EAST", "SUP_OUT")
identify_block_connections()
