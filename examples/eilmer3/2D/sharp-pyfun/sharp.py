# sharp-pyfun/sharp.py
# PJ, 14-Dec-2006
#     16-Sep-2008 ported to Elmer3
#     29-Apr-2009 PyPath used instead of spline.
#
gdata.title = "Mach 3.0 flow over a curved 2D-planar body."
print gdata.title

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
# Define flow conditions.
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)
inflow  = FlowCondition(p=95.84e3, u=2000.0, v=0.0, T=1103.0)
# One can get access to the details of the FlowCondition.
print "inflow M=", inflow.flow.u / inflow.flow.gas.a

# Geometry of flow domain.
def y(x):
    "(x,y)-space path for x>=0"
    if x <= 3.291:
        return -0.008333 + 0.609425*x - 0.092593*x*x
    else:
        return 1.0

def xypath(t):
    "Parametric path with 0<=t<=1."
    global y
    x = 10.0 * t
    yval = y(x)
    if yval < 0.0:
        yval = 0.0
    return (x, yval, 0.0)

a = Node(-1.0, 0.0, label="A")
b = Node( 0.0, 0.0, label="B")
c = Node(10.0, 1.0, label="C")
d = Node(10.0, 7.0, label="D")
e = Node( 0.0, 7.0, label="E")
f = Node(-1.0, 7.0, label="F")

north0 = Line(f, e)
e0w1   = Line(b, e)
south0 = Line(a, b)
west0  = Line(a, f)
south1 = PyFunctionPath(xypath)
north1 = Line(e, d)
east1  = Line(c, d)

# Define the blocks, grid resolution and boundary conditions.
ny = 60
clustery = RobertsClusterFunction(1, 0, 1.3)
clusterx = RobertsClusterFunction(1, 0, 1.2)
blk_0 = Block2D(make_patch(north0, e0w1, south0, west0),
                nni=16, nnj=ny,
                cf_list=[None,clustery,None,clustery],
                fill_condition=initial)
blk_1 = Block2D(make_patch(north1, east1, south1, e0w1),
                nni=80, nnj=ny,
                cf_list=[clusterx,None,clusterx,clustery],
                fill_condition=initial)
identify_block_connections()
blk_0.bc_list[WEST]=SupInBC(inflow)
blk_1.bc_list[EAST]=ExtrapolateOutBC()

# Do a little more setting of global data.
gdata.flux_calc = ADAPTIVE
gdata.max_time = 15.0e-3  # seconds
gdata.max_step = 2500
gdata.dt = 1.0e-6

sketch.xaxis(0.0,10.0, 2.0, -0.6)
sketch.yaxis(0.0, 8.0, 2.0, -1.6)
sketch.window(0.0, 0.0, 10.0, 10.0, 0.05, 0.05, 0.17, 0.17)
