# sharp.py
# PJ, 14-Dec-2006
#     16-Sep-2008 ported to Eilmer3
job_title = "Mach 3.0 flow over a curved 2D-planar body."
print job_title
gdata.title = job_title

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
# Define flow conditions
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)
inflow  = FlowCondition(p=95.84e3, u=2000.0, v=0.0, T=1103.0)

# Geometry
def shape(x):
    return -0.008333 + 0.609425*x - 0.092593*x*x

a = Node(-1.0, 0.0, label="A")
b = Node( 0.0, 0.0, label="B")
x_list = [0.5, 1.5, 2.5, 3.291]
b_list = [b,] # to accumulate points in the spline
for x in x_list:
    b_list.append(Node(x, shape(x)))
c = Node(10.0, b_list[-1].y, label="C") # extend at same y-value
d = Node(10.0, 7.0, label="D")
e = Node( 0.0, 7.0, label="E")
f = Node(-1.0, 7.0, label="F")

north0 = Line(f, e)
e0w1 = Line(b, e)
south0 = Line(a, b)
west0 = Line(a, f)
south1 = Polyline([Spline(b_list), Line(b_list[-1], c)])
north1 = Line(e, d)
east1 = Line(c, d)

# Define the blocks, boundary conditions and set the discretisation.
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
