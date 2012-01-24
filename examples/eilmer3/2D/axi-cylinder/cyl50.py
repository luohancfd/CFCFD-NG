## \file cyl50.py
## \author PJ
## \version 14-Aug-2006 updated from Tcl script
## \version 18-Jan-2010 updated for Eilmer3

gdata.title = "Mach 2 flow along the axis of a 5mm cylinder."
print gdata.title
gdata.axisymmetric_flag = 1

# Accept defaults giving perfect air (R=287 J/kg.K, gamma=1.4)
select_gas_model(model='ideal gas', species=['air'])

inflow  = FlowCondition(p=257.3, u=597.3, v=0.0, T=222.0)

# Set up a quadrilateral in the (x,y)-plane.
#  y     c
#  ^   / |
#  | /   |
#  d     |
#  a-----b
#  0------------> x
a = Node(0.0,0.005); b = Node(1.0,0.005); c = Node(1.0,0.7); d = Node(0.0,0.06)
south = Line(a,b); north = Line(d,c); west = Line(a,d); east = Line(b,c)

# The following lists are in order [N, E, S, W]
bndry_list = [SupInBC(inflow), ExtrapolateOutBC(), FixedTBC(222.0), SupInBC(inflow)]
cluster_list = [RobertsClusterFunction(1,0,1.1), RobertsClusterFunction(1,0,1.01),
                RobertsClusterFunction(1,0,1.1), RobertsClusterFunction(1,0,1.01)]

# Assemble the block from the geometry, discretization and boundary data.
blk = Block2D(psurf=make_patch(north, east, south, west, grid_type="AO"),
              nni=50, nnj=50,
              bc_list=bndry_list, cf_list=cluster_list,
              fill_condition=inflow)

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 8.0e-3  # seconds
gdata.max_step = 230000
gdata.dt = 3.0e-8
gdata.dt_plot = 4.0e-3

sketch.xaxis(0.0, 1.0, 0.2, -0.05)
sketch.yaxis(0.0, 1.0, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)
