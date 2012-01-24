## \file lam_flat_plate.py
## \brief Laminar flow over a flat plate 
##
##  This script creates and runs a simulation of a worked
##  example (Example 5-11 on page 155) in Joseph Schetz's
##  book "Boundary Layer Analysis".
##
##  The example involves a laminar Mach 4.0 flow over a
##  1.0 m long flat plate. The boundary layer will grow on
##  the NORTH boundary. Simulation results will be compared 
##  with those from the CLBL code from Schetz's book.
##
## \author Peter Jacobs & Wilson Chan, 18 Jan 2010

gdata.title = "Schetz's Mach 4 laminar flow over a flat plate"
print gdata.title
gdata.viscous_flag = 1
gdata.flux_calc = AUSMDV
gdata.max_time = 2.4e-3  # will allow 3 flow lengths   
gdata.dt_plot =  0.4e-3
gdata.dt_history = 1.0e-5
gdata.max_step = 300000
gdata.cfl = 0.4
gdata.cfl_count = 3
gdata.dt = 1.0e-9  # only an initial guess - the simulation will take over

select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Schetz's worked example 5-1
p_inf = 1.013e3  # Pa
u_inf = 1390.0   # m/s
T_inf = 300.0    # degrees K

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,])

# Geometry of plate and flow domain.
# Have set up the flat plate 0.1 m longer than the actual plate
# so that we don't have to use the simulation profiles right on 
# the boundary.
L = 1.1       # Length of flat plate (in metres)
H = 0.4 * L   # Height of flow domain (in metres)

#         wall
#        c---------b
# flow=> |         |
#        d         |
#          -\-     |
#    flow=>    -\- |
#        0         a ----> x
# 
a = Node(L, 0.0, label='a'); b = Node(L, H, label='b')
c = Node(0.0, H, label='c'); d = Node(0.0, 3.0*H/4.0, label='d')
north = Line(c,b); east = Line(a,b); south = Line(d,a); west = Line(d,c)

# Split the domain into a grid of blocks so that we can make use
# of our cluster computer.
blk = SuperBlock2D(make_patch(north, east, south, west), 
                   nni=220, nnj=192, nbi=4, nbj=4, 
                   fill_condition=inflow,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.016),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.05)],
                   bc_list=[FixedTBC(300.0), ExtrapolateOutBC(),
                            SupInBC(inflow), SupInBC(inflow)])
# Set to make a nice 2D rendering of the blocks.
sketch.xaxis(0.0, 1.2, 0.2, -0.05)
sketch.yaxis(0.0, 0.6, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)

