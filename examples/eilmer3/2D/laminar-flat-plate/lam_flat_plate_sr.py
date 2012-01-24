## \file lam_flat_plate_sr.py
## \brief Laminar flow over a flat plate, singularity removed. 
##
##  This script creates and runs a simulation of a worked
##  example (Example 5-11 on page 155) in Joseph Schetz's
##  book "Boundary Layer Analysis".
##
##  The example involves a laminar Mach 4.0 flow over a
##  1.0 m long flat plate. The boundary layer will grow on
##  the NORTH boundary but the first 10cm of the plate has
##  been removed and the inflow boundary now includes the
##  profile for a self-similar boundary layer.
##
##  Again, Simulation results will be compared 
##  with those from the CLBL code from Schetz's book.
##
## \author PJ 07-mar-2011 adapted from lam_flat_plate.py

gdata.title = "Schetz's Mach 4 laminar flow over a flat plate, singularity removed"
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

select_gas_model(fname='air_pr_0p71.lua')

# Define flow conditions to match Schetz's worked example 5-1
p_inf = 1.013e3  # Pa
u_inf = 1390.0   # m/s
T_inf = 300.0    # degrees K

external_flow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=[1.0,])
# Read lists defining the boundary-layer profile so that we can use the profile
# as the starting flow throughout the domain.
execfile('profile.py')
# When defining the flow function, we pass through the recently read lists.
def fill_gas(x, y, z, extflow=external_flow, ylist=y, plist=p, Tlist=T, ulist=u):
    if y < ylist[-1]:
        return extflow.to_dict()
    else:
        nearest = 0
        d_min = abs(y - ylist[0])
        for i in range(1,len(ylist)):
            d = abs(y - ylist[i])
            if d < d_min:
                d_min = d
                nearest = i
        flow = FlowCondition(p=plist[nearest], u=ulist[nearest], v=0.0,
                             T=Tlist[nearest], massf=[1.0,], add_to_list=0)
        return flow.to_dict()


# Geometry of plate and flow domain.
# Have set up the flat plate 0.1 m longer than the actual plate
# so that we don't have to use the simulation profiles right on 
# the boundary.
L = 1.1       # Length of original flat plate (in metres)
H = 0.4 * L   # Height of flow domain (in metres)

#         wall
#        c---------b
# flow=> |         |
#        d         |
#          -\-     |
#    flow=>    -\- |
#        0         a ----> x
#
# Relative to the original geometry, we will truncate 1/11 off the
# upstream end of the flow domain but try to keep the rest
# of the boundaries the same.
a = Node(L, 0.0, label='a'); b = Node(L, H, label='b')
c = Node(0.1, H, label='c'); d = Node(0.1, 3.0*H/4.0*10.0/11, label='d')
north = Line(c,b); east = Line(a,b); south = Line(d,a); west = Line(d,c)

# Split the domain into a grid of blocks so that we can make use
# of our cluster computer.
bc_list=[FixedTBC(300.0), ExtrapolateOutBC(), SupInBC(external_flow),
         UserDefinedBC(filename='udf-boundary-layer-in.lua')]
blk = SuperBlock2D(make_patch(north, east, south, west), 
                   nni=220, nnj=192, nbi=4, nbj=4, 
                   fill_condition=fill_gas,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.016),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.05)],
                   bc_list=bc_list)
# Set to make a nice 2D rendering of the blocks.
sketch.xaxis(0.0, 1.2, 0.2, -0.05)
sketch.yaxis(0.0, 0.6, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)

