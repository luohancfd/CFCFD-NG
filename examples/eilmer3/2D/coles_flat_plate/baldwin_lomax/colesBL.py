## \file colesBL.py
## \brief Turbulent flow over a flat plate 
##        -- Coles test case
##        -- Baldwin-Lomax turbulence model
## \author PJ (30-Sep-2008) and Jan-Pieter Nap 

job_title = "Coles Mach 3.7 flow over a flat plate (Baldwin-Lomax)"
print job_title
gdata.dimensions = 2
# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define flow conditions to match Coles' data set 53010801
p_inf = 1.358e3  # Pa
u_inf = 677.4    # m/s
T_inf = 83.34    # degrees K

inflow = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf)
if 1:
    # It seems that everyone has been starting the flow field with
    # inflow conditions throughout and then letting the boundary layer
    # grow out into the main stream.
    initial = inflow
else:
    # Here is a low-pressure initial state more like a shock tunnel.
    # Unfortunately, it seems to play havoc with the turbulence.
    initial = FlowCondition(p= 0.1*p_inf, u=0.0, v=0.0, T=296.0)


# Geometry of plate and flow domain.
L = 0.60 # metres
H = 0.4 * L
NB = 4 	 # number of blocks

#         wall
#        c---------b
# flow=> |         |
#        d         |
#          -\-     |
#    flow=>    -\- |
#        0         a ----> x
# 
a = Node(L, 0.0); b = Node(L, H); c = Node(0.0, H); d = Node(0.0, 3.0*H/4.0)
north = Line(c,b); east = Line(a,b); south = Line(d,a); west = Line(d,c)

# Define the blocks, boundary conditions and set the discretisation.
blk = SuperBlock2D(make_patch(north, east, south, west), 
                   nni=100, nnj=90, nbi=NB, nbj=1, 
                   fill_condition=initial,
                   cf_list=[RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.01),
                            RobertsClusterFunction(1,0, 1.05),
                            RobertsClusterFunction(0,1, 1.05)],
                   bc_list=[AdiabaticBC(), ExtrapolateOutBC(),
                            SupInBC(inflow), SupInBC(inflow)])


# Do a little more setting of global data.
gdata.turbulence_flag = 1 # to activate our turbulence model
gdata.turbulence_model = "baldwin_lomax" # not the default k_omega
gdata.title = job_title
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
 	
gdata.max_time = 5.0e-3  # should allow a few flow lengths   
gdata.dt_plot =  1.0e-3
gdata.dt_history = 1.0e-5
gdata.max_step = 3000000 

gdata.cfl = 0.4	
gdata.cfl_count = 3
gdata.stringent_cfl = 0 # 1 is more robust
gdata.dt = 1.0e-9	# only an initial guess, the simulation will take this over

# The following scales provide a reasonable picture.
sketch.xaxis(0.0, 0.6, 0.1, -0.05)
sketch.yaxis(0.0, 0.5, 0.1, -0.04)
sketch.window(0.0, 0.0, 0.6, 0.6, 0.05, 0.05, 0.17, 0.17)

