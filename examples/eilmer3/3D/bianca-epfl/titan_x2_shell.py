# titan_x2_shell.py
# A job description file for Bianca's Titan Aeroshell used in X2.
# PJ
# 30-Oct-2006: Elmer2 original
# 07-Feb-2010: Eilmer3 port


# First, set the global data
gdata.title = "Titan Aeroshell used in X2."
gdata.dimensions = 3
# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
gdata.max_time = 300.0e-3
gdata.dt = 1.0e-7
gdata.max_step = 5000
gdata.dt_plot = 30.0e-3

# Second, set up flow conditions
from math import pi, sin, cos
alpha = 20.0*pi/180.0 # angle of attack in radians
initialCond = FlowCondition(p=1000.0, u=0.0, T=300.0)
M_inf = 7.0
u_inf = M_inf * initialCond.flow.gas.a
inflowCond  = FlowCondition(p=50.0e3, u=-u_inf*cos(alpha), 
                            v=u_inf*sin(alpha), T=300.0)

# Third, set up the blocks from the ICEM-generated grids.
# The discretization is just a fraction of the original ICEM grids.
# block      0  1  2  3  4  5  6  7  8  9 10 11 12
nni_list = [10,10,10,24,10,10,10,10,10, 7,10,10, 3]
nnj_list = [10,24,24,10,10,10,10, 7, 7,10, 3, 3,10]
nnk_list = [30,30,30,30,30,30,30,30,30,30,30,30,30]
pv_list = []
blk_list = []
for ib in range(13):
    pv_list.append( MeshVolume("icem_grid."+str(ib)+".g.vtk") )
    blk_list.append( Block3D(nni=nni_list[ib], 
                             nnj=nnj_list[ib], 
                             nnk=nnk_list[ib],
                             parametric_volume=pv_list[ib], 
                             fill_condition=initialCond) )
identify_block_connections()

# Apply boundary conditions.
# The appropriate surfaces were determined by loading each block
# with MayaVi, then putting on a gridplane, and fiddling with the
# index directions to find out which surface was which.
blk_list[0].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[0].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[1].set_BC("BOTTOM", "SUP_IN", inflow_condition=inflowCond)
blk_list[1].set_BC("TOP", "FIXED_T", Twall=300.0)
blk_list[1].set_BC("SOUTH", "SUP_OUT")
blk_list[2].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[2].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[2].set_BC("SOUTH", "SUP_OUT")
blk_list[3].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[3].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[3].set_BC("WEST", "SUP_OUT")
blk_list[4].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[4].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[5].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[5].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[6].set_BC("BOTTOM", "SUP_IN", inflow_condition=inflowCond)
blk_list[6].set_BC("TOP", "FIXED_T", Twall=300.0)
blk_list[7].set_BC("BOTTOM", "SUP_IN", inflow_condition=inflowCond)
blk_list[7].set_BC("TOP", "FIXED_T", Twall=300.0)
blk_list[8].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[8].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[9].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[9].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[10].set_BC("BOTTOM", "SUP_IN", inflow_condition=inflowCond)
blk_list[10].set_BC("TOP", "FIXED_T", Twall=300.0)
blk_list[11].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[11].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)
blk_list[12].set_BC("BOTTOM", "FIXED_T", Twall=300.0)
blk_list[12].set_BC("TOP", "SUP_IN", inflow_condition=inflowCond)

