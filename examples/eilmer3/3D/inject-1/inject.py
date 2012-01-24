# inject.py -- single discrete-hole injection.
# PJ
# Elmer2 original: Nov-2006
# Eilmer3 port: 06-Feb-2010

# -------------------------------------------------------------------
# Some handy definitions for later.
import math
# from cfpylib.geom.box3d import makeSimpleBox

# ---------------- First, set the global data ----------------------
gdata.title = "Single-hole injection."
gdata.dimensions = 3
gdata.dt = 1.0e-8
gdata.t_order = 1
gdata.max_time = 3.0e-3
gdata.max_step = 60000
gdata.reacting_flag = 0
gdata.dt_plot = 0.5e-3
gdata.dt_history = 1.0e-5

# ------------ Second, set up flow conditions -------------------
# These will be used for fill and boundary conditions.
species_list = select_gas_model(model='ideal gas', species=['N2', 'H2'])
initialCond = FlowCondition(p=5.955e3, u=0.0, T=304.0, massf={'N2':1.0})
inflowCond = FlowCondition(p=95.84e3, u=1000.0, T=1103.0, massf={'N2':1.0})
injectCond = FlowCondition(p=95.84e3, w=1000.0, T=300.0, massf={'H2':1.0})
 
# ------------ Third, set up the blocks ---------------------
# Parameters defining the duct...
L0 = 20.0e-2    # length of duct in flow direction
L1 = 5.0e-2     # distance from leading edge to injector
L2 = 1.0e-2     # streamwise length of injector
Whalf0 = 5.0e-2 # half-width of duct
Whalf1 = 5.0e-3 # half-width of injector
H = 5.0e-2      # height of duct

# Plan of blocks
#                 NORTH BNDRY
#         +--------+---+---------------+
#         |        |   |               |
#         |   01   | 11|      21       |
# inflow> |        |   |               | outflow>
# (WEST)  +--------+---+---------------+ (EAST)
#         |   00   | 10|      20       |
#         +--------+---+---------------+
#                 SOUTH BNDRY
#
#                    ^
#                 injector
cluster_k  = RobertsClusterFunction(1, 0, 1.2) # cluster down, toward the bottom surface
cluster_i0 = RobertsClusterFunction(0, 1, 1.2) # cluster streamwise toward injector
cluster_i2 = RobertsClusterFunction(1, 0, 1.2)
cluster_j1 = RobertsClusterFunction(1, 0, 1.2) # cluster cross-stream toward injector
# upstream pair of blocks
pv = makeSimpleBox(xPos=0.0, yPos=0.0, xSize=L1, ySize=Whalf1, zSize=H)
cflist = [cluster_i0,None,cluster_i0,None]*2 + [cluster_k,]*4;
# 12 edges is a full complement; see elmer_prep.py for the order of edges
blk00 = Block3D(nni=40, nnj=10, nnk=30, parametric_volume=pv,
                cf_list=cflist, fill_condition=initialCond)
pv = makeSimpleBox(xPos=0.0,yPos=Whalf1, xSize=L1,ySize=Whalf0-Whalf1, zSize=H)
cflist = [cluster_i0,cluster_j1,cluster_i0,cluster_j1]*2 + [cluster_k,]*4;
blk01 = Block3D(nni=40, nnj=30, nnk=30, parametric_volume=pv,
                cf_list=cflist, fill_condition=initialCond)

# injector and part of plate beside it
pv = makeSimpleBox(xPos=L1,yPos=0.0, xSize=L2, ySize=Whalf1, zSize=H)
cflist = [None,None,None,None]*2 + [cluster_k,]*4;
blk10 = Block3D(nni=10, nnj=10, nnk=30, parametric_volume=pv,
                cf_list=cflist, fill_condition=initialCond)
pv = makeSimpleBox(xPos=L1,yPos=Whalf1, xSize=L2, ySize=Whalf0-Whalf1, zSize=H)
cflist = [None,cluster_j1,None,cluster_j1]*2 + [cluster_k,]*4;
blk11 = Block3D(nni=10, nnj=30, nnk=30, parametric_volume=pv,
                cf_list=cflist, fill_condition=initialCond)

# blocks downstream of injector
pv = makeSimpleBox(xPos=L1+L2,yPos=0.0, xSize=L0-(L1+L2), ySize=Whalf1, zSize=H)
cflist = [cluster_i2,None,cluster_i2,None]*2 + [cluster_k,]*4;
blk20 = Block3D(nni=50, nnj=10, nnk=30, parametric_volume=pv,
                cf_list=cflist, fill_condition=initialCond)
pv = makeSimpleBox(xPos=L1+L2,yPos=Whalf1, xSize=L0-(L1+L2), ySize=Whalf0-Whalf1, zSize=H)
cflist = [cluster_i2,cluster_j1,cluster_i2,cluster_j1]*2 + [cluster_k,]*4;
blk21 = Block3D(nni=50, nnj=30, nnk=30, parametric_volume=pv,
                cf_list=cflist, fill_condition=initialCond)


identify_block_connections()
blk00.set_BC("WEST", "SUP_IN", inflow_condition=inflowCond)
blk01.set_BC("WEST", "SUP_IN", inflow_condition=inflowCond)
blk10.set_BC("BOTTOM", "SUP_IN", inflow_condition=injectCond)
blk20.set_BC("EAST", "SUP_OUT")
blk21.set_BC("EAST", "SUP_OUT")
