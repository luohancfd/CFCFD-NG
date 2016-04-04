# ffs.py -- Forward-facing-step example
# PJ & RG  2016-04-04 back-ported from the Eilmer4 example

# We can set individual attributes of the global data object.
gdata.title = "Forward-facing step with supersonic flow."
print(gdata.title)
gdata.dimensions = 2

# Gas model and flow conditions.
select_gas_model(model='ideal gas', species=['air'])
initial = FlowCondition(p=101.325e3, T=300.0)
sspeed = initial.to_dict()['a']
print("sound speed=", sspeed)
inflow = FlowCondition(p=101.325e3, T=300.0, u=3.0*sspeed, v=0.0)

# Geometry of the flow domain.
a0 = Vector3(0.0, 0.0)
a1 = Vector3(0.0, 0.2)
a2 = Vector3(0.0, 1.0)
b0 = Vector3(0.6, 0.0)
b1 = Vector3(0.6, 0.2)
b2 = Vector3(0.6, 1.0)
c1 = Vector3(3.0, 0.2)
c2 = Vector3(3.0, 1.0)
surf0 = CoonsPatch(a0, b0, b1, a1)
surf1 = CoonsPatch(a1, b1, b2, a2)
surf2 = CoonsPatch(b1, c1, c2, b2)

# Mesh the patches, with particular discretisation.
dx = 10.0e-3
nab = int(0.6/dx); nbc = int(2.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = int(0.2/dx); n12 = int(0.8/dx)
print("n01=", n01, "n12=", n12)

# Set boundary conditions that we care about.
bcList0 = [SlipWallBC(), SlipWallBC(), SlipWallBC(), SupInBC(inflow)]
bcList1 = [SlipWallBC(), SlipWallBC(), SlipWallBC(), SupInBC(inflow)]
bcList2 = [SlipWallBC(), ExtrapolateOutBC(), SlipWallBC(), SlipWallBC()]

# Define the flow-solution blocks and stitch them together.
blk0 = SuperBlock2D(surf0, nni=nab, nnj=n01, nbi=1, nbj=1, 
                    fill_condition=inflow, bc_list=bcList0, label="BLOCK-0")
blk1 = SuperBlock2D(surf1, nni=nab, nnj=n12, nbi=1, nbj=4,
                    fill_condition=inflow, bc_list=bcList1, label="BLOCK-1")
blk2 = SuperBlock2D(surf2, nni=nbc, nnj=n12, nbi=4, nbj=4,
		   fill_condition=inflow, bc_list=bcList2, label="BLOCK-2")
identify_block_connections()

# Do a little more setting of global data.
gdata.max_time = 5.0e-3  # seconds
gdata.max_step = 6000
gdata.dt = 1.0e-6
gdata.cfl = 0.5
gdata.dt_plot = 1.0e-3
gdata.dt_history = 10.0e-6
