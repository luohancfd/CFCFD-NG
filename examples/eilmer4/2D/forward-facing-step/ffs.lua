-- ffs.lua -- Forward-facing-step example
-- PJ & RG  2015-03-08

job_title = "Forward-facing step with supersonic flow."
print(job_title)

-- We can set individual attributes of the global data object.
config.dimensions = 2
config.title = job_title

nsp, nmodes = setGasModel('ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=101.325e3, T=300.0}
sspeed = initial:toTable().a
print("sound speed=", sspeed)
inflow = FlowState:new{p=101.325e3, T=300.0, velx=3.0*sspeed, vely=0.0}

a0 = Vector3:new{0.0,0.0}; a1 = Vector3:new{0.0,0.2}; a2 = Vector3:new{0.0,1.0}
b0 = Vector3:new{0.6,0.0}; b1 = Vector3:new{0.6,0.2}; b2 = Vector3:new{0.6,1.0}
c1 = Vector3:new{3.0,0.2}; c2 = Vector3:new{3.0,1.0}
surf0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
surf1 = CoonsPatch:new{p00=a1, p10=b1, p11=b2, p01=a2}
surf2 = CoonsPatch:new{p00=b1, p10=c1, p11=c2, p01=b2}

-- Mesh the patches, with particular discretisation.
dx = 10.0e-3
nab = math.floor(0.6/dx); nbc = math.floor(2.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = math.floor(0.2/dx); n12 = math.floor(0.8/dx)
print("n01=", n01, "n12=", n12)
grid0 = StructuredGrid2D:new{surf=surf0, niv=nab+1, njv=n01+1}
grid1 = StructuredGrid2D:new{surf=surf1, niv=nab+1, njv=n12+1}
grid2 = StructuredGrid2D:new{surf=surf2, niv=nbc+1, njv=n12+1}

-- Define the flow-solution blocks.
blk0 = SBlock:new{grid=grid0, fillCondition=inflow, label="BLOCK-0"}
blk1 = SBlock:new{grid=grid1, fillCondition=inflow, label="BLOCK-1"}
blk2 = SBlock:new{grid=grid2, fillCondition=inflow, label="BLOCK-2"}

-- Set boundary conditions.
identifyBlockConnections()
blk0.bcList[west] = SupInBC:new{flowCondition=inflow, label="inflow-boundary"}
blk1.bcList[west] = SupInBC:new{flowCondition=inflow, label="inflow-boundary"}
blk2.bcList[east] = ExtrapolateOutBC:new{label="outflow-boundary"}

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 30
config.dt_init = 1.0e-6
config.cfl_value = 0.5
-- config.dt_max = 10.0e-6
config.dt_plot = 1.0e-3
config.dt_history = 10.0e-6
