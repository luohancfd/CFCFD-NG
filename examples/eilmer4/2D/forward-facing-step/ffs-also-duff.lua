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
speed_of_sound = initial:toTable().a
print("speed_of_sound=", speed_of_sound)
inflow = FlowState:new{p=101.325e3, T=300.0, velx=3.0*speed_of_sound, vely=0.0}

a0 = Vector3:new{0.0,0.0}; a1 = Vector3:new{0.0,0.2}; a2 = Vector3:new{0.0,1.0}
b0 = Vector3:new{0.6,0.0}; b1 = Vector3:new{0.6,0.2}; b2 = Vector3:new{0.6,1.0}
c1 = Vector3:new{3.0,0.2}; c2 = Vector3:new{3.0,1.0}
surf0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
surf1 = CoonsPatch:new{p00=a1, p10=b1, p11=b2, p01=a2}
surf2 = CoonsPatch:new{p00=b1, p10=c1, p11=c2, p01=b2}

-- Mesh the patches, with particular discretisation.
dx = 10.0e-3
nab = math.floor(0.6/dx); nbc = math.floor(3.4/dx)
print("nab=", nab, "nbc=", nbc)
n01 = math.floor(0.2/dx); n12 = math.floor(0.8/dx)
print("n01=", n01, "n12=", n12)
print("hello 1")

-- Define the flow-solution blocks.
blk0 = SBlock:new{grid=StructuredGrid2D:new{surf=CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1},
					    niv=nab+1, njv=n01+1},
		  fillCondition=FlowState:new{p=101.325e3, T=300.0, velx=3.0*speed_of_sound, vely=0.0},
		  label="BLOCK-0"}
print("hello 2")
blk1 = SBlock:new{grid=StructuredGrid2D:new{surf=CoonsPatch:new{p00=a1, p10=b1, p11=b2, p01=a2},
					    niv=nab+1, njv=n12+1},
		  fillCondition=FlowState:new{p=101.325e3, T=300.0, velx=3.0*speed_of_sound, vely=0.0},
		  label="BLOCK-1"}
blk2 = SBlock:new{grid=StructuredGrid2D:new{surf=CoonsPatch:new{p00=b1, p10=c1, p11=c2, p01=b2},
					    niv=nbc+1, njv=n12+1},
		  fillCondition=FlowState:new{p=101.325e3, T=300.0, velx=3.0*speed_of_sound, vely=0.0},
		  label="BLOCK-2"}
print("hello 3")

-- Set boundary conditions.
identifyBlockConnections()
print("hello 4")
blk0.bcList[west] = SupInBC:new{flowCondition=FlowState:new{p=101.325e3, T=300.0, velx=3.0*speed_of_sound, vely=0.0}, label="inflow-boundary"}
blk1.bcList[west] = SupInBC:new{flowCondition=FlowState:new{p=101.325e3, T=300.0, velx=3.0*speed_of_sound, vely=0.0}, label="inflow-boundary"}
blk2.bcList[east] = ExtrapolateOutBC:new{label="outflow-boundary"}

-- Do a little more setting of global data.
config.max_time = 5.0e-3  -- seconds
config.max_step = 3000
config.dt_init = 1.0e-6
config.cfl_value = 0.5
-- config.dt_max = 10.0e-6
config.dt_plot = 1.0e-3
config.dt_history = 10.0e-6
print("hello 5")