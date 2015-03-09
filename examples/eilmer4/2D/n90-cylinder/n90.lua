-- n90.lua -- Cylinder in dissociating nitrogen flow
-- Rg & PJ  2015-03-09

job_title = "Cylinder in dissociating nitrogen flow."
print(job_title)

config.dimensions = 2
config.title = job_title

nsp, nmodes = setGasModel('nitrogen-2sp.lua')
print("setting inflow")
inflow = FlowState:new{p=500.0, T=700.0, u=5000.0, massf={1.0, 0.0}}
print "setting initial"
initial = FlowState:new{p=5.0, T=300.0, massf={1.0, 0.0}}


print "Reading Gridpro grid."


-- Grid has been built earlier in GridPro
gproName = 'n90.gpro'
grid = importGridproGrid(gproName)

blk0 = SBlock:new{grid=grid[1], fillCondition=inflow, label="blk-0"}
-- We can leave east and south as SlipWalls
blk0.bcList[west] = SupInBC:new{flowCondition=inflow}
blk0.bcList[north] = ExtrapolateOutBC:new{}

-- Set a few more config options
config.flux_calc = ADAPTIVE
config.max_time = 100.0e-6
config.max_step = 40000
config.dt_init = 1.0e-8
config.cfl = 0.5
config.dt_plot = 20.0e-6
