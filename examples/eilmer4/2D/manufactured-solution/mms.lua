-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-03-17
-- History: Ported from eilmer3 example
--          Presently only Euler case on regular grid

config.title = "Method of Manufactured Solutions, Euler case."
print(config.title)
config.dimensions = 2

setGasModel('ideal-air-gas-model.lua')
p0 = 1.0e5; T0 = p0/287.10325; u0 = 800.0; v0 = 800.0
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0}

p00 = Vector3:new{0.0, 0.0}
p10 = Vector3:new{1.0, 0.0}
p01 = Vector3:new{0.0, 1.0}
p11 = Vector3:new{1.0, 1.0}
nicell = 16; njcell = 16
grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01},
			  niv=nicell+1, njv=njcell+1}

bcList = {}
bcList[north] = ExtrapolateOutBC:new{xOrder=1}
bcList[east] = ExtrapolateOutBC:new{xOrder=1}
bcList[south] = UserDefinedBC:new{fileName='udf-bc.lua'}
bcList[west] = UserDefinedBC:new{fileName='udf-bc.lua'}
blks = SBlockArray{grid=grid, fillCondition=initial, bcList=bcList, 
		   nib=2, njb=2, label="blk"}

config.interpolation_order = 2
config.gasdynamic_update_scheme = "predictor-corrector"
config.flux_calculator = 'ausmdv'
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.dt_init = 1.0e-6
config.max_time = 60.0e-3
config.dt_plot = config.max_time/20.0
config.max_step = 1000000
config.cfl = 0.5
config.stringent_cfl = 1
config.apply_limiter = false
config.extrema_clipping = false

				 
		
