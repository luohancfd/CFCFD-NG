-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-03-17
-- History: Ported from eilmer3 example
--          Presently only Euler case on regular grid

title = "Method of Manufactured Solutions, Euler case."
print(title)
config.dimensions = 2
config.title = title

setGasModel('ideal-air-gas-model.lua')

p0 = 1.0e5
T0 = p0/287.10325
u0 = 800.0
v0 = 800.0

initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0}

a = Vector3:new{0.0, 0.0}
b = Vector3:new{1.0, 0.0}
c = Vector3:new{0.0, 1.0}
d = Vector3:new{1.0, 1.0}

nx0 = 16
ny0 = 16

grid = StructuredGrid:new{psurface=makePatch{Line:new{c,d}, Line:new{b,d}, Line:new{a,b}, Line:new{a,c}},
			  niv=nx0+1, njv=ny0+1}

blk = SBlock:new{grid=grid, fillCondition=initial, label="blk"}
blk.bcList[north] = ExtrapolateOutBC:new{xOrder=1}
blk.bcList[east] = ExtrapolateOutBC:new{xOrder=1}
blk.bcList[south] = UserDefinedBC:new{fileName='udf-bc.lua'}
blk.bcList[west] = UserDefinedBC:new{fileName='udf-bc.lua'}

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

				 
		
