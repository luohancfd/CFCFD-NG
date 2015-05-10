-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-05-05
--

config.title = "Method of Manufactured Solutions, Navier-Stokes case."
print(config.title)
config.dimensions = 2

setGasModel('very-viscous-air.lua')
R = 287.0; p0 = 1.0e5; T0 = p0/R; u0 = 70.0; v0 = 90.0;
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0}

p00 = Vector3:new{0.0, 0.0}
p10 = Vector3:new{1.0, 0.0}
p01 = Vector3:new{0.0, 1.0}
p11 = Vector3:new{1.0, 1.0}
nicell = 16; njcell = 16
grid = StructuredGrid:new{psurface=CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01},
			  niv=nicell+1, njv=njcell+1}

bcList = {}
for _,face in ipairs{north, east, south, west} do
   bcList[face] = BoundaryCondition:new{
      preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
      preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
				UpdateThermoTransCoeffs:new() }
   }
end
blk = SBlockArray{grid=grid, fillCondition=initial, bcList=bcList, nib=2, njb=2, label="blk"}

config.interpolation_order = 2
config.gasdynamic_update_scheme = "predictor-corrector"
config.flux_calculator = 'ausmdv'
config.viscous = true
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.dt_init = 1.0e-7
config.max_time = 150.0e-3
config.dt_plot = config.max_time/20.0
config.max_step = 3000000
config.cfl = 0.5
config.stringent_cfl = 1
config.apply_limiter = false
config.extrema_clipping = false

				 
		
