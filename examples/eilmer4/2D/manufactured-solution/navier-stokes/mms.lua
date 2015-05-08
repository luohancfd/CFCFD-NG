-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-05-05
--

title = "Method of Manufactured Solutions, Navier-Stokes case."
print(title)
config.dimensions = 2
config.title = title

setGasModel('very-viscous-air.lua')

R = 287.0

p0 = 1.0e5
T0 = p0/R
u0 = 70.0
v0 = 90.0

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
blk.bcList[north] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
blk.bcList[east] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
blk.bcList[south] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
blk.bcList[west] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}

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

				 
		
