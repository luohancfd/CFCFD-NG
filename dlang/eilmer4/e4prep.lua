-- e4prep.lua
-- A place to put helper functions and classes for the preparation activities.
-- 
print("Loading e4prep.lua...")
-- --------------------------------------------------------------------------

gdata = {
   title = "An Eilmer4 Simulation.",
   dimensions = 2,
   axisymmetric_flag = false,

   gas_model_file = "gas_model.lua",

   gasdynamic_update_scheme = "predictor-corrector",
   x_order = 2,
   interpolation_type = "rhoe",
   interpolate_in_local_frame = true,
   apply_limiter_flag = true,
   extrema_clipping_flag = true,
   flux_calc = "adaptive",
   compression_tolerance = -0.30,
   shear_tolerance = 0.20,
   M_inf = 0.01,
   cfl = 0.5,
   stringent_cfl = false,
   fixed_time_step = false,
   dt_reduction_factor = 0.2,
   t0 = 0.0, -- may be useful to change t0 if we are restarting from another job
   dt = 1.0e-6,
   dt_max = 1.0e-3,

   viscous_flag = false,
   separate_update_for_viscous_flag = false,
   implicit_flag = false,

   turbulence_model = "none",
   turbulence_prandtl_number = 0.89,
   turbulence_schmidt_number = 0.75,
   max_mu_t_factor = 300.0,
   transient_mu_t_factor = 1.0,
   separate_update_for_k_omega_source = false,

   reacting_flag = false,
   dt_chem = -1.0,

   print_count = 20,
   control_count = 10,
   cfl_count = 10,
   max_invalid_cells = 10,
   max_time = 1.0e-3,
   max_step = 10,
   dt_plot = 1.0e-3,
   dt_history = 1.0e-3,
   write_at_step = 0,
   halt_now = 0,
} -- end gdata

-- -------------------------------------------------------------------------------

-- Symbolic names for identifying boundaries
north = "north"; NORTH = "north"
east = "east"; EAST = "east"
south = "south"; SOUTH = "south"
west = "west"; WEST = "west"
top = "top"; TOP = "top"
bottom = "bottom"; BOTTOM = "bottom"

-- Class for BoundaryCondition
-- For the classes below, we just follow the prototype pattern
-- given in Ierusalimchy's book "Programming in Lua"

BoundaryCondition = {
   label = "",
   myType = ""
}

function BoundaryCondition:new (o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

SlipWallBC = BoundaryCondition:new()
SlipWallBC.myType = "SlipWall"

SupInBC = BoundaryCondition:new()
SupInBC.myType = "SupIn"
SupInBC.flowCondition = nil

ExtrapolateOutBC = BoundaryCondition:new()
ExtrapolateOutBC.myType = "ExtrapolateOut"

FullFaceExchangeBC = BoundaryCondition:new()
FullFaceExchangeBC.myType = "FullFaceExchange"
FullFaceExchangeBC.otherBlock = nil
FullFaceExchangeBC.otherFace = nil

-- Class for ClusterFunction
UnivariateFunction = {
   myType = ""
}

function UnivariateFunction:new (o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

LinearFunction = UnivariateFunction:new()
LinearFunction.t0 = 0.0
LinearFunction.t1 = 1.0
LinearFunction.myType = "Linear"

RobertsFunction = UnivariateFunction:new()
RobertsFunction.end0 = true
RobertsFunction.end1 = false
RobertsFunction.beta = 1.1
RobertsFunction.myType = "Roberts"


-- Class for Block construction (for a StructuredGrid).
SBlock = {
   psurf = nil, -- ParametricSurface object for 2D simulation
   pvol = nil,  -- ParametricVolume object for 3D simulation
   grid = nil,  -- The StructuredGrid object
   nic = nil,   -- number of cells i-index direction
   njc = nil,   -- number of cells j-index direction
   nkc = nil,   -- number of cells k-index direction
   fillCondition = nil, -- expects a FlowState object
   bcList = {
      [north]=SlipWallBC:new(), [east]=SlipWallBC:new(), 
      [south]=SlipWallBC:new(), [west]=SlipWallBC:new(),
      [top]=SlipWallBC:new(), [bottom]=SlipWallBC:new()
   },
   cfList = {}, -- Clustering functions
   label = "",
   hcellList = {},
   xforceList = {},
   myType = "SBlock",
   -- The following names are for the corner locations,
   -- to be used for testing for block connections.
   p000 = nil, p100 = nil, p110 = nil, p010 = nil,
   p001 = nil, p101 = nil, p111 = nil, p011 = nil,
} -- end Block

blocks = {}

function SBlock:new (o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   self.id = #blocks
   blocks[self.id+1] = self
   return o
end

-- ---------------------------------------------------------------------------

function connectBlocks(blkA, faceA, blkB, faceB, orientation)
   blkA.bcList[faceA] = FullFaceExchangeBC:new{otherBlock=blkB.id, otherFace=faceB,
					       orientation=orientation}
   blkB.bcList[faceB] = FullFaceExchangeBC:new{otherBlock=blkA.id, otherFace=faceA,
					       orientation=orientation}
   -- [TODO] need to test for matching corner locations and 
   -- concistent numbers of cells
end

function identifyBlockConnections()
   -- [TODO] identify block connections by trying to match corner points.
   -- Hard-code the cone20 connection for the moment.
   connectBlocks(blocks[1], west, blocks[2], east, 0)
end

-- ---------------------------------------------------------------------------

function write_control_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"x_order": %d,\n', gdata.x_order))
   f:write(string.format('"gasdynamic_update_scheme": "%s",\n',
			 gdata.gasdynamic_update_scheme))
   f:write(string.format('"separate_update_for_viscous_flag": %s,\n',
			 tostring(gdata.separate_update_for_viscous_flag)))
   f:write(string.format('"implicit_flag": %s,\n', tostring(gdata.implicit_flag)))
--[[
   dt = 1.0e-6,
   dt_max = 1.0e-3,
   cfl = 0.5,
   stringent_cfl = false,
   fixed_time_step = false,
   dt_reduction_factor = 0.2,
   print_count = 20,
   cfl_count = 10,
   max_time = 1.0e-3,
   max_step = 10,
   dt_plot = 1.0e-3,
   dt_history = 1.0e-3,
   write_at_step = 0,
   halt_now = 0,
--]]
   f:write('"last_entry": 0\n') -- no comma on last entry
   f:write("}\n")
   f:close()
end

function write_config_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"title": "%s",\n', gdata.title))
   f:write(string.format('"dimensions": %d,\n', gdata.dimensions))
   f:write(string.format('"axisymmetric_flag": %s,\n',
			 tostring(gdata.axisymmetric_flag)))
   f:write(string.format('"gas_model_file": "%s",\n', gdata.gas_model_file))
--[[
   interpolation_type = "rhoe",
   interpolate_in_local_frame = true,
   apply_limiter_flag = true,
   extrema_clipping_flag = true,
   flux_calc = "adaptive",
   compression_tolerance = -0.30,
   shear_tolerance = 0.20,
   M_inf = 0.01,
   t0 = 0.0, -- may be useful to change t0 if we are restarting from another job

   viscous_flag = false,
   turbulence_model = "none",
   turbulence_prandtl_number = 0.89,
   turbulence_schmidt_number = 0.75,
   max_mu_t_factor = 300.0,
   transient_mu_t_factor = 1.0,
   separate_update_for_k_omega_source = false,

   reacting_flag = false,
   dt_chem = -1.0,

   control_count = 10,
   max_invalid_cells = 10,
--]]
   for i,blk in pairs(blocks) do
      f:write("[TODO] write JSON config for block")
   end
   f:write('"last_entry": 0\n') -- no comma on last entry
   f:write("}\n")
   f:close()
end

function write_times_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# tindx sim_time dt_global\n");
   f:write(string.format("%04d %12.6e %12.6e\n", 0, gdata.t0, gdata.dt))
   f:close()
end

function write_block_list_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# indx label\n")
   for i,blk in pairs(blocks) do
      f:write("blk.indx blk.label [TODO] finish me")
   end
   f:close()
end

function build_job_files(jobName)
   print("Build job files for ", jobName)
   write_config_file(jobName .. ".config")
   write_control_file(jobName .. ".control")
   write_times_file(jobName .. ".times")
   write_block_list_file(jobName .. ".list")
end

print("Done loading e4prep.lua")
