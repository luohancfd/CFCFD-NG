-- prep.lua
-- A place to put helper functions and classes for the preparation activities.
-- 
print("Loading prep.lua...")
-- --------------------------------------------------------------------------

blocks = {} -- storage for later definitions of Block objects

gdata = {
   title = "An Eilmer4 Simulation.",
   dimensions = 2,
   axisymmetric_flag = false,

   gas_model_file = "gas_model.lua",
   mhd_flag = false,
   radiation_flag = false,

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

function faceList(dimensions)
   local myList = {north, east, south, west}
   if dimensions == 3 then 
      table.insert(myList, top)
      table.insert(myList, bottom)
   end
   return myList
end

-- Class for BoundaryCondition
-- For the classes below, we just follow the prototype pattern
-- given in Ierusalimchy's book "Programming in Lua"

BoundaryCondition = {
   label = "",
   myType = ""
}
function BoundaryCondition:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   return o
end

SlipWallBC = BoundaryCondition:new()
SlipWallBC.myType = "SlipWall"
function SlipWallBC:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", ', self.label)
   str = str .. string.format('"bc": "%s"', self.myType)
   str = str .. '}'
   return str
end

SupInBC = BoundaryCondition:new{flowCondition=nil}
SupInBC.myType = "SupIn"
function SupInBC:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", ', self.label)
   str = str .. string.format('"bc": "%s", ', self.myType)
   str = str .. string.format('"inflow_condition": %s',
			      self.flowCondition:toJSONString())
   str = str .. '}'
   return str
end

ExtrapolateOutBC = BoundaryCondition:new{x_order=1}
ExtrapolateOutBC.myType = "ExtrapolateOut"
function ExtrapolateOutBC:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", ', self.label)
   str = str .. string.format('"bc": "%s", ', self.myType)
   str = str .. string.format('"x_order": %d', self.x_order)
   str = str .. '}'
   return str
end

FullFaceExchangeBC = BoundaryCondition:new{otherBlock=nil, otherFace=nil, orientation=0}
FullFaceExchangeBC.myType = "FullFaceExchange"
function FullFaceExchangeBC:tojson()
   local str = '{'
   str = str .. string.format('"label": "%s", ', self.label)
   str = str .. string.format('"bc": "%s", ', self.myType)
   str = str .. string.format('"other_block": %d, ', self.otherBlock)
   str = str .. string.format('"other_face": "%s", ', self.otherFace)
   str = str .. string.format('"orientation": %d', self.orientation)
   str = str .. '}'
   return str
end

-- Class for Block construction (based on a StructuredGrid).
SBlock = {
   myType = "SBlock",
   active = true,
   label = "",
   omegaz = 0.0,
   grid = nil,  -- The StructuredGrid object
   fillCondition = nil, -- expects a FlowState object
   bcList = nil, -- boundary conditions
   hcellList = nil,
   xforceList = nil,
   -- The following names are for the corner locations,
   -- to be used for testing for block connections.
   p000 = nil, p100 = nil, p110 = nil, p010 = nil,
   p001 = nil, p101 = nil, p111 = nil, p011 = nil,
} -- end Block

function SBlock:new(o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self
   -- Fill in default values, if already not set
   o.bcList = o.bcList or {}
   for _,face in ipairs(faceList(gdata.dimensions)) do
      o.bcList[face] = o.bcList[face] or SlipWallBC:new()
   end
   o.hcellList = o.hcellList or {}
   o.xforceList = o.xforceList or {}
   -- Extract some information from the StructuredGrid
   -- Note 0-based indexing for vertices and cells.
   o.nic = o.grid:get_niv() - 1
   o.njc = o.grid:get_njv() - 1
   if gdata.dimensions == 3 then
      o.nkc = o.grid:get_nkv() - 1
      o.p000 = o.grid:get_vtx(0, 0, 0)
      o.p100 = o.grid:get_vtx(o.nic, 0, 0)
      o.p110 = o.grid:get_vtx(o.nic, o.njc, 0)
      o.p010 = o.grid:get_vtx(0, o.nic, 0)
      o.p001 = o.grid:get_vtx(0, 0, o.nkc)
      o.p101 = o.grid:get_vtx(o.nic, 0, o.nkc)
      o.p111 = o.grid:get_vtx(o.nic, o.njc, o.nkc)
      o.p011 = o.grid:get_vtx(0, o.nic, o.nkc)
   else
      o.nkc = 1
      o.p000 = o.grid:get_vtx(0, 0)
      o.p100 = o.grid:get_vtx(o.nic, 0)
      o.p110 = o.grid:get_vtx(o.nic, o.njc)
      o.p010 = o.grid:get_vtx(0, o.nic)
   end
   -- Make a record of the new block, for later constructio of the config file.
   -- Note that we want block id to start at zero for the D code.
   o.id = #(blocks)
   blocks[#(blocks)+1] = o
   return o
end

function SBlock:tojson()
   local str = string.format('"block_%d": {\n', self.id)
   str = str .. string.format('    "label": "%s",\n', self.label)
   str = str .. string.format('    "active": %s,\n', tostring(self.active))
   str = str .. string.format('    "nic": %d,\n', self.nic)
   str = str .. string.format('    "njc": %d,\n', self.njc)
   str = str .. string.format('    "nkc": %d,\n', self.nkc)
   str = str .. string.format('    "omegaz": %f,\n', self.omegaz)
   str = str .. string.format('    "nhcell": %d,\n', #(self.hcellList))
   for i = 1, #(self.hcellList) do
      local hcell = self.hcellList[i]
      if gdata.dimensions == 3 then
	 str = str .. string.format('    "history-cell-%d": [%d, %d, %d],\n', 
				    i-1, hcell[1], hcell[2], hcell[3])
      else
	 str = str .. string.format('    "history-cell-%d": [%d, %d],\n',
				    i-1, hcell[1], hcell[2])
      end
   end
   -- Boundary conditions
   for _,face in ipairs(faceList(gdata.dimensions)) do
      str = str .. string.format('    "face_%s": ', face) ..
	 self.bcList[face]:tojson() .. ',\n'
   end
   str = str .. '    "dummy_entry_without_trailing_comma": 0\n'
   str = str .. '},\n'
   return str
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
   connectBlocks(blocks[1], east, blocks[2], west, 0)
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
   f:write(string.format('"dt": %e,\n', gdata.dt))
   f:write(string.format('"dt_max": %e,\n', gdata.dt_max))
   f:write(string.format('"cfl": %e,\n', gdata.cfl))
   f:write(string.format('"stringent_cfl": %s,\n', tostring(gdata.stringent_cfl)))
   f:write(string.format('"fixed_time_step": %s,\n', tostring(gdata.fixed_time_step)))
   f:write(string.format('"dt_reduction_factor": %e,\n', gdata.dt_reduction_factor))
   f:write(string.format('"print_count": %d,\n', gdata.print_count))
   f:write(string.format('"cfl_count": %d,\n', gdata.cfl_count))
   f:write(string.format('"max_time": %e,\n', gdata.max_time))
   f:write(string.format('"max_step": %d,\n', gdata.max_step))
   f:write(string.format('"dt_plot": %e,\n', gdata.dt_plot))
   f:write(string.format('"dt_history": %e,\n', gdata.dt_history))
   f:write(string.format('"write_at_step": %d,\n', gdata.write_at_step))
   f:write('"halt_now": 0\n') -- presumably, we want the simulation to proceed
   -- Note, also, no comma on last entry in JSON object.
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

   f:write(string.format('"interpolation_type": "%s",\n', 
			 string.lower(gdata.interpolation_type)))
   f:write(string.format('"interpolate_in_local_frame": %s,\n', 
			 tostring(gdata.interpolate_in_local_frame)))
   f:write(string.format('"apply_limiter_flag": %s,\n', tostring(gdata.apply_limiter_flag)))
   f:write(string.format('"extrema_clipping_flag": %s,\n',
			 tostring(gdata.extrema_clipping_flag)))

   f:write(string.format('"flux_calc": "%s",\n', gdata.flux_calc))
   f:write(string.format('"compression_tolerance": %e,\n', gdata.compression_tolerance))
   f:write(string.format('"shear_tolerance": %e,\n', gdata.shear_tolerance))
   f:write(string.format('"M_inf": %e,\n', gdata.M_inf))
   f:write(string.format('"viscous_flag": %s,\n', tostring(gdata.viscous_flag)))

   f:write(string.format('"turbulence_model": "%s",\n',
			 string.lower(gdata.turbulence_model)))
   f:write(string.format('"turbulence_prandtl_number": %g,\n',
			 gdata.turbulence_prandtl_number))
   f:write(string.format('"turbulence_schmidt_number": %g,\n',
			 gdata.turbulence_schmidt_number))
   f:write(string.format('"max_mu_t_factor": %e,\n', gdata.max_mu_t_factor))
   f:write(string.format('"transient_mu_t_factor": %e,\n', gdata.transient_mu_t_factor))
   f:write(string.format('"separate_update_for_k_omega_source": %s,\n', 
			 tostring(gdata.separate_update_for_k_omega_source)))

   f:write(string.format('"reacting_flag": %s,\n', tostring(gdata.reacting_flag)))
   f:write(string.format('"max_invalid_cells": %d,\n', gdata.max_invalid_cells))
   f:write(string.format('"control_count": %d,\n', gdata.control_count))

   f:write(string.format('"nblock": %d,\n', #(blocks)))
   for i = 1, #blocks do
      f:write(blocks[i]:tojson())
   end
   f:write('"dummy_entry_without_trailing_comma": 0\n') -- no comma on last entry
   f:write("}\n")
   f:close()
end

function write_times_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# tindx sim_time dt_global\n");
   f:write(string.format("%04d %12.6e %12.6e\n", 0, 0.0, gdata.dt))
   f:close()
end

function write_block_list_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("# indx label\n")
   for i = 1, #(blocks) do
      f:write(string.format("%4d %s\n", blocks[i].id, blocks[i].label))
   end
   f:close()
end

function build_job_files(job)
   print("Build job files for ", job)
   write_config_file(job .. ".config")
   write_control_file(job .. ".control")
   write_times_file(job .. ".times")
   write_block_list_file(job .. ".list")
   os.execute("mkdir -p grid/t0000")
   os.execute("mkdir -p flow/t0000")
   for i = 1, #blocks do
      local id = blocks[i].id
      print("Block id=", id)
      local fileName = "grid/t0000/" .. job .. string.format(".grid.b%04d.t0000", id)
      blocks[i].grid:write_to_text_file(fileName)
      os.execute("gzip -f " .. fileName)
      local fileName = "flow/t0000/" .. job .. string.format(".flow.b%04d.t0000", id)
      write_initial_flow_file(fileName, blocks[i].grid, blocks[i].fillCondition, 0.0)
      os.execute("gzip -f " .. fileName)
   end
   print("Done building files.")
end

print("Done loading e4prep.lua")
