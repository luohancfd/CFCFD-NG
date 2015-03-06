-- prep.lua
-- A place to put helper functions and classes for the preparation activities.
-- 
print("Loading prep.lua...")

-- Storage for later definitions of Block objects
blocks = {}

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
   for _,face in ipairs(faceList(config.dimensions)) do
      o.bcList[face] = o.bcList[face] or SlipWallBC:new()
   end
   o.hcellList = o.hcellList or {}
   o.xforceList = o.xforceList or {}
   -- Extract some information from the StructuredGrid
   -- Note 0-based indexing for vertices and cells.
   o.nic = o.grid:get_niv() - 1
   o.njc = o.grid:get_njv() - 1
   if config.dimensions == 3 then
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
      if config.dimensions == 3 then
	 str = str .. string.format('    "history-cell-%d": [%d, %d, %d],\n', 
				    i-1, hcell[1], hcell[2], hcell[3])
      else
	 str = str .. string.format('    "history-cell-%d": [%d, %d],\n',
				    i-1, hcell[1], hcell[2])
      end
   end
   -- Boundary conditions
   for _,face in ipairs(faceList(config.dimensions)) do
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
   f:write(string.format('"interpolation_order": %d,\n', config.interpolation_order))
   f:write(string.format('"gasdynamic_update_scheme": "%s",\n',
			 config.gasdynamic_update_scheme))
   f:write(string.format('"separate_update_for_viscous_terms": %s,\n',
			 tostring(config.separate_update_for_viscous_terms)))
   --f:write(string.format('"implicit_flag": %s,\n', tostring(config.implicit_flag)))
   f:write(string.format('"dt_init": %e,\n', config.dt_init))
   f:write(string.format('"dt_max": %e,\n', config.dt_max))
   f:write(string.format('"cfl_value": %e,\n', config.cfl_value))
   f:write(string.format('"stringent_cfl": %s,\n', tostring(config.stringent_cfl)))
   f:write(string.format('"fixed_time_step": %s,\n', tostring(config.fixed_time_step)))
   f:write(string.format('"dt_reduction_factor": %e,\n', config.dt_reduction_factor))
   f:write(string.format('"print_count": %d,\n', config.print_count))
   f:write(string.format('"cfl_count": %d,\n', config.cfl_count))
   f:write(string.format('"max_time": %e,\n', config.max_time))
   f:write(string.format('"max_step": %d,\n', config.max_step))
   f:write(string.format('"dt_plot": %e,\n', config.dt_plot))
   f:write(string.format('"dt_history": %e,\n', config.dt_history))
   f:write(string.format('"write_at_step": %d,\n', config.write_at_step))
   f:write(string.format('"halt_now": %d\n', config.halt_now))
   -- Note, also, no comma on last entry in JSON object.
   f:write("}\n")
   f:close()
end

function write_config_file(fileName)
   local f = assert(io.open(fileName, "w"))
   f:write("{\n")
   f:write(string.format('"title": "%s",\n', config.title))
   f:write(string.format('"dimensions": %d,\n', config.dimensions))
   f:write(string.format('"axisymmetric": %s,\n',
			 tostring(config.axisymmetric)))
   f:write(string.format('"gas_model_file": "%s",\n', config.gas_model_file))

   f:write(string.format('"thermo_interpolator": "%s",\n', 
			 string.lower(config.thermo_interpolator)))
   f:write(string.format('"interpolate_in_local_frame": %s,\n', 
			 tostring(config.interpolate_in_local_frame)))
   f:write(string.format('"apply_limiter": %s,\n', tostring(config.apply_limiter)))
   f:write(string.format('"extrema_clipping": %s,\n', tostring(config.extrema_clipping)))

   f:write(string.format('"flux_calculator": "%s",\n', config.flux_calculator))
   f:write(string.format('"compression_tolerance": %e,\n', config.compression_tolerance))
   f:write(string.format('"shear_tolerance": %e,\n', config.shear_tolerance))
   f:write(string.format('"M_inf": %e,\n', config.M_inf))
   f:write(string.format('"viscous": %s,\n', tostring(config.viscous)))

   f:write(string.format('"turbulence_model": "%s",\n',
			 string.lower(config.turbulence_model)))
   f:write(string.format('"turbulence_prandtl_number": %g,\n',
			 config.turbulence_prandtl_number))
   f:write(string.format('"turbulence_schmidt_number": %g,\n',
			 config.turbulence_schmidt_number))
   f:write(string.format('"max_mu_t_factor": %e,\n', config.max_mu_t_factor))
   f:write(string.format('"transient_mu_t_factor": %e,\n', config.transient_mu_t_factor))
   f:write(string.format('"separate_update_for_k_omega_source": %s,\n', 
			 tostring(config.separate_update_for_k_omega_source)))

   f:write(string.format('"reacting": %s,\n', tostring(config.reacting)))
   f:write(string.format('"max_invalid_cells": %d,\n', config.max_invalid_cells))
   f:write(string.format('"control_count": %d,\n', config.control_count))

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
   f:write(string.format("%04d %12.6e %12.6e\n", 0, 0.0, config.dt_init))
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
