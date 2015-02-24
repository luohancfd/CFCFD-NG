-- test.lua
-- Simple job-specification file for e4prep -- for use with Eilmer4
-- PJ & RG
-- 2015-02-24 -- adapted from the Python version of cone20

job_title = "Mach 1.5 flow over a 20 degree cone."
print(job_title)

-- We can set individual attributes of the global data object.
gdata.dimensions = 2
gdata.title = job_title
gdata.axisymmetric_flag = true

nsp, nmodes = setGasModel('sample-data/ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
initial = FlowState:new{p=5955.0, T=304.0, u=0.0, v=0.0}
inflow = FlowState:new{p=95.84e3, T=1103.0, u=1000.0, v=0.0}

-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.2, 1.0}
c = Vector3:new{1.0, 0.29118}
d = Vector3:new{1.0, 1.0}
e = Vector3:new{0.2, 1.0}
f = Vector3:new{0.0, 1.0}
ab = Line:new{a, b}; bc = Line:new{b, c} -- lower boundary including cone surface
fe = Line:new{f, e}; ed = Line:new{e, d} -- upper boundary
af = Line:new{a, f}; be = Line:new{b, e}; cd = Line:new{c, d} -- vertical lines

-- Define the blocks, with particular discretisation.
nx0 = 10; nx1 = 30; ny = 40
surf1 = makePatch{fe, be, ab, af}
surf2 = makePatch{ed, cd, bc, be, gridType="ao"}
--[[
blk_0 = Block2D(make_patch(fe, be, ab, af), nni=nx0, nnj=ny,
                fill_condition=inflow, label="BLOCK-0")
blk_1 = Block2D(make_patch(ed, cd, bc, be, "AO"), nni=nx1, nnj=ny,
                fill_condition=initial, label="BLOCK-1",
                hcell_list=[(9,0)], xforce_list=[0,0,1,0])

-- Set boundary conditions.
identify_block_connections()
blk_0.bc_list[WEST] = SupInBC(inflow, label="inflow-boundary")
blk_1.bc_list[EAST] = ExtrapolateOutBC(label="outflow-boundary")
--]]

-- Do a little more setting of global data.
gdata.max_time = 5.0e-3  -- seconds
gdata.max_step = 3000
gdata.dt = 1.0e-6
gdata.cfl = 0.5
-- gdata.dt_max = 10.0e-6
gdata.dt_plot = 1.5e-3
gdata.dt_history = 10.0e-5
