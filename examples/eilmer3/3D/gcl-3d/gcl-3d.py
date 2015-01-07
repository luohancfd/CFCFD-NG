# gcl-3d.py
# Jason (Kan) Qin, Jan 2015

gdata.dimensions = 3
select_gas_model(model='ideal gas', species=['air'])

initial =  FlowCondition(u=0.0, v=0.0, w=0.0)

def makeSimpleBox(ini_x0, ini_x1, ini_y0, ini_y1, ini_z0, ini_z1):
   x0 = ini_x0 ; x1 = ini_x1 ;
   y0 = ini_y0 ; y1 = ini_y1 ;
   z0 = ini_z0 ; z1 = ini_z1 ;
   p0 = Vector(x0, y0, z0)
   p1 = Vector(x1, y0, z0)
   p2 = Vector(x1, y1, z0)
   p3 = Vector(x0, y1, z0)
   p4 = Vector(x0, y0, z1)
   p5 = Vector(x1, y0, z1)
   p6 = Vector(x1, y1, z1)
   p7 = Vector(x0, y1, z1)
   p01 = Line(p0, p1) 
   p12 = Line(p1, p2)
   p32 = Line(p3, p2)
   p03 = Line(p0, p3)
   p45 = Line(p4, p5)
   p56 = Line(p5, p6)
   p76 = Line(p7, p6)
   p47 = Line(p4, p7)
   p04 = Line(p0, p4)
   p15 = Line(p1, p5)
   p26 = Line(p2, p6)
   p37 = Line(p3, p7)
   return WireFrameVolume(p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37)  
   
# Geometry
edge = 0.03  

# Define the blocks, boundary conditions and set the discretisaztion.
nx0 = 5 ; ny0 = 5 ; nz0 = 5 ; 
c_0 = RobertsClusterFunction(1,1,1.0)
pvolume0 = makeSimpleBox(0.0,edge,0.0,edge,0.0,edge)
cflist0 = [c_0,]*12 ; 
blk = Block3D(label="plate", nni=nx0, nnj=ny0, nnk=nz0,
               parametric_volume=pvolume0,
               cf_list=cflist0,
               fill_condition=initial)
               
blk.bc_list[NORTH] = SlipWallBC()
blk.bc_list[EAST] = SlipWallBC()
blk.bc_list[SOUTH] = SlipWallBC()
blk.bc_list[WEST] = SlipWallBC()
blk.bc_list[TOP] = SlipWallBC()
blk.bc_list[BOTTOM] = SlipWallBC()               

identify_block_connections()

gdata.title = "random grid motion"
gdata.viscous_flag = 0
gdata.flux_calc = AUSMDV
gdata.max_time = 10e-3
gdata.max_step = 10000
gdata.dt = 1.0e-6
gdata.dt_plot = 1e-3
gdata.dt_history = 1e-3
gdata.gasdynamic_update_scheme = "pc"
print "The gas-dynamic update scheme is", gdata.gasdynamic_update_scheme
gdata.x_order = 2
gdata.interpolate_in_local_frame = True

# moving grid flag
gdata.dt_moving = 0
gdata.moving_grid_flag = 1

# udf file to prescribe vertex velocity
gdata.udf_file = "udf-vtx.lua"
gdata.udf_vtx_velocity_flag = 1

# The following scales provide a reasonable picture.
sketch.xaxis(0.0, 0.040, 0.01, -0.005)
sketch.yaxis(0.0, 0.010, 0.01, -0.004)
sketch.window(0.0, 0.0, 0.040, 0.010, 0.05, 0.05, 0.15, 0.075)

