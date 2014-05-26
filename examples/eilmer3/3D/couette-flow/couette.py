# couette.py
# Jason (Kan) Qin, November 2013

from math import pi, sin, cos

gdata.dimensions = 3
gdata.title = "pressure distribution in a thrust bearing chamber 3D"
print gdata.title

select_gas_model(model='ideal gas', species=['air'])
gdata.viscous_flag = 1 
gdata.turbulence_model = "k_omega"
gdata.flux_calc = ADAPTIVE
gdata.max_time = 5.0e-3 # seconds
gdata.max_step = 30000
gdata.dt = 1.0e-10     
gdata.dt_plot = 1.0e-3 
gdata.dt_history = 1.0e-3

# Define flow conditions
p_exit = 0.1e6
v_trans = 130.0 ;

# Geometry
h_1 = 0.0 ;
h_2 = 3e-3 ;
r_1 = 0.0  ;
r_2 = 1e-1 ;
l_1 = 0.0 ;
l_2 = 1.5e-1 ;

def initial_flow(x, y, z):
    global h_2, p_exit, v_trans
    v = v_trans * z / h_2 # linear velocity profile
    return FlowCondition(p=p_exit, u=0.0, v=v, w=0.0).to_dict()

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

# Define the blocks, boundary conditions and set the discretisaztion.
nx0 = 5 ; ny0 = 30 ; nz0 = 20 ; 
c_0 = RobertsClusterFunction(1,1,1.0)
pvolume0 = makeSimpleBox(r_1,r_2,l_1,l_2,h_1,h_2)
cflist0 = [c_0,]*12 ; 
blk_0 = Block3D(label="plate", nni=nx0, nnj=ny0, nnk=nz0,
               parametric_volume=pvolume0,
               cf_list=cflist0,
               fill_condition=initial_flow)

blk_0.set_BC("TOP", "MOVING_WALL", r_omega=[0.0,0.0,0.0],v_trans=[0.0,v_trans,0.0]) 
blk_0.bc_list[BOTTOM] = AdiabaticBC()
blk_0.set_BC("WEST","SLIP_WALL")
blk_0.set_BC("EAST","SLIP_WALL")

# the south face is connected with the north face
connect_blocks_3D(blk_0,blk_0,[(1,2),(5,6),(4,7),(0,3)],
                  reorient_vector_quantities=True,
                  nA=[0.0,1.0,0.0], t1A=[1.0,0.0,0.0],
                  nB=[0.0,1.0,0.0], t1B=[1.0,0.0,0.0],
                  check_corner_locations=False)

identify_block_connections()




