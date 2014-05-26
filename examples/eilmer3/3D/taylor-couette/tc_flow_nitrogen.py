# taylor_couette.py
# Jason (Kan) Qin, December 2013

from math import pi, sin, cos, sqrt

gdata.dimensions = 3
gdata.title = "taylor couette flow"
print gdata.title

select_gas_model(model='ideal gas', species=['N2'])
gdata.viscous_flag = 1 
gdata.turbulence_model = "k_omega"
gdata.flux_calc = ADAPTIVE
gdata.max_time = 500.0e-3           # seconds
gdata.max_step = 400000
gdata.dt = 1.0e-11    
gdata.dt_plot = 1.0e-4 
gdata.dt_history = 1.0e-4

# Define flow conditions
p_exit = 1000 ;
r_omega = 2*pi*27600.0/60.0 ;
T_1 = 351.0 ;
T_2 = 366.0 ;
theta = 60.0*pi/180 ;

# Geometry
r_1 = 0.2125 ;
g_width = 0.0031 ;
r_2 = r_1 + g_width ;
h_1 = 0.0 ;
h_2 = 10.0*g_width ;

initial = FlowCondition(p=p_exit, u=0.0, v=0.0, w=0.0, T=T_1)

def makeSimpleBox(ini_angular1, ini_angular2, ini_h1, ini_h2):
   from math import pi, sin, cos
   inih1 = ini_h1 ;
   inih2 = ini_h2 ;
   ini1 = ini_angular1 ;
   ini2 = ini_angular2 ;
   center_b = Node(0.0, 0.0, inih1)
   center_t = Node(0.0, 0.0, inih2)
   p0 = Vector(r_1*cos(ini1), r_1*sin(ini1), inih1)
   p1 = Vector(r_2*cos(ini1), r_2*sin(ini1), inih1)
   p2 = Vector(r_2*cos(ini2), r_2*sin(ini2), inih1)
   p3 = Vector(r_1*cos(ini2), r_1*sin(ini2), inih1)
   p4 = Vector(r_1*cos(ini1), r_1*sin(ini1), inih2)
   p5 = Vector(r_2*cos(ini1), r_2*sin(ini1), inih2)
   p6 = Vector(r_2*cos(ini2), r_2*sin(ini2), inih2)
   p7 = Vector(r_1*cos(ini2), r_1*sin(ini2), inih2)
   p01 = Line(p0, p1) 
   p12 = Arc(p1, p2, center_b)
   p32 = Line(p3, p2)
   p03 = Arc(p0, p3, center_b)
   p45 = Line(p4, p5)
   p56 = Arc(p5, p6, center_t)
   p76 = Line(p7, p6)
   p47 = Arc(p4, p7, center_t)
   p04 = Line(p0, p4)
   p15 = Line(p1, p5)
   p26 = Line(p2, p6)
   p37 = Line(p3, p7)
   return WireFrameVolume(p01, p12, p32, p03, p45, p56, p76, p47, p04, p15, p26, p37)    
  
nx = 25 ; ny = 80 ; nz = 200 ;
nbx = 1 ; nby = 40 ; nbz = 1 ;  
c_0 = RobertsClusterFunction(1,1,1.0)

# North, East, South, West, Top, Bottom
mv = MovingWallBC(r_omega=[0.0,0.0,r_omega],v_trans=[0.0,0.0,0.0],Twall_flag=True,Twall=T_1)
ft = FixedTBC(Twall=T_2)
slip = SlipWallBC()

pvolume0 = makeSimpleBox(0.0*pi/180, 60.0*pi/180, h_1, h_2)
bclist0 = [None,ft,None,mv,slip,slip]
cflist0 = [c_0,]*12 ; 
blk0 = SuperBlock3D(label="check", nni=nx, nnj=ny, nnk=nz, 
               nbi=nbx, nbj=nby, nbk=nbz,
               parametric_volume=pvolume0,
               bc_list=bclist0,
               cf_list=cflist0,
               fill_condition=initial)

# South and North
connect_blocks_3D(blk0.blks[0][0][0],blk0.blks[0][-1][0],[(1,2),(5,6),(4,7),(0,3)],
                  reorient_vector_quantities=True,
                  nA=[0.0,1.0,0.0], t1A=[1.0,0.0,0.0],
                  nB=[-sin(theta),cos(theta),0.0], t1B=[cos(theta),sin(theta),0.0],
                  check_corner_locations=False)


identify_block_connections()




