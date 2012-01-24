# !/usr/bin/env python
# Name: testcase.py
#        Test model for the mixing plane boundary condition. 
#        Two models are available: a full annulus or a quarter annulus.
# Author: M.P. van der Laan
# Starting Date: 06/10/2009
# Last modification: 19/11/2009

# ------------------------------------------------
# ------------------Setup ------------------------
# ------------------------------------------------

job_directory = "/home/paulvdlaan/work/Eilmer3/testcase"
job_title = "Testcase for Radial Turbine flow with mixing plane"
print " "; print job_title; print " "

viscous_effect = 0

# Import library
from math import sin, cos, tan, pi
from numpy import *

# Global data
gdata.title = job_title
gdata.dimensions = 3
gdata.dt = 1.0e-6                   # Initial time step
gdata.flux_calc = ADAPTIVE          # AUSMDV.
gdata.t_order = 2
gdata.x_order = 2
gdata.cfl=0.5                       # 0.5 by default
gdata.viscous_flag = viscous_effect
gdata.max_time = 0.1                # Seconds
gdata.max_step = 4000               # Test runned with 4000 step
gdata.dt_plot = 0.0001              # Plot time
##gdata.apply_limiter_flag = 0       # Turn off limiter 

# Simulation setup
annulus     = "off"         # Annulus or quarter model: "on" or "off"   
flow        = "sub"         # Subsonic or supersonic flow: "sub" or "sup"
n           = 3             # Number of blocks; n=1: only stator; n=2:
                            # stator and rotor; n=3: stator, rotorand outlet 
mp_sr       = "on"          # Mixing Plane "on" or "off"
mp_ro       = "on"          # Mixing Plane "on" or "off"
periodic    = "on"          # Set periodic BC "on" or "off"
omega_z     = 2000.0        # Angular velocity rotor  
##def version():
##    save_input  = "051"     # Version number for saving data with testcase_output.py
##    return save_input

print ""
print "annulus model = ", annulus
print " "

   
# ------------------------------------------------
# -------------Flow conditions -------------------
# ------------------------------------------------
print ""
print "flow = ",flow
print ""
select_gas_model(model='ideal gas', species=['Ar']) # Argon

def initial_gas_state(x, y, z):
    global flow
    if flow == "sub":
        p_ini = 80.0e3      # Initial pressure (Pa)
        T_ini = 274.383031  # Initial temperature (K), value based on 1D
                            # calculation to speed up convergence.
        w_ini = 163.274995  # Initial velocity z-direction, value based on 1D 
                            # calculation to speed up convergence.
        w_x = 0.0
        w_y = 0.0
    elif flow == "sup":
        p_ini = 80.0e3      
        T_ini = 250.0      
        w_ini = 400.0      
        omega_ini = 2000.0      # rad/s
        r = sqrt(x*x + y*y)
        c_theta = omega_ini*r
        c_r = 0.0
        w_x = -y/r*c_theta + x/r*c_r 
        w_y =  x/r*c_theta + y/r*c_r  
    return FlowCondition(p=p_ini, u=w_x, v=w_y, w=w_ini, T=T_ini, add_to_list=0).to_dict()


# ------------------------------------------------
# -----------Quarter model------------------------
# ------------------------------------------------

if annulus == "off":    
    
    # -----------Geometrie flow domain -----------
    r1 = 0.01           # Radius of inner cylinder (hub)
    r2 = 0.03           # Radius of outer cylinder
    L = 0.05            # Length of one cylinder block    
    for i in range(1,n+1):
        # Coordinates
        o0 = Node(0,  0,  (i-1)*L)      # Origin
        o1 = Node(0,  0,  i*L)          # Origin xy-plane at z=L
        p0 = Node(0,  r1, (i-1)*L)
        p1 = Node(0,  r1, i*L)
        p2 = Node(r1, 0,  i*L)
        p3 = Node(r1, 0,  (i-1)*L)
        p4 = Node(0,  r2, (i-1)*L)
        p5 = Node(0,  r2, i*L) 
        p6 = Node(r2, 0,  i*L) 
        p7 = Node(r2, 0,  (i-1)*L)
        # Edges
        p01 = Line(p0, p1)
        p12 = Arc(p1, p2, o1)
        p32 = Line(p3, p2)
        p03 = Arc(p0, p3, o0)
        p45 = Line(p4, p5)
        p56 = Arc(p5, p6, o1)
        p76 = Line(p7, p6)
        p47 = Arc(p4, p7, o0)
        p04 = Line(p0, p4)
        p15 = Line(p1, p5)
        p26 = Line(p2, p6)
        p37 = Line(p3, p7)
        if i == 1:
            # Define stator
            volume = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
            stator = Block3D(label="stator", nni=10, nnj=30, nnk=30,
                        parametric_volume=volume, fill_condition=initial_gas_state)
        elif i == 2:   
            # Define rotor
            volume = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
            rotor = Block3D(label="rotor", nni=10, nnj=30, nnk=30,
                            parametric_volume=volume, fill_condition=initial_gas_state,
                            omegaz=omega_z)
        elif i == 3:   
            # Define outlet
            volume = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
            outlet = Block3D(label="outlet", nni=10, nnj=30, nnk=30,
                            parametric_volume=volume, fill_condition=initial_gas_state)

    # Connect blocks automatically
    identify_block_connections()
            
    # --------------Boundary conditions---------------    
    print ""
    print "Set boundary conditions:"
    print ""

    # Boundary conditions will override the identify_block_connections()
    
    # Supersonic
    if flow == "sup":
        # Supersonic inflow
        stator.set_BC(WEST, USER_DEFINED, filename="udf-supersonic-in-bc.lua")
        print "UDF supersonic inflow = on"
        # Extrapolate out and connect blocks in z direction
        if n == 1:
            stator.set_BC(EAST, EXTRAPOLATE_OUT)
        elif n == 2:
            rotor.set_BC(EAST, EXTRAPOLATE_OUT)
        elif n == 3:
            outlet.set_BC(EAST, EXTRAPOLATE_OUT)
    # Subsonic  
    elif flow == "sub":    
        # Subsonic inflow, using T0 and P0, alpha and beta
        print "UDF subsonic inflow = on"
        stator.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-in-bc.lua")
        # Mixing Plane rotor-stator
        if mp_sr == "on" and (n == 2 or n == 3):
            stator.set_BC(EAST, USER_DEFINED, filename="udf-mp-sr-up-bc.lua")
            rotor.set_BC(WEST, USER_DEFINED,  filename="udf-mp-sr-down-bc.lua")          
        # Mixing Plane rotor-outlet   
        if mp_ro == "on" and n == 3:
            rotor.set_BC(EAST, USER_DEFINED,  filename="udf-mp-ro-up-bc.lua")
            outlet.set_BC(WEST, USER_DEFINED, filename="udf-mp-ro-down-bc.lua")
        print "Mixing Plane =", mp_sr, "Mixing Plane outlet =", mp_ro
        print ""
        
        # Forced exit pressure BC
        print "Forced exit pressure = on"
        if n == 1:
            stator.set_BC(EAST, USER_DEFINED, filename="udf-subsonic-out-bc.lua")
        elif n == 2:
            rotor.set_BC(EAST, USER_DEFINED, filename="udf-subsonic-out-bc.lua")
        elif n == 3:
            outlet.set_BC(EAST, USER_DEFINED, filename="udf-subsonic-out-bc.lua") 
      
    # Perdiodic BC's
    if periodic == "on":
        print "Periodic BSs =", periodic
        stator.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")    
        stator.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")
        if n == 2 or n == 3:
            rotor.set_BC(SOUTH,  USER_DEFINED, filename="udf-periodic-bc.lua")
            rotor.set_BC(NORTH,  USER_DEFINED, filename="udf-periodic-bc.lua")
        if n == 3:
            outlet.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")
            outlet.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua") 
    elif periodic == "off":
        print "Periodic BSs =", periodic
        
# ------------------------------------------------
# -----------Annulus model------------------------
# ------------------------------------------------

elif annulus == "on":
    
    # -----------Geometrie flow domain -----------
    r1 = 0.01           # Radius of inner cylinder (hub)
    r2 = 0.03           # Radius of outer cylinder
    L = 0.05            # Length of one cylinder block
    for i in range(1,n+1):
        # Cylinder               
        # Coordinates
        o0 = Node(0.0,  0.0,  (i-1)*L)  
        o1 = Node(0.0,  0.0,  i*L)
        p0 = Node(0.0,  r1,   (i-1)*L)
        p1 = Node(0.0,  r1,   i*L)
        p2 = Node(r1,   0.0,  i*L)
        p3 = Node(r1,   0.0,  (i-1)*L)
        p4 = Node(0.0,  r2,   (i-1)*L)
        p5 = Node(0.0,  r2,   i*L) 
        p6 = Node(r2,   0.0,  i*L) 
        p7 = Node(r2,   0.0,  (i-1)*L)
        #Edges
        p01 = Line(p0, p1)
        p12 = Arc(p1, p2, o1)
        p32 = Line(p3, p2)
        p03 = Arc(p0, p3, o0)
        p45 = Line(p4, p5)
        p56 = Arc(p5, p6, o1)
        p76 = Line(p7, p6)
        p47 = Arc(p4, p7, o0)
        p04 = Line(p0, p4)
        p15 = Line(p1, p5)
        p26 = Line(p2, p6)
        p37 = Line(p3, p7)
        if i == 1:                          
            volume_s0 = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
            stator = []
            for j in range(0,4):
                new_volume = volume_s0.clone()
                new_volume.rotate_about_zaxis(j*pi/2)
                blk = Block3D(label="stator%d" % j,
                              nni=10, nnj=15, nnk=2*15,
                              parametric_volume=new_volume, 
                              fill_condition=initial_gas_state)                          
                stator.append(blk)             
        elif i == 2:
            volume_r0 = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
            rotor = []
            for j in range(0,4):
                new_volume = volume_s0.clone()
                new_volume.rotate_about_zaxis(j*pi/2)
                blk = Block3D(label="rotor%d" % j,
                              nni=10, nnj=15, nnk=15,
                              parametric_volume=new_volume, 
                              fill_condition=initial_gas_state,
                              omegaz = omega_z)                          
                rotor.append(blk)                            
        elif i == 3:
            volume_o0 = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
            outlet = []
            for j in range(0,4):
                new_volume = volume_o0.clone()
                new_volume.rotate_about_zaxis(j*pi/2)
                blk = Block3D(label="outlet%d" % j,
                              nni=10, nnj=15, nnk=15,
                              parametric_volume=new_volume, 
                              fill_condition=initial_gas_state)                          
                outlet.append(blk)
                
    #Connect blocks automatically
    identify_block_connections()
  
##    south_to_north = [(0,3),(4,7),(5,6),(1,2)]
##    connect_blocks_3D(A = stator[0], B = stator[1], vtx_pairs = south_to_north)
##    connect_blocks_3D(A = stator[1], B = stator[2], vtx_pairs = south_to_north)
##    connect_blocks_3D(A = stator[2], B = stator[3], vtx_pairs = south_to_north)
##    connect_blocks_3D(A = stator[3], B = stator[0], vtx_pairs = south_to_north)
    # --------------Boundary conditions---------------
    print ""
    print "Set boundary conditions:"
    print ""

    # Boundary conditions will override the identify_block_connections()
    
    for j in range(0,4):
        # Supersonic
        if flow == "sup":
            # Supersonic inflow
            stator[j].set_BC(WEST, USER_DEFINED, filename="udf-supersonic-in-bc.lua")
            print "UDF supersonic inflow = on"
            # Extrapolate out
            if n == 1:
                stator[j].set_BC(EAST, EXTRAPOLATE_OUT)            
            elif n == 2:
                rotor[j].set_BC(EAST, EXTRAPOLATE_OUT)
         
            elif n == 3:
                outlet[j].set_BC(EAST, EXTRAPOLATE_OUT)              
        # Subsonic  
        elif flow == "sub":    
            # Subsonic inflow, using T0 and P0, alpha and beta
            print "UDF subsonic inflow = on"
            stator[j].set_BC(WEST, USER_DEFINED, filename="udf-subsonic-in-bc.lua")
            # Mixing Plane stator-rotor
            if mp_sr == "on" and (n == 2 or n == 3):
                stator[j].set_BC(EAST, USER_DEFINED, filename="udf-mp-annulus-sr-up-bc.lua")
                rotor[j].set_BC(WEST, USER_DEFINED,  filename="udf-mp-annulus-sr-down-bc.lua")   
            # Mixing Plane rotor-outlet   
            if mp_ro == "on" and n == 3:
                rotor[j].set_BC(EAST, USER_DEFINED,  filename="udf-mp-annulus-ro-up-bc.lua")         
                outlet[j].set_BC(WEST, USER_DEFINED, filename="udf-mp-annulus-ro-down-bc.lua")
            print "Mixing Plane =", mp_sr, "Mixing Plane outlet =", mp_ro
            print ""
            
            ### Forced exit pressure BC
            print "Forced exit pressure = on"
            if n == 1:
                stator[j].set_BC(EAST, USER_DEFINED, filename="udf-subsonic-out-bc.lua")
            elif n == 2:
                rotor[j].set_BC(EAST, USER_DEFINED, filename="udf-subsonic-out-bc.lua")
            elif n == 3:
                outlet[j].set_BC(EAST, USER_DEFINED, filename="udf-subsonic-out-bc.lua")
   
print ""
print "done"
