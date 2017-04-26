## \Sliding_Mesh_test.py
#
""" 
Script Example Case for setting up a sliding mesh interface in Eilmer 3.
The interface passes convective fluxes, however viscous effects are ignored at 
the interface.

The sliding mesh interface is created by four udf files:
udf-config_e3.lua
    This file is used to set-up the sliding mesh interface. See details in the 
    file. During a normal slidign mesh simualtuion you would only ned to modify 
    this file. The other lua files will update automatcially. The exception are 
    additional functions tha may be executed as part of the at_timestep_start() 
    or at_timestep_end() routines. These are set in udf-process_e3.lua

udf-process_e3.lua
    This file contains thew function at_timestep_start(args) and 
    at_timestep_end(args), which are called before and after each timestp. These 
    functions are used to complete the following functions.
 (1) At first timestep, set-up inforamtion for the sliding interface is 
    generated and saved in the two files UP_tbl.lua and DOWN_tbl.lua. Details on
    the setup of this 
 (2) Every 50 times-steps the angular position that has been reached is reported.
 (3) Every 100 time-setps the mass flow rate across the mesh inlet and outlet is 
    reported.

udf-rotor_in_e3.lua
    This file contains the convective_flux(args) function, which is used to 
    define the moving mesh boundary condition on the rotor inlet face.

udf-stator_out_e3.lua
    This file contains the convective_flux(args) function, which is used to 
    define the moving mesh boundary condition on the stator outlet face.


Author: Ingo Jahn
Last modified: 26/04/2017 
"""

import numpy as np

####################################
### Setting up Basic Information ###
####################################
# For grid development, set gdata.dimensions = 2, this will create the 2-D projection of the mesh.
gdata.dimensions = 3
gdata.axisymmetric_flag = 0

# Set some fluid propertied to allow e3prep to solve
# These only need to be correct if using Eilmer as solver.
select_gas_model(model='ideal gas', species=['CO2'])

# Do a little more setting of global data.
gdata.max_time = 20.e-4  # seconds
gdata.max_step = 100000
gdata.dt = 1.0e-8
gdata.dt_plot = 1.e-5
#gdata.dt_history = 10.0e-5

# call lua to run the "udf-config_e3.lua" script to assist with the set up of the BC.
import os
os.system('lua udf-config_e3.lua')

if 0:
    initial1 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    initial2 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    initial3 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    stagnation1 = FlowCondition(p=1e3, u=0., T=800.)
    stagnation2 = FlowCondition(p=1e3, u=0., T=800.)

if 0:
    initial1 = FlowCondition(p=0.99e5,  u=0.,    v=0.0, T=200.)
    initial2 = FlowCondition(p=0.99e5,  u=0.,    v=0.0, T=200.)
    initial3 = FlowCondition(p=0.99e5,  u=0.,    v=0.0, T=200.)
    stagnation1 = FlowCondition(p=1.0e5, u=0., T=833.)
    stagnation2 = FlowCondition(p=1.0e5, u=0., T=400.)
    Pout = 0.99e5

if 0:
    initial1 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    initial2 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    initial3 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    stagnation1 = FlowCondition(p=1e5, u=-100., T=800.)

if 0:
    initial1 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    initial2 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    initial3 = FlowCondition(p=1e3,  u=0.,    v=0.0, T=100.)
    stagnation1 = FlowCondition(p=1e5, u=-600., T=400.)
    stagnation2 = FlowCondition(p=1e5, u=-600., T=500.)


if 1:
    initial1 = FlowCondition(p=1e5,  u=0.,    v=0.0, T=100.)
    initial2 = FlowCondition(p=1e5,  u=0.,    v=0.0, T=100.)
    initial3 = FlowCondition(p=1e5,  u=0.,    v=0.0, T=100.)
    stagnation1 = FlowCondition(p=1.1e5, T=500.)
    stagnation2 = FlowCondition(p=1e5, T=500.)
    Pout = 0.99e5


# Create Stator and Rotor Mesh 
# Define Geometry
R1 = 0.05
R2 = 0.045 
R3 = 0.04
theta = 2.*np.pi/16  # 16 blades
theta_a = -0.5*theta 
theta_b = 0.0*theta
theta_c = 0.5*theta

C = Node(0.,0.,0.,label="C")
R1a = Node(R1*np.cos(theta_a),R1*np.sin(theta_a),0.,label="R1a")
R1b = Node(R1*np.cos(theta_b),R1*np.sin(theta_b),0.,label="R1b")
R1c = Node(R1*np.cos(theta_c),R1*np.sin(theta_c),0.,label="R1c")
R2a = Node(R2*np.cos(theta_a),R2*np.sin(theta_a),0.,label="R2a")
R2b = Node(R2*np.cos(theta_b),R2*np.sin(theta_b),0.,label="R2b")
R2c = Node(R2*np.cos(theta_c),R2*np.sin(theta_c),0.,label="R2c")
R3a = Node(R3*np.cos(theta_a),R3*np.sin(theta_a),0.,label="R3a")
R3b = Node(R3*np.cos(theta_b),R3*np.sin(theta_b),0.,label="R3b")
R3c = Node(R3*np.cos(theta_c),R3*np.sin(theta_c),0.,label="R3c")

R1_1 = Arc(R1a,R1b,C); R1_2 = Arc(R1b,R1c,C)
R2_1 = Arc(R2a,R2b,C); R2_2 = Arc(R2b,R2c,C)
R3_1 = Arc(R3a,R3b,C); R3_2 = Arc(R3b,R3c,C)
R2R1a = Line(R2a,R1a); R2R1b = Line(R2b,R1b); R2R1c = Line(R2c,R1c)
R2R3a = Line(R2a,R3a); R2R3b = Line(R2b,R3b); R2R3c = Line(R2c,R3c)

H_1 = 0.01
H_2 = 0.02
H1 = Line(R2a, Node(R2a.x,R2a.y,R2a.z+H_2) )
H2 = Line(R2b, Node(R2b.x,R2b.y,R2b.z+H_2) )
H3 = Line(R2a, Node(R2a.x,R2a.y,R2a.z+H_1) )
H4 = Line(Node(R2a.x,R2a.y,R2a.z+H_1), Node(R2a.x,R2a.y,R2a.z+H_2) )

# Create Blocks
blk0 = Block3D(WireFrameVolume(make_patch(R2R1b,R1_1,R2R1a,R2_1), H1), 
                nni=10, nnj=6, nnk=7,
                fill_condition=initial1, label="blk0",omegaz=0.)
blk1 = Block3D(WireFrameVolume(make_patch(R2R1c,R1_2,R2R1b,R2_2), H2), 
                nni=10, nnj=6, nnk=7,
                fill_condition=initial2, label="blk1",omegaz=0.)

blk2 = Block3D(WireFrameVolume(make_patch(Polyline([R3_1,R3_2]), R2R3c, Polyline([R2_1,R2_2]),R2R3a), H3), 
                nni=10, nnj=10, nnk=5,
                fill_condition=initial3, label="blk2",omegaz=0.)
blk3 = Block3D(WireFrameVolume(make_patch(Polyline([R3_1,R3_2]), R2R3c, Polyline([R2_1,R2_2]),R2R3a), H4), 
                nni=10, nnj=10, nnk=4,
                fill_condition=initial3, label="blk3",omegaz=0.)


"""
Looking radially outwards at upstream patch faces
   +------+------+
   | blk1 | blk0 | ^
   |      |      | | Z,k
   +------+------+ |
              <----+
                 T,,j

Looking radially outwards at downstream patch faces
   +-------------+
   |    blk3     |  
   |             |  
   +-------------+     
   |    blk2     |  ^
   |             |  | Z,k
   +-------------+  |   
              <-----+ 
                 T,i
"""

# connect blocks
identify_block_connections()


# periodically connect downstream blocks.
connect_blocks_3D(blk2,blk2,[(0,1),(3,2),(7,6),(4,5)],
        reorient_vector_quantities=True, 
        nA=[0.,1.,0.],t1A=[1.,0.,0.],
        nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
        check_corner_locations=False)
connect_blocks_3D(blk3,blk3,[(0,1),(3,2),(7,6),(4,5)],
        reorient_vector_quantities=True, 
        nA=[0.,1.,0.],t1A=[1.,0.,0.],
        nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
        check_corner_locations=False)

# periodically connect upstream blocks.
connect_blocks_3D(blk0,blk1,[(1,2),(0,3),(4,7),(5,6)],
        reorient_vector_quantities=True, 
        nA=[0.,1.,0.],t1A=[1.,0.,0.],
        nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
        check_corner_locations=False)


# Set Boundary Connections
if 0:
    # supersonic flow in the radial inwards direction
    blk0.bc_list[EAST] = SupInBC(stagnation1,label='OF_inlet_00')
    blk1.bc_list[EAST] = SupInBC(stagnation2,label='OF_inlet_00')

    blk2.bc_list[NORTH] = ExtrapolateOutBC(label='OF_outlet_00')
    blk3.bc_list[NORTH] = ExtrapolateOutBC(label='OF_outlet_00')

if 0:
    # supersonic flow in the radial outwards direction
    blk0.bc_list[EAST] = ExtrapolateOutBC(label='OF_outlet_00')
    blk1.bc_list[EAST] = ExtrapolateOutBC(label='OF_outlet_00')

    blk2.bc_list[NORTH] = SupInBC(stagnation1,label='OF_inlet_00')
    blk3.bc_list[NORTH] = SupInBC(stagnation1,label='OF_inlet_00')

if 1:
    # subsonic flow in the radial inwards direction
    blk0.bc_list[EAST] = SubsonicInBC(stagnation1,label='OF_inlet_00') 
    blk1.bc_list[EAST] = SubsonicInBC(stagnation2,label='OF_inlet_00')

    blk2.bc_list[NORTH] = FixedPOutBC(Pout, label='OF_outlet_00')
    blk3.bc_list[NORTH] = FixedPOutBC(Pout, label='OF_outlet_00')


if 0:
    blk0.bc_list[EAST] = FixedPOutBC(Pout, label='OF_outlet_00') 
    blk1.bc_list[EAST] = FixedPOutBC(Pout, label='OF_outlet_00') 

    blk2.bc_list[NORTH] = SubsonicInBC(stagnation2,label='OF_inlet_00')
    blk3.bc_list[NORTH] = SubsonicInBC(stagnation2,label='OF_inlet_00')



#######################################
### Couple Meshes and Overwrite BCs ###
#######################################
if 0: # MappedCell_BC 
    # Stator Outlet faces: E0.bc_list[WEST]; E1.bc_list[WEST]; E2.bc_list[WEST]; E3.bc_list[WEST]; E4.bc_list[WEST]
    blk0.bc_list[WEST] = MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))
    blk1.bc_list[WEST] = MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))
   
    # Rotor Inlet faces: BLI_IT0.bc_list[SOUTH]; BLI_IT1.bc_list[SOUTH]; BLI_IB0.bc_list[SOUTH]; BLI_IB1.bc_list[SOUTH]; BLI_IC0.bc_list[SOUTH]; BLI_IC1.bc_list[SOUTH]; 
    blk2.bc_list[SOUTH] = MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))
    blk3.bc_list[SOUTH] = MappedCellBC(ghost_cell_trans_fn=lambda x, y, z: (x, y, z))

if 1: # Sliding Mesh BC
    print "#######################################################"
    print "Creating simulation environment for Mapped (sliding mesh) "
    print "         WARNING: FURTHER VALIDATION ADVISED"
    print "#######################################################"

    gdata.udf_file = "udf-process_e3.lua"

    # UserDefinedBC(filename, is_wall=0, sets_conv_flux=0, sets_visc_flux=0,label='')

    blk0.bc_list[WEST] = UserDefinedBC("udf-stator_out_e3.lua", sets_conv_flux=1, sets_visc_flux=0)
    blk1.bc_list[WEST] = UserDefinedBC("udf-stator_out_e3.lua", sets_conv_flux=1, sets_visc_flux=0)

    blk2.bc_list[SOUTH] = UserDefinedBC("udf-rotor_in_e3.lua", sets_conv_flux=1, sets_visc_flux=0)
    blk3.bc_list[SOUTH] = UserDefinedBC("udf-rotor_in_e3.lua", sets_conv_flux=1, sets_visc_flux=0)

    # by set_conv_flux=1 --> need convective_flux() in lua
    # by set_visc_flux=0 --> need interface() in lua  if using viscous simulation  


