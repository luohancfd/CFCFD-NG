# file diaphragm_test.py
#
# Author: David Gildfind, 28-Feb-2011
#
# This is a cut-down model of the X2 expansion tube as run for producing
# low-enthalpy, high total-pressure scramjet flow conditions.
# The simulation starts at the end of the piston stroke, with a compressed
# slug of He+Ar driver gas rupturing the primary diaphragm.
# There are 4 slugs of gas, separated by three diaphragms:
# 1 Primary diaphragm holds back driver gas.
#   Ruptures immediately, but opens in iris fashion over a fixed time.
# 2 Secondary diaphragm at downstream end of secondary (helium) driver.
#   Holds until there is sufficient delta_p across one of its block-to-block 
#   cell interfaces. Continues to hold for a period of time, 
#   then opens like iris but faster than the primary diaphragm.
# 3 Tertiary diaphragm holds the test gas in place until it is shock processed.
#   Holds like Diaphragm 2, except it opens fully after its hold period.
#
# Edited for User-Guide, 14-Jun-2011, PJ.

gdata.title = "Various diaphragm configuration test cases."
print gdata.title
gdata.dimensions = 2
gdata.axisymmetric_flag = 1

# General wall temp
T_start = 300.0

# Block 1 Compressed driver gas
p_init_1 = 30.0e6 # Initial fill pressure, [Pa].
T_init_1 = 2500.0 # Initial fill temperature, [T].
u_init_1 = 0.0 # Initial velocity, [m/s].
# Driver mass fractions are determined from partial pressures.
perc_He=80.0
perc_Ar=20.0
R_Ar=208.0
R_He=2077.0
R_m=(perc_He+perc_Ar)/(perc_He/R_He+perc_Ar/R_Ar)
mf_He=(perc_He*R_m)/(R_He*100.0)
mf_Ar=(perc_Ar*R_m)/(R_Ar*100.0)

# Block 2 Secondary (helium) driver.
p_init_2 = 100.0e3  	      # Initial fill pressure, [Pa].
T_init_2 = T_start            # Initial fill temperature, [T].
u_init_2 = 0.0                # Initial velocity, [m/s].

# Block 3 Test gas (air, LUT)
p_init_3 = 150.0e3  	      # Initial fill pressure, [Pa].
T_init_3 = T_start            # Initial fill temperature, [T].
u_init_3 = 0.0                # Initial velocity, [m/s].

# Block 4 Accelerator gas (air, LUT)
p_init_4 = 20.0  	      # Initial fill pressure, [Pa].
T_init_4 = T_start            # Initial fill temperature, [T].
u_init_4 = 0.0                # Initial velocity, [m/s].

# Gas compositions (mass fractions)
create_gas_file(model="ideal gas", species=['Ar', 'He', 'N2', 'air'], 
                fname="gas-model.lua", lut_file="cea-lut-air.lua.gz")
species_list = select_gas_model(fname="gas-model.lua")
print "species_list=", species_list

# species_list =       [  LUT_Air  Ar      He      N2    Perf_air] 
mfs_1 =                [  0.0,     mf_Ar,  mf_He,  0.0,  0.0     ]
mfs_2 =                [  0.0,     0.0,    1.0,    0.0,  0.0     ]
mfs_3 =                [  1.0,     0.0,    0.0,    0.0,  0.0     ]
mfs_4 =                [  1.0,     0.0,    0.0,    0.0,  0.0     ]

tube_1_init = FlowCondition(p=p_init_1, u=0.0, v=0.0, T=T_init_1, massf=mfs_1)
tube_2_init = FlowCondition(p=p_init_2, u=0.0, v=0.0, T=T_init_2, massf=mfs_2)
tube_3_init = FlowCondition(p=p_init_3, u=0.0, v=0.0, T=T_init_3, massf=mfs_3)
tube_4_init = FlowCondition(p=p_init_4, u=0.0, v=0.0, T=T_init_4, massf=mfs_4)


# Variables for geometry to change
d_tube = 0.085 # Diameter of compression tube (driver #1).

grid_scale_x=1.0
grid_scale_y=1.0

# number of cells per block - lengthwise
nnx1 = int(20*grid_scale_x+0.5)
nnx2 = int(60*grid_scale_x+0.5)

# number of cells per block - heightwise
nny = int(20*grid_scale_y+0.5)

# Nodes to define the segments of the tube. 
n1  = Node(0.0, 0.0, label="n1")
n2  = Node(0.0, d_tube/2, label="n2")
n3  = Node(0.1, 0.0, label="n3")
n4  = Node(0.1, d_tube/2, label="n4")
n5  = Node(0.4, 0.0, label="n5")
n6  = Node(0.4, d_tube/2, label="n6")
n7  = Node(0.7, 0.0, label="n7")
n8  = Node(0.7, d_tube/2, label="n8")
n9  = Node(1.0, 0.0, label="n9")
n10 = Node(1.0, d_tube/2, label="n10")

def make_box(SE, SW, NW, NE):
    """
    Make blocks from corner nodes.
              N
       NW___________NE
        |           |
        |           |
      W |           | E
        |___________|
       SW           SE   
              S
    """
    return make_patch(Line(NW, NE), Line(SE, NE), Line(SW,SE), Line(SW,NW))

pri_driver = Block2D(make_box(n3, n1, n2, n4), nni=nnx1, nnj=nny,
		     fill_condition=tube_1_init, label="pri_driver")
pri_driver.bc_list[WEST] = FixedTBC(T_start)
pri_driver.bc_list[NORTH] = FixedTBC(T_start)

sec_driver = Block2D(make_box(n5, n3, n4, n6), nni=nnx2, nnj=nny,
		     fill_condition=tube_2_init, label="sec_driver")
sec_driver.bc_list[NORTH] = FixedTBC(T_start)

shock_tube = Block2D(make_box(n7, n5, n6, n8), nni=nnx2, nnj=nny,
		     fill_condition=tube_3_init, label="shock_tube")
shock_tube.bc_list[NORTH] = FixedTBC(T_start)

accel_tube = Block2D(make_box(n9, n7, n8, n10), nni=nnx2, nnj=nny,
		     fill_condition=tube_4_init, label="accel_tube")
accel_tube.bc_list[EAST] = FixedTBC(T_start)
accel_tube.bc_list[NORTH] = FixedTBC(T_start)

identify_block_connections()
connect_blocks_2D(pri_driver, EAST, sec_driver, WEST, with_udf=1, 
		  filename="diaphragm_1.lua", is_wall=0, sets_conv_flux=0)
connect_blocks_2D(sec_driver, EAST, shock_tube, WEST, with_udf=1,
		  filename="diaphragm_2.lua", is_wall=0, sets_conv_flux=0)
connect_blocks_2D(shock_tube, EAST, accel_tube, WEST, with_udf=1,
		  filename="diaphragm_3.lua", is_wall=0, sets_conv_flux=0)

# Do a little more setting of global data.
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.shear_tolerance = 0.0
gdata.turbulence_flag = 0 # inactive turbulence model
gdata.turbulence_model = "baldwin_lomax" # not the default k_omega

gdata.max_time = 3.5e-4  # seconds
gdata.max_step = 5000000
gdata.dt = 1.0e-12
gdata.cfl = 0.25
gdata.dt_plot = 1.0e-6 # (seconds)
gdata.dt_history = 1.0e-6

# SVG Model Sketch Parameters
sketch.xaxis(0.0, 1.0, 0.2, -0.06)
sketch.yaxis(0.0, 0.2, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)
