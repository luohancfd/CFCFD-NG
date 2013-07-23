# jones-stator.py
# Import Jones' turbine grid and set up a bit of the Eilmer3 sim.
# PJ, 26-Feb-2010 adapted from the nasaRotor file.

import math

gdata.title = "Nasa Radial Turbine"
print "\n%s\n" % gdata.title
gdata.dimensions = 3
gdata.dt = 1.0e-8
gdata.flux_calc = ADAPTIVE
gdata.t_order = 2
gdata.x_order = 1
gdata.cfl=0.3
gdata.viscous_flag = 0
gdata.turbulence_model = 'k-omega'
gdata.max_time = 0.1
gdata.max_step = 2000
gdata.dt_plot = 0.01

initialFlow = 1

def stator_flow_condition(x, y, z):
    """
    The initial flow conditions into the outer part of the stator.
    Returns a dictionary of flow properties.
    """
    # The following conditions are for the nasaRotor and
    # need to be updated for the Jones turbine conditions.
    initMach = 0.3
    speedSound = 566.0
    initFlowAngle = 76.0
    initPres = 400000.0
    velMag = initMach * speedSound
    r = sqrt( x * x + y * y)
    alpha = initFlowAngle * 2.0 * math.pi / 360.0
    vel_t = velMag * math.sin(alpha)
    vel_r = -velMag * math.cos(alpha)
    u_x = (-y * vel_t + x * vel_r) / r
    v_y = (x * vel_t + y * vel_r) / r
    return FlowCondition (p=initPres, u = u_x, v=v_y, T=800.0, add_to_list=0).to_dict()


def getLUAfunction(bcType, row):
    if (bcType == "inlet"):
        lua = "udf-subsonic-in-bc.lua"
    elif (bcType == "outlet"):
        lua = "udf-subsonic-out-bc.lua"
    elif (bcType == "periodic"):
        lua = "udf-stator-periodic-bc.lua"

    return lua
    


select_gas_model(model='ideal gas', species=['air'])

# Flow domain is defined by Carlos' ICEM grid written in CGNS format.
from cgns_grid import read_ICEM_CGNS_grids

cgns_data = read_ICEM_CGNS_grids('Stator.cgns', labelStem='stator', gridScale=0.001)
nb = cgns_data['nblock']
print "Stator blocks read:", nb
stator_blks = []

for ib in range(nb):
    blk = {}
    if (initialFlow == 1):
        blk = Block3D(grid=cgns_data['grids'][ib], fill_condition=stator_flow_condition) 
    else:
        stator_fc = ExistingSolution("jones-stator", '.', nb, 9999, assume_same_grid=1)
        blk = Block3D(grid=cgns_data['grids'][ib], fill_condition=stator_fc) 

    stator_blks.append(blk)

identify_block_connections()

# Apply boundary conditions

for bc in cgns_data['bcs']:
    luaFunction = getLUAfunction(bc['type'], 'stator')   
    face = bc['face']
    block = bc['block']
    stator_blks[block].set_BC(face, USER_DEFINED, filename=luaFunction)



