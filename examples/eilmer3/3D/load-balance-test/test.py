# Author: Rowan J. Gollan
# Date: 17-Aug-2012
# Place: University of Queensland, Brisbane, Australia
#
# This Python code is an input file for Eilmer3.
#
# This input file describes as simulation of a
# ducted flow with porthole injection. The
# main stream in nitrogen and the injected
# stream is hydrogen.
#

# --- 1: Set some global data
gdata.title = "Duct with porthole injection"
gdata.dimensions = 3

gpro_grid = 'blk.tmp'
gpro_conn = 'blk.tmp.conn'
gpro_pty = 'blk.tmp.pty'

duct_length = 50.0e-3 # m
T_wall = 293.0 # K
# Inflow stream
u_inf = 2250.0 # m/s
p_inf = 16.0e3 # Pa
T_inf = 650.0 # K
f_inf = {'N2':1.0, 'H2':0.0}

# Injected stream
u_inj = 1100.0 # m/s
rho_inj = 0.5325 # kg/m^3
T_inj = 255.0 # K
f_inj = {'N2':0.0, 'H2':1.0}

# Initially, set to 0 for inviscid flow
gdata.viscous_flag = 0

gdata.max_time = 5.0*duct_length/u_inf
gdata.dt = 1.0e-9
gdata.dt_plot = gdata.max_time/20

# --- 2: Set up flow conditions
select_gas_model(model='thermally perfect gas',
                 species=['N2', 'H2'])
# Compute pressure from temperature and density in the H2 stream
gmodel = get_gas_model_ptr()
Q = Gas_data(gmodel)
set_massf(Q, gmodel, f_inj)
Q.rho = rho_inj
Q.T[0] = T_inj
gmodel.eval_thermo_state_rhoT(Q)
p_inj = Q.p

inflow = FlowCondition(p=p_inf, T=[T_inf], u=u_inf, v=0.0, massf=f_inf)
initial = FlowCondition(p=p_inf/3.0, T=[T_inf], u=0.0, v=0.0, massf=f_inf)
injection = FlowCondition(p=p_inj, T=[T_inj], u=u_inj, v=0.0, massf=f_inj)

# --- 3: Set up blocks

# Read in grids from gridpro
grids = read_gridpro_grid(gpro_grid)

# Set up blocks
blk_list = []
for ib in range(len(grids)):
    blk_list.append(Block3D(grid=grids[ib], fill_condition=initial))

apply_gridpro_connectivity(gpro_conn, blk_list) 

bc_map = {'SUP_IN': inflow,
          'FIXED_T': T_wall}

apply_gridpro_bcs(gpro_pty, blk_list, bc_map)

# --- 4: Set numerics
gdata.flux_calc = AUSMDV
gdata.cfl = 0.4

