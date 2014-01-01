##couette.py

gdata.dimensions = 2
select_gas_model(model='ideal gas', species=['air'])

x_max = 0.040
y_max = 0.010
nx = 20
ny = 10

def simple_rectangle(r, s, t=0.0):
    global x_max, y_max
    return (x_max*r, y_max*s, 0.0)

p_inf = 100.0e3  # Pa
u_max = 100.0    # m/s

initial = FlowCondition(p=p_inf, u=u_max, v=0.0)
def initial_flow(x, y, z):
    global y_max, T_inf, p_inf, u_max
    u = u_max * y / y_max
    return FlowCondition(p=p_inf, u=u, v=0.0, add_to_list=0).to_dict()

blk = Block2D(PyFunctionSurface(simple_rectangle), 
              nni=nx, nnj=ny,
              fill_condition=initial_flow,
              cf_list=4*[None,])

blk.set_BC("NORTH", "MOVING_WALL", r_omega=[0.0,0.0,0.0],v_trans=[u_max,0.0,0.0]) 
blk.bc_list[SOUTH] = AdiabaticBC()

# the WEST face is connected with the EAST face
connect_blocks_2D(blk,WEST,blk,EAST,check_corner_locations=False)

identify_block_connections()

gdata.title = "Couette flow (Just at start-up)"
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 50e-3
gdata.max_step = 20000
gdata.dt = 1.0e-9
gdata.dt_plot = 1e-3
gdata.dt_history = 1e-3
# The following scales provide a reasonable picture.
sketch.xaxis(0.0, x_max, 0.002, -0.05)
sketch.yaxis(0.0, y_max, 0.002, -0.04)
sketch.window(0.0, 0.0, x_max, y_max, 0.05, 0.05, 0.17, 0.17)

