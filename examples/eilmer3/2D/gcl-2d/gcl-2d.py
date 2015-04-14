# gcl-2d.py
# Jason (Kan) Qin, Jan 2015

gdata.dimensions = 2
select_gas_model(model='ideal gas', species=['air'])

x_max = 0.1
y_max = 0.1
nx = 10
ny = 10

def simple_rectangle(r, s, t=0.0):
    global x_max, y_max
    return (x_max*r, y_max*s, 0.0)

p_inf = 100.0e3  # Pa
u_inf = 1.0

initial = FlowCondition(p=p_inf, u=u_inf, v=0.0)

blk = Block2D(PyFunctionSurface(simple_rectangle), 
              nni=nx, nnj=ny,
              fill_condition=initial,
              cf_list=4*[None,])

blk.bc_list[NORTH] = SlipWallBC()
blk.bc_list[EAST] = SlipWallBC()
blk.bc_list[SOUTH] = SlipWallBC()
blk.bc_list[WEST] = SlipWallBC()

identify_block_connections()

gdata.title = "random grid motion"
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 20e-3
gdata.max_step = 10000
gdata.dt = 1.0e-9
gdata.dt_plot = 1e-3
gdata.dt_history = 1e-3
gdata.gasdynamic_update_scheme = "pc"
print "The gas-dynamic update scheme is", gdata.gasdynamic_update_scheme
gdata.x_order = 2
gdata.interpolate_in_local_frame = True

# moving grid flag
gdata.dt_moving = 1e-3
gdata.moving_grid_flag = 1

# udf file to prescribe vertex velocity
gdata.udf_file = "udf-vtx.lua"
gdata.udf_vtx_velocity_flag = 1

# The following scales provide a reasonable picture.
sketch.xaxis(0.0, 0.040, 0.01, -0.005)
sketch.yaxis(0.0, 0.010, 0.01, -0.004)
sketch.window(0.0, 0.0, 0.040, 0.010, 0.05, 0.05, 0.15, 0.075)

