# plate_2d.py
# Ingo Jahn, Apr 2015

gdata.dimensions = 2
#select_gas_model(model='ideal gas', species=['air'])


R = 20e-3
Y = 20e-3
YY = 5.2e-3
nx0 = 5
nx1 = 50
nx2 = 50
nx3 = 5
ny = 50
ny4 = 10

p_inf = 35.8e3  # Pa
u_inf = 3257
v_inf = 0
T_inf = 268

#p_inf = 15e6  # Pa
#u_inf = 5e3
#v_inf = 0
#T_inf = 4000.0


# gas model selection and setting of inflow and initial flow conditions
gas_model_flag = 0
if gas_model_flag == 0:
    select_gas_model(model='ideal gas', species=['air'])
    inflow_condition = FlowCondition(p=p_inf, u=u_inf, v=v_inf, T=T_inf)
    initial = FlowCondition(p=p_inf, u=0.0, v=0.0, T=T_inf)
    print "Running Ideal gas simulation"
else:
    select_gas_model(model='thermally perfect gas', 
                 species=['N2','O2','N','O','NO'])
    # species index             0    1   2   3    4 
    f_inf = {"N2":0.78, "O2":0.22}  # nominal mass fractions for air
    set_reaction_scheme("gupta_etal_air_reactions.lua",reacting_flag=1)
    inflow_condition = FlowCondition(p=p_inf, u=u_inf, T=T_inf, massf=f_inf)
    initial = FlowCondition(p=p_inf, u=0.0, v=0.0, T=T_inf)
    M_inf = u_inf / inflow_condition.flow.gas.a
    print "Free-stream Mach number=", M_inf
    print "Running gupta_etal_air_reaction.lua"



rcf1 = RobertsClusterFunction(1,0,1.1) 
rcf2 = RobertsClusterFunction(0,1,1.15)

def simple_rectangle0(r, s, t=0.0):
    global R, Y
    return (R/10.*(r - 1.), Y*s, 0.0)

def simple_rectangle1(r, s, t=0.0):
    global R, Y
    return (R*r, Y*s, 0.0)

def simple_rectangle2(r, s, t=0.0):
    global R, Y, YY
    return (R + (R**2 - YY**2)**0.5 *r, YY*r + (Y-YY*r)*s, 0.0)

def simple_rectangle3(r, s, t=0.0):
    global R, Y, YY
    return (R + (R**2 - YY**2)**0.5 + R/10.*r, YY + (Y-YY)*s, 0.0)

def simple_rectangle4(r, s, t=0.0):
    global R, Y, YY
    return (R + (R**2 - YY**2)**0.5 + R/10.*r, YY + (Y-YY)*(s-1.)/10., 0.0)


blk0 = Block2D(PyFunctionSurface(simple_rectangle0), 
              nni=nx0, nnj=ny,
              fill_condition=initial,
              cf_list=[None,rcf1,None,rcf1])

blk1 = Block2D(PyFunctionSurface(simple_rectangle1), 
              nni=nx1, nnj=ny,
              fill_condition=initial,
              cf_list=[None,rcf1,None,rcf1])

blk2 = Block2D(PyFunctionSurface(simple_rectangle2), 
              nni=nx2, nnj=ny,
              fill_condition=initial,
              cf_list=[None,rcf1,None,rcf1])

blk3 = Block2D(PyFunctionSurface(simple_rectangle3), 
              nni=nx3, nnj=ny,
              fill_condition=initial,
              cf_list=[None,rcf1,None,rcf1])

blk4 = Block2D(PyFunctionSurface(simple_rectangle4), 
              nni=nx3, nnj=ny4,
              fill_condition=initial,
              cf_list=[None,rcf2,None,rcf2])

blk0.bc_list[NORTH] = ExtrapolateOutBC()
blk0.bc_list[WEST] = SupInBC(inflow_condition)
blk0.bc_list[SOUTH] = SlipWallBC()

blk1.bc_list[NORTH] = ExtrapolateOutBC()
blk1.bc_list[SOUTH] = FixedTBC(300.)
#blk1.bc_list[EAST] = ExtrapolateOutBC()

blk2.bc_list[NORTH] = ExtrapolateOutBC()
blk2.bc_list[SOUTH] = FixedTBC(300.)
#blk2.bc_list[EAST] = ExtrapolateOutBC()

blk3.bc_list[NORTH] = ExtrapolateOutBC()
blk3.bc_list[EAST] = ExtrapolateOutBC()

blk4.bc_list[WEST] = SlipWallBC()
blk4.bc_list[SOUTH] = SlipWallBC(300.)
blk4.bc_list[EAST] = ExtrapolateOutBC()

identify_block_connections()

gdata.title = "2-D moving flat plate"
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 100e-6
gdata.max_step = 10000
gdata.dt = 1.0e-9
gdata.dt_plot = 5e-6
gdata.dt_history = 1e-6
gdata.x_order = 2
gdata.interpolate_in_local_frame = True

# moving grid flag
gdata.dt_moving = 1e-6
gdata.moving_grid_flag = 1

# udf file to prescribe vertex velocity
gdata.udf_file = "udf-vtx.lua"
gdata.udf_vtx_velocity_flag = 1

# The following scales provide a reasonable picture.
sketch.xaxis(0.0, 0.1, 0.01, -0.005)
sketch.yaxis(0.0, 0.1, 0.01, -0.004)
sketch.window(0.0, 0.0, 0.1, 0.1, 0.05, 0.05, 0.15, 0.15)

