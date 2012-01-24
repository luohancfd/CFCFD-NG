## \file couette.py
## \brief Laminar flow over a flat plate -- test derivatives
##
## This is set up something like incompressible Couette flow.
## With compressibility and dodge boundary conditions it will
## evolve into something else but we're only interested in the
## field derivatives at the first time step.
##
## BEWARE: this is for DEBUG ONLY
##
## \author PJ 04-May-2009 

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
T_inf = 1000.0   # degrees K

upperflow = FlowCondition(p=p_inf, u=u_max, v=0.0, T=T_inf)
def initial_flow(x, y, z):
    global y_max, T_inf, p_inf, u_max
    u = u_max * y / y_max # linear velocity profile
    T = T_inf * y / y_max
    return FlowCondition(p=p_inf, u=u, v=u, T=T, add_to_list=0).to_dict()

blk = Block2D(PyFunctionSurface(simple_rectangle), 
              nni=nx, nnj=ny,
              fill_condition=initial_flow,
              cf_list=4*[None,],
              bc_list=4*[UserDefinedBC("inflow.lua"),])

gdata.title = "Couette flow (Just at start-up)"
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_step = 1 
gdata.dt = 1.0e-9

# The following scales provide a reasonable picture.
sketch.xaxis(0.0, x_max, x_max/4, -0.05)
sketch.yaxis(0.0, y_max, y_max, -0.04)
sketch.window(0.0, 0.0, x_max, x_max, 0.05, 0.05, 0.17, 0.17)

