## \file couette.py
## \brief Laminar flow over a flat plate -- test derivatives
##
## This is set up something like incompressible Couette flow.
## With compressibility and dodge boundary conditions it will
## evolve into something else but we're only interested in the
## field derivatives at the first time step.
##
## Beware -- this is for DEBUG ONLY
##
## \author PJ 04-May-2009 

gdata.dimensions = 3
select_gas_model(model='ideal gas', species=['air'])

x_max = 0.040
y_max = 0.010
z_max = 0.010
nx = 20
ny = 10
nz = 10

def simple_hexagon(r, s, t):
    global x_max, y_max, z_max
    return (x_max*r, y_max*s, z_max*t)

p_inf = 100.0e3  # Pa
u_max = 100.0    # m/s
T_inf = 1000.0   # degrees K

upperflow = FlowCondition(p=p_inf, u=u_max, v=0.0, w=0.0, T=T_inf)
def initial_flow(x, y, z):
    global y_max, T_inf, p_inf, u_max
    u = u_max * y / y_max # linear velocity profile
    T = T_inf * y / y_max + 300.0
    return FlowCondition(p=p_inf, u=u, v=u, w=u, T=T_inf, add_to_list=0).to_dict()

blk = Block3D(PyFunctionVolume(simple_hexagon), 
              nni=nx, nnj=ny, nnk=nz,
              fill_condition=initial_flow,
              cf_list=12*[None,],
              bc_list=6*[UserDefinedBC("inflow.lua"),]
              )

gdata.title = "Couette flow (Just at start-up)"
gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_step = 1 
gdata.dt = 1.0e-9
