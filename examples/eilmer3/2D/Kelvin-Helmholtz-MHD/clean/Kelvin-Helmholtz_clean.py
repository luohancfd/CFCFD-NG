## \file Kelvin-Helmholtz_clean.py
## \brief Simple setup file for testing MHD implementation in Eilmer3.
## \author DB, 17-Nov-2014

from math import sin, sqrt, pi, tanh, exp

job_title = "Kelvin-Helmholtz test"
print job_title

# We can set individual attributes of the global data object.
gdata.dimensions = 2
gdata.title = job_title
gdata.mhd_flag = 1
gdata.axisymmetric_flag = 0
gdata.stringent_cfl = 1  # to match the old mb_cns behaviour

gdata.flux_calc = HLLE

#==============================================================================
# FLOW DEFINITION
#==============================================================================

# Accept defaults for He
select_gas_model(model='ideal gas', species=['air'])

rho0 = 1.2
p0 = 101.325e3
B0 = sqrt(p0)
R = 287.0
gam = 7.0/5.0

L_inf = 1.0
u_inf = sqrt(p0/rho0)
t_inf = L_inf/u_inf

gdata.divB_damping_length = 1.0
gdata.div_clean_flag = 1

rho = 1.0*rho0
p = 50.0*p0
T = p/(rho*R)

flow = FlowCondition(p=p, T=T)

def UDF_flow(x,y,z):
    """
    user defined flow condition
    """
    
    out = flow.to_dict()
    
    out['vel.x'] = 5.0*(tanh(20*(y+0.5))-(tanh(20*(y-0.5))+1))*u_inf
    out['vel.y'] = 0.25*sin(2*pi*x)*(exp(-100*(y+0.5)**2)-exp(-100*(y-0.5)**2))*u_inf
    out['vel.z'] = 0.0
    
    out['B.x'] = 1.0*B0
    out['B.y'] = 0.0
    out['B.z'] = 0.0
    
    return out
    
#==============================================================================
# RUN DEFINITION
#==============================================================================

gdata.max_time = 0.5*t_inf
gdata.max_step = 10000
gdata.dt = 1.0e-9
gdata.dt_plot = gdata.max_time/100.0

#==============================================================================
# DOMAIN DEFINITION
#==============================================================================

start_x = 0.0
mid_x = 0.5
stop_x = 1.0

start_y = -1.0
mid_y = 0.0
stop_y = 1.0

# nodes
a = Node(start_x, start_y, 0.0, label="a")
b = Node(mid_x, start_y, 0.0, label="b")
c = Node(stop_x, start_y, 0.0, label="c")

d = Node(start_x, mid_y, 0.0, label="d")
e = Node(mid_x, mid_y, 0.0, label="e")
f = Node(stop_x, mid_y, 0.0, label="f")

g = Node(start_x, stop_y, 0.0, label="g")
h = Node(mid_x, stop_y, 0.0, label="h")
i = Node(stop_x, stop_y, 0.0, label="i")

# lines

ab = Line(a, b)
bc = Line(b, c)
de = Line(d, e)
ef = Line(e, f)
gh = Line(g, h)
hi = Line(h, i)

ad = Line(a, d)
dg = Line(d, g)
be = Line(b, e)
eh = Line(e, h)
cf = Line(c, f)
fi = Line(f, i)

# patches

p0 = make_patch(de, be, ab, ad)
p1 = make_patch(ef, cf, bc, be)
p2 = make_patch(hi, fi, ef, eh)
p3 = make_patch(gh, eh, de, dg)

# blocks


# Define the blocks, boundary conditions and set the discretisation.

N = 64
blk_0 = Block2D(p0, nni=N/2, nnj=N,
                fill_condition=UDF_flow, label="flow")
                
blk_1 = Block2D(p1, nni=N/2, nnj=N,
                fill_condition=UDF_flow, label="flow")
                
blk_2 = Block2D(p2, nni=N/2, nnj=N,
                fill_condition=UDF_flow, label="flow")
                
blk_3 = Block2D(p3, nni=N/2, nnj=N,
                fill_condition=UDF_flow, label="flow")
                
                
identify_block_connections()

connect_blocks_2D(blk_0, SOUTH, blk_3, NORTH, check_corner_locations=False)
connect_blocks_2D(blk_1, SOUTH, blk_2, NORTH, check_corner_locations=False)
connect_blocks_2D(blk_0, WEST, blk_1, EAST, check_corner_locations=False)
connect_blocks_2D(blk_3, WEST, blk_2, EAST, check_corner_locations=False)
                     
                     
#==============================================================================
# SVG
#==============================================================================

sketch.xaxis(start_x, stop_x, 0.5, -0.05)
sketch.yaxis(start_y, stop_y, 0.5, -0.05)
sketch.window(start_x, start_y, stop_x, stop_y, 0.02, 0.02, 0.18, 0.18)
