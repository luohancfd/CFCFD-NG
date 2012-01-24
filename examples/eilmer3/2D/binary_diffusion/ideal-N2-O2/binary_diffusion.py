# binary_diffusion.py
#
# A Python script to setup a simple
# binary diffusion test case.
#
# This Python file prepared by...
# Rowan J Gollan
# 20-Jun-2008
#
# Updated for eilmer3 on 17-Apr-2009

gdata.title = "Binary diffusion of N2 into O2 (and vice versa)"
select_gas_model(model='ideal gas', species=['N2', 'O2'])
gdata.diffusion_model = "FicksFirstLaw"
gdata.diffusion_flag = 1

left = FlowCondition( p=100.0e3, T=273.2, u=0.0, v=0.0, massf={'N2':1.0} )
right = FlowCondition( p=100.0e3, T=273.2, u=0.0, v=0.0, massf={'O2':1.0} )
#
# Fill function
def tube_gas(x, y, z, lfs=left, rfs=right):
    if x < 0.0:
        return lfs.to_dict()
    else:
        return rfs.to_dict()
#
# Geometry
xL = -2.0e-5
xR = 2.0e-5
ymin = 0.0
ymax = 0.1

a = Node(xL, ymin, label="a")
b = Node(xR, ymin, label="b")
c = Node(xL, ymax, label="c")
d = Node(xR, ymax, label="d")

south0 = Line(a,b)
west0 = Line(a,c)
north0 = Line(c,d)
east0 = Line(b,d)

nnx = 100
nny = 10

blk = Block2D(psurf=make_patch(north0, east0, south0, west0),
              nni=nnx, nnj=nny,
              bc_list=[AdiabaticBC()]*4,
              fill_condition=tube_gas,
              label="blk" )

gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 1.0e-6  # seconds
gdata.max_step = 200000
gdata.dt = 1.0e-10
gdata.dt_plot = gdata.max_time / 10.0
