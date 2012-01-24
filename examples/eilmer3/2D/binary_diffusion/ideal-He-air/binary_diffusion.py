# binary_diffusion.py
#
# A Python script to setup a simple
# binary diffusion test case.
#
# This file is an adaption of the N2-O2 test case.
#
# This Python file prepared by...
# Rowan J Gollan
# 20-Jun-2008

gdata.title = "Binary diffusion of He into air (and vice versa)"

#
# Gas model setup
#

select_gas_model(model='ideal gas',
                 species=['He', 'air'])

#
# Diffusion model
#

gdata.diffusion_model = "FicksFirstLaw"
gdata.diffusion_flag = 1

#
# Flow conditions
#

# mostly air, with a little bit of helium at left end
left = FlowCondition( p=100.0e3, T=300.0, u=0.0, v=0.0, massf=[0.1, 0.9] )
# air at right end
right = FlowCondition( p=100.0e3, T=300.0, u=0.0, v=0.0, massf=[0.0, 1.0] )

#
# Fill function
#

def tube_gas(x, y, z, lfs=left, rfs=right):
    if x < 0.0:
        return lfs.to_dict()
    else:
        return rfs.to_dict()

#
# Geometry
#

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
gdata.max_time = 5.0e-7  # seconds
gdata.max_step = 200000
gdata.dt = 1.0e-11
gdata.dt_plot = gdata.max_time / 10.0

