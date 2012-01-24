# binary_diffusion.py
#
# A Python script to setup a simple
# binary diffusion test case for 3D.
#
# This Python file prepared by...
# Rowan J Gollan 20-Jun-2008, 17-Apr-2009
# and PJ (3D) 09-Feb-2010
#
gdata.title = "Binary diffusion of N2 into O2 (and vice versa)"
gdata.dimensions = 3
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
H = 0.1
# from cfpylib.geom.box3d import makeSimpleBox
nnx = 100
nny = 10

blk = Block3D(parametric_volume=makeSimpleBox(xL, 0.0, 0.0, xR-xL, H, H),
              nni=nnx, nnj=nny, nnk=nny,
              bc_list=[AdiabaticBC()]*6,
              fill_condition=tube_gas,
              label="blk" )

gdata.viscous_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 1.0e-6  # seconds
gdata.max_step = 200000
gdata.dt = 1.0e-10
gdata.dt_plot = gdata.max_time / 10.0
