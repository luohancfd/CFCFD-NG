#
# This file can be used to simulate the
# Method of Manufactured Solutions test case.
#
# Author: Rowan J. Gollan
# Updated: 05-Feb-2008
#

gdata.title = "Method of Manufactured Solutions: Euler test case."
gdata.viscous_flag = 0
gdata.stringent_cfl = 1

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas',
                 species=['air'])

p0 = 1.0e5
u0 = 800.0
v0 = 800.0
T0 = p0 / 287.1

initial = FlowCondition(p=p0,  u=u0,    v=v0, T=T0,  massf=[1.0,])

a = Node(0.0, 0.0, label="a")
b = Node(1.0, 0.0, label="b")
c = Node(0.0, 1.0, label="c")
d = Node(1.0, 1.0, label="d")

ab = Line(a, b)
ac = Line(a, c)
cd = Line(c, d)
bd = Line(b, d)

nx = 16
ny = 16

blk_0 = Block2D(make_patch(cd, bd, ab, ac),
                nni=nx, nnj=ny,
                fill_condition=initial, label="blk-0")
blk_0.set_BC(NORTH, EXTRAPOLATE_OUT)
blk_0.set_BC(EAST, EXTRAPOLATE_OUT)
blk_0.set_BC(SOUTH, USER_DEFINED, filename="udf-bc.lua")
blk_0.set_BC(WEST, USER_DEFINED, filename="udf-bc.lua")

gdata.udf_file = "udf-source.lua"
gdata.udf_source_vector_flag = 1
gdata.flux_calc = AUSM
gdata.max_time = 20.0e-3
gdata.max_step = 2000
gdata.dt = 1.0e-6
gdata.fixed_time_step = False
gdata.cfl = 0.5
gdata.dt_plot = gdata.max_time/20.0


