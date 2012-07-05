# mms.py
# This file can be used to simulate the
# Method of Manufactured Solutions test case.
#
# Author: Rowan J. Gollan
# Updated: 05-Feb-2008
# Generalized to the viscous case by PJ, June 2011.
#

# Read some case parameters from a fixed file format.
fp = open('case.txt', 'r');
case_str = fp.readline().strip()
case = int(case_str)
flux_calc_str = fp.readline().strip()
flux_calc = fluxcalcIndexFromName[flux_calc_str]
x_order_str = fp.readline().strip()
x_order = int(x_order_str)
blocking = fp.readline().strip()
nn_str = fp.readline()
nn = int(nn_str)
fp.close()

gdata.title = "Method of Manufactured Solutions, Case=%d." % case

select_gas_model(fname='very-viscous-air.lua')
p0 = 1.0e5
T0 = p0 / 287.0  # rho0 = 1.0
if case == 1 or case == 3:
    # Supersonic inviscid flow
    u0 = 800.0; v0 = 800.0
    gdata.viscous_flag = 0
elif case == 2 or case == 4:
    # Subsonic viscous flow
    u0 = 70.0; v0 = 90.0
    gdata.viscous_flag = 1
else:
    print "UNKNOWN CASE"
    sys.exit()

initial = FlowCondition(p=p0, u=u0, v=v0, T=T0, massf=[1.0,])

a = Node(0.0, 0.0, label="a")
b = Node(1.0, 0.0, label="b")
c = Node(0.0, 1.0, label="c")
d = Node(1.0, 1.0, label="d")

ab1 = Node(0.2, 0.0, label="ab1")
ab2 = Node(0.5, -0.15, label="ab2")
ab3 = Node(0.85, 0.0, label="ab3")
BezAB = Bezier([a, ab1, ab2, ab3, b])

ac1 = Node(0.0, 0.1, label="ac1")
ac2 = Node(0.25, 0.2, label="ac2")
ac3 = Node(-0.3, 0.7, label="ac3")
ac4 = Node(0.0, 0.8, label="ac4")
BezAC = Bezier([a, ac1, ac2, ac3, ac4, c])

cd1 = Node(0.3, 1.0, label="cd1")
cd2 = Node(0.4, 1.3, label="cd2")
cd3 = Node(0.9, 1.0, label="cd3")
BezCD = Bezier([c, cd1, cd2, cd3, d])

bd1 = Node(1.0, 0.15, label="bd1")
bd2 = Node(0.75, 0.7, label="bd2")
bd3 = Node(1.0, 0.75, label="bd3")
BezBD = Bezier([b, bd1, bd2, bd3, d])

if case == 1 or case == 3:
    bc_list = [ExtrapolateOutBC(x_order=1), ExtrapolateOutBC(x_order=1),
               UserDefinedBC("udf-bc.lua"), UserDefinedBC("udf-bc.lua")]
elif case == 2 or case == 4:
    bc_list = [UserDefinedBC("udf-bc.lua"),]*4

if blocking == 'single':
    blk = Block2D(make_patch(BezCD, BezBD, BezAB, BezAC), 
                  nni=nn, nnj=nn,
                  bc_list=bc_list,
                  fill_condition=initial, label="blk")
elif blocking == 'multi':
    blk = SuperBlock2D(make_patch(BezCD, BezBD, BezAB, BezAC), 
                       nni=nn, nnj=nn, nbi=4, nbj=4,
                       bc_list=bc_list,
                       fill_condition=initial, label="blk")
else:
    print "UNKOWN BLOCKING SELECTION:", blocking
    sys.exit()

gdata.udf_file = "udf-source.lua"
gdata.udf_source_vector_flag = 1
gdata.flux_calc = flux_calc
gdata.x_order = x_order
if case == 1 or case == 3:
    gdata.max_time = 60.0e-3
    gdata.max_step = 1000000
    gdata.dt = 1.0e-6
    gdata.cfl = 0.5
elif case == 2 or case == 4:
    gdata.max_time = 150.0e-3
    gdata.max_step = 3000000
    gdata.dt = 1.0e-7
    gdata.cfl = 0.5
# For the verification tests,
# do NOT use the limiters
gdata.apply_limiter_flag = 0
gdata.extrema_clipping_flag = 0
gdata.stringent_cfl = 1
gdata.dt_plot = gdata.max_time/20.0


