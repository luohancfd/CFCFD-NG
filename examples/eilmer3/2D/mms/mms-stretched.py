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

cf_list = [RobertsClusterFunction(1, 1, 1.5), RobertsClusterFunction(1, 0, 1.4),]*2

if case == 1 or case == 3:
    bc_list = [ExtrapolateOutBC(x_order=1), ExtrapolateOutBC(x_order=1),
               UserDefinedBC("udf-bc.lua"), UserDefinedBC("udf-bc.lua")]
elif case == 2 or case == 4:
    bc_list = [UserDefinedBC("udf-bc.lua"),]*4

if blocking == 'single':
    blk = Block2D(make_patch(Line(c,d), Line(b,d), Line(a,b), Line(a,c)), 
                  nni=nn, nnj=nn,
                  bc_list=bc_list,
                  cf_list=cf_list,
                  fill_condition=initial, label="blk")
elif blocking == 'multi':
    blk = SuperBlock2D(make_patch(Line(c,d), Line(b,d), Line(a,b), Line(a,c)), 
                       nni=nn, nnj=nn, nbi=4, nbj=4,
                       bc_list=bc_list,
                       cf_list=cf_list,
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
gdata.stringent_cfl = 1
gdata.dt_plot = gdata.max_time/20.0


