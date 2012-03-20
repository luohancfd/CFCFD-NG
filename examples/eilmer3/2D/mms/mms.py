# mms.py
# This file can be used to simulate the
# Method of Manufactured Solutions test case.
#
# Author: Rowan J. Gollan
# Updated: 05-Feb-2008
# Generalized to the viscous case by PJ, June 2011.
#
fp = open('case.txt', 'r'); case_str = fp.read().strip(); fp.close()
case = int(case_str)
gdata.title = "Method of Manufactured Solutions, Case=%d." % case

select_gas_model(fname='very-viscous-air.lua')
p0 = 1.0e5
T0 = p0 / 287.1  # rho0 = 1.0
if case == 1:
    # Supersonic inviscid flow
    u0 = 800.0; v0 = 800.0
    gdata.viscous_flag = 0
elif case == 2:
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

blk = SuperBlock2D(make_patch(Line(c,d), Line(b,d), Line(a,b), Line(a,c)), 
                   nni=16, nnj=16, nbi=2, nbj=2,
                   bc_list=[UserDefinedBC("udf-bc.lua"),]*4,
                   fill_condition=initial, label="blk")

gdata.udf_file = "udf-source.lua"
gdata.udf_source_vector_flag = 1
gdata.flux_calc = AUSM
if case == 1:
    gdata.max_time = 20.0e-3
    gdata.max_step = 2000
    gdata.dt = 1.0e-6
    gdata.cfl = 0.5
elif case == 2:
    gdata.max_time = 80.0e-3
    gdata.max_step = 20000
    gdata.dt = 1.0e-6
    gdata.cfl = 0.5
gdata.stringent_cfl = 1
gdata.dt_plot = gdata.max_time/20.0


