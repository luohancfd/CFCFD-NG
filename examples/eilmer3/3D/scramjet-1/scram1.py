# A job description file for a simplified scramjet combustor and nozzle.
# PJ, Feb, Mar 2006
#     Apr 2007, updated for Elmer2
#     Feb 2010, updated for Eilmer3
# -------------------------------------------------------------------
from math import pi, sin, cos
def deg2rad(d):
    from math import pi
    return d/180.0*pi

# Define some geometric parameters that will be useful for specifying control points.
# See Katsu's sketch and workbook sketch on page 37 for labelling of points.
r1 = r2 = r3 = 34.0e-3
r4 = 41.5e-3
r5 = 14.6e-3
r6 = 18.0e-3
r7 = 28.0e-3
r8 = 30.0e-3
x1 = x8 = -95.981e-3
x2 = x7 = -65.981e-3
x3 = -60.564e-3
x4 = x5 = 9.019e-3
x6 = -15.858e-3
th1 = th8 = deg2rad(11.6)
th2 = th7 = deg2rad(14.0)
th3 = deg2rad(16.0)
th4 = th5 = deg2rad(29.0)
th6 = deg2rad(18.5)

# Create the collection of points for use in defining the surfaces.
p0  = Vector( 0.0, 0.0, 0.0  )
p1  = Vector( x1, r1*cos(th1), -r1*sin(th1) )
p2  = Vector( x2, r2*cos(th2), -r2*sin(th2) )
p3  = Vector( x3, r3*cos(th3), -r3*sin(th3) )
p4  = Vector( x4, r4*cos(th4), -r4*sin(th4) )
p5  = Vector( x5, r5*cos(th5), -r5*sin(th5) )
p6  = Vector( x6, r6*cos(th6), -r6*sin(th6) )
p7  = Vector( x7, r7*cos(th7), -r7*sin(th7) )
p8  = Vector( x8, r8*cos(th8), -r8*sin(th8) )
# Define the plane of symmetry
p9  = Vector( x1, r1, 0.0 )
p10 = Vector( x2, r2, 0.0 )
p11 = Vector( x3, r3, 0.0 )
p12 = Vector( x4, r4, 0.0 )
p13 = Vector( x5, r5, 0.0 )
p14 = Vector( x6, r6, 0.0 )
p15 = Vector( x7, r7, 0.0 )
p16 = Vector( x8, r8, 0.0 )
# A few more points along the x-axis for later generation of circular arcs.
p1_0 = Vector( x1, 0.0, 0.0 )
p2_0 = Vector( x2, 0.0, 0.0 )
p3_0 = Vector( x3, 0.0, 0.0 )
p4_0 = Vector( x4, 0.0, 0.0 )
p6_0 = Vector( x6, 0.0, 0.0 )

# North and south surfaces are defined directly as TrianglePatches
# In preparation, gather the control points into a single list.
# Note that a list is indexed from 0.
p = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
     p11, p12, p13, p14, p15, p16]
north = TrianglePatch(p, [1,8,7, 1,7,2, 2,7,3, 7,6,3, 3,6,4, 6,5,4],
                      [8,7,6,5], [1,2,3,4], [8,1], [5,4], "NORTH")
south = TrianglePatch(p, [9,16,15, 9,15,10, 10,15,11, 15,14,11, 11,14,12, 14,13,12],
                      [16,15,14,13], [9,10,11,12], [16,9], [13,12], "SOUTH")

# The top and bottom surfaces are somewhat curved.
# The surfaces in the physical model were cut on a lathe.
c9_1 = Arc(p9, p1, p1_0)
c10_2 = Arc(p10, p2, p2_0)
c9_10 = Line(p9, p10)
c1_2 = Line(p1, p2)
top = TrianglePatch(CoonsPatch(c9_10, c1_2, c9_1, c10_2), 1, 5, "COWL")
c11_3 = Arc(p11, p3, p3_0)
c10_11 = Line(p10, p11)
c2_3 = Line(p2, p3)
top.add(TrianglePatch(CoonsPatch(c10_11, c2_3, c10_2, c11_3), 1, 5))
c12_4 = Arc(p12, p4, p4_0)
c11_12 = Line(p11, p12)
c3_4 = Line(p3, p4)
top.add(TrianglePatch(CoonsPatch(c11_12, c3_4, c11_3, c12_4), 1, 5))

c16_8 = Arc(p16, p8, p1_0)
c15_7 = Arc(p15, p7, p2_0)
c16_15 = Line(p16, p15)
c8_7 = Line(p8, p7)
bottom = TrianglePatch(CoonsPatch(c16_15, c8_7, c16_8, c15_7), 1, 5, "CENTRE-BODY")
c14_6 = Arc(p14, p6, p6_0)
c15_14 = Line(p15, p14)
c7_6 = Line(p7, p6)
bottom.add(TrianglePatch(CoonsPatch(c15_14, c7_6, c15_7, c14_6), 1, 5))
c13_5 = Arc(p13, p5, p4_0)
c14_13 = Line(p14, p13)
c6_5 = Line(p6, p5)
bottom.add(TrianglePatch(CoonsPatch(c14_13, c6_5, c14_6, c13_5), 1, 5))

# The west and east faces are built to close the ends of the duct.
f_zero = LinearFunction(0.0, 0.0)
f_one = LinearFunction(0.0, 1.0)
f_linear = LinearFunction(1.0, 0.0)
cA = PathOnSurface(bottom, f_zero, f_linear)
cB = PathOnSurface(top, f_zero, f_linear)
cC = PathOnSurface(south, f_zero, f_linear)
cD = PathOnSurface(north, f_zero, f_linear)
west = CoonsPatch(cA, cB, cC, cD, "INLET")
cA = PathOnSurface(bottom, f_one, f_linear)
cB = PathOnSurface(top, f_one, f_linear)
cC = PathOnSurface(south, f_one, f_linear)
cD = PathOnSurface(north, f_one, f_linear)
east = CoonsPatch(cA, cB, cC, cD, "OUTLET")

# then assemble the surfaces into a volume.
pvolume = ParametricVolume([south, bottom, west, east, north, top],
                           "Simplified-scramjet")

# -------------------------------------------------------------------
# Now, to the flow part of the simulation definition...
gdata.title = "Simplified scramjet duct -- Katsu."
gdata.dimensions = 3
# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
gdata.max_time = 5.0e-3
gdata.dt = 1.0e-9
gdata.max_step = 1000
initialCond = FlowCondition(p=1000.0, u=0.0,    T=304.0)
inflowCond  = FlowCondition(p=50.0e3, u=2000.0, T=300.0)
nblocks = 3
blk = MultiBlock3D(label="duct",
                   parametric_volume=pvolume,
                   nbi=nblocks,
                   nni=40, nnj=20, nnk=20,
                   fill_condition=initialCond)
# Inlet to the nozzle is the first block.
blk.blks[0][0][0].set_BC("WEST", "SUP_IN", inflow_condition=inflowCond)
# Exit from the nozzle is in the last block.
blk.blks[nblocks-1][0][0].set_BC("EAST", "SUP_OUT")

# We are done with definitions; e3prep.py will do its work...

