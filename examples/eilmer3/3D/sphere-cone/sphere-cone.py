# A job description file for a spherically-blunted cone.
# PJ
# Elmer2 original: 13-Feb-2006
# Eilmer3 port: 06-Feb-2010

from math import *
from cfpylib.gasdyn.billig import x_from_y

# First, set the global data
gdata.title = "Sphere-cone."
gdata.dimensions = 3
select_gas_model(model='ideal gas', species=['air'])
gdata.max_time = 5.0e-3
gdata.dt = 1.0e-7
gdata.max_step = 1000

# Second, set up flow conditions
initialCond = FlowCondition(p=1000.0, u=0.0, T=300.0)
M_inf = 4.0
u_inf = M_inf * initialCond.flow.gas.a
inflowCond  = FlowCondition(p=50.0e3, u=u_inf, T=300.0)

# Third, set up the block

# The vehicle surface is defined as a path that is revolved about the x-axis.
Rnose = 1.0                  # radius of spherical nose
Angle = 45.0 * pi / 180.0    # angle of cone wrt x-axis 
Dmax = 4.0                   # diameter of base
# The conical section extends from nose to base radius. 
Length = (Dmax / 2.0 - Rnose * sin(Angle)) / sin(Angle)
c = Vector(0.0, 0.0, 0.0)    # centre of radius
a = Vector(-Rnose, 0.0, 0.0) # tip of nose
b = Vector(-Rnose*cos(Angle), Rnose*sin(Angle), 0.0) # join between sphere and cone
d = b + Length * Vector(cos(Angle), sin(Angle), 0.0) # skirt of cone
path = Polyline([Arc(a,b,c), Line(b,d)])
surf1 = RevolvedSurface(path, "vehicle_surface")
# To put a mesh onto this revolved surface, we define a query surface with
# a better outline for the block grid.
L2 = Dmax / 2.0 / sqrt(2.0)
# We have made sure that our query surface is within the bounds of the original.
q0 = Vector3(0.0, -L2, L2)
q1 = Vector3(0.0, -L2, 0.0)
q2 = Vector3(0.0,  L2, 0.0)
q3 = Vector3(0.0,  L2, L2)
qsurf2 = CoonsPatch(q0, q1, q2, q3, "query_surface")
east = MappedSurface(qsurf2, surf1)

# The outer mesh surface is derived from Billig's shock-shape correlation.
# In preparation for defining nodes, generate a few sample points
# along the expected shock position.
e = [] # will use a list to keep the nodes for the shock boundary
for y in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]:
    y *= Dmax/2.0  # scale up to cover the base of the vehicle
    # Note that we lie about the cone angle.  Detached shock.
    x = x_from_y(y, M_inf, theta=20.0/180.0*pi, axi=1, R_nose=Rnose)
    # print "x=", x, "y=", y
    # the outer boundary should be a little further than the shock itself
    e.append( Vector(-1.2*x, 1.2*y, 0.0) )
shock = Spline(e)
# print "shock=", shock
surf2 = RevolvedSurface(shock, "shock_surface")
L3 = e[-1].y / sqrt(2.0)
qs0 = Vector3(0.0, -L3, L3)
qs1 = Vector3(0.0, -L3, 0.0)
qs2 = Vector3(0.0,  L3, 0.0)
qs3 = Vector3(0.0,  L3, L3)
qsurf2 = CoonsPatch(qs0, qs1, qs2, qs3, "query_surface_shock")
west = MappedSurface(qsurf2, surf2)

p0 = west.eval(0.0, 0.0)
p1 = east.eval(0.0, 0.0)
p2 = east.eval(1.0, 0.0)
p3 = west.eval(1.0, 0.0)
p4 = west.eval(0.0, 1.0)
p5 = east.eval(0.0, 1.0)
p6 = east.eval(1.0, 1.0)
p7 = west.eval(1.0, 1.0)
# print "p0=", p0, "p1=", p1

# We shall assemble the other surfaces as CoonsPatch surfaces with
# their relevant bounding edges lying on the shock and body surfaces.
c76 = Line(p7, p6)
c32 = Line(p3, p2)
c37 = PathOnSurface(west, LinearFunction(0.0,1.0), LinearFunction(1.0,0.0))
c26 = PathOnSurface(east, LinearFunction(0.0,1.0), LinearFunction(1.0,0.0))
north = CoonsPatch(c32, c76, c37, c26, "symmetry-plane")

c45 = Line(p4, p5)
c01 = Line(p0, p1)
c04 = PathOnSurface(west, LinearFunction(0.0,0.0), LinearFunction(1.0,0.0))
c15 = PathOnSurface(east, LinearFunction(0.0,0.0), LinearFunction(1.0,0.0))
south = CoonsPatch(c01, c45, c04, c15, "south-outflow")

c47 = PathOnSurface(west, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))
c56 = PathOnSurface(east, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))
top = CoonsPatch(c45, c76, c47, c56)

c03 = PathOnSurface(west, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))
c12 = PathOnSurface(east, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))
bottom = CoonsPatch(c01, c32, c03, c12)

if 0:
    # Here is a bit of debug...
    # You can look to see that the surfaces are reasonable and join at the edges.
    # print surf2
    # print "surf2.eval(0.25,0.75)=", surf2.eval(0.25,0.75)
    # print "west=", west
    # print "west.eval(0.25,0.75)=", west.eval(0.25,0.75)
    print "Render to VRML"
    outfile = open("sphere-cone.wrl", "w")
    outfile.write("#VRML V2.0 utf8\n")
    outfile.write(east.vrml_str() + "\n")
    outfile.write(west.vrml_str() + "\n")
    outfile.write(north.vrml_str() + "\n")
    outfile.write(south.vrml_str() + "\n")
    outfile.write(top.vrml_str() + "\n")
    outfile.write(bottom.vrml_str() + "\n")
    outfile.close()
    import sys; sys.exit()
    
# Assemble the surfaces into a volume.
pvolume = ParametricVolume(north, east, south, west, top, bottom, "Sphere-cone")

blk = Block3D(label="first-block", nni=20, nnj=20, nnk=40,
              parametric_volume=pvolume,
              fill_condition=initialCond)
blk.set_BC("WEST", "SUP_IN", inflow_condition=inflowCond)
blk.set_BC("SOUTH", "SUP_OUT")
blk.set_BC("TOP", "SUP_OUT")
blk.set_BC("BOTTOM", "SUP_OUT")

