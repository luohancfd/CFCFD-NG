# back.py
# Conical nozzle from Back, Massier and Gier (1965)
gdata.title = "Flow through a conical nozzle."
print gdata.title

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
# The stagnation gas represents a reservoir condition.
stagnation_gas = FlowCondition(p=500.0e3, T=300.0)
low_pressure_gas = FlowCondition(p=30.0, T=300.0)

# Define geometry.
# The original paper specifies sizes in inches, Eilmer3 works in metres.
inch = 0.0254 # metres
L_subsonic = 3.0 * inch
L_nozzle = 3.0 * inch
R_tube = 1.5955 * inch
R_throat = 0.775 * inch
R_curve = 1.55 * inch # radius of curvature of throat profile
theta = 15.0 * math.pi / 180.0 # radians

# Compute the centres of curvature for the contraction profile.
height = R_throat + R_curve
hypot = R_tube + R_curve
base = math.sqrt(hypot*hypot - height*height)
centre_A = Node(0.0, height, label="centre_A")
centre_B = Node(-base, 0.0, label="centre_B")
fraction = R_tube/hypot
intersect_point = centre_B + Vector(fraction*base, fraction*height)

# The following Nodes will be rendered in the SVG file.
z0 = Node(-L_subsonic, 0.0)  # assemble from coordinates
p0 = Node(-L_subsonic, R_tube)
z1 = Node(centre_B)  # initialize from a previously defined Node
p1 = Node(centre_B + Vector(0.0,R_tube))  # vector sum
p2 = Node(intersect_point)
z2 = Node(p2.x, 0.0)  # on the axis, below p2
z3 = Node(0.0, 0.0)
p3 = Node(0.0, R_throat)
# Compute the details of the conical nozzle
p4 = Node(R_curve*math.sin(theta), height - R_curve*math.cos(theta))
z4 = Node(p4.x, 0.0)
L_cone = L_nozzle - p4.x
p5 = Node(p4 + Vector(L_cone, L_cone*math.tan(theta)))
z5 = Node(p5.x, 0.0)

north0 = Polyline([Line(p0,p1),Arc(p1,p2,centre_B),Arc(p2,p3,centre_A)])
east0west1 = Line(z3, p3)
south0 = Line(z0, z3)
west0 = Line(z0, p0)
north1 = Polyline([Arc(p3,p4,centre_A), Line(p4,p5)])
east1 = Line(z5, p5)
south1 = Line(z3, z5)

# Define the blocks, boundary conditions and set the discretisation.
nx0 = 50; nx1 = 60; ny = 30
subsonic_region = Block2D(make_patch(north0, east0west1, south0, west0),
                          nni=nx0, nnj=ny,
                          fill_condition=stagnation_gas,
                          label="subsonic-region")
supersonic_region = Block2D(make_patch(north1, east1, south1, east0west1),
                            nni=nx1, nnj=ny,
                            fill_condition=low_pressure_gas,
                            label="supersonic-region")
identify_block_connections()
subsonic_region.bc_list[WEST] = SubsonicInBC(stagnation_gas)
supersonic_region.bc_list[EAST] = ExtrapolateOutBC()

# Flow-history to be recorded at the following points.
HistoryLocation(0.001, 0.002, label="nozzle-throat")
HistoryLocation(L_nozzle-0.001, 0.002, label="nozzle-exit")

# Do a little more setting of global data.
gdata.axisymmetric_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 4.0e-3  # seconds
gdata.max_step = 50000
gdata.dt = 1.0e-7
gdata.dt_plot = 0.2e-3
gdata.dt_history = 10.0e-6

sketch.xaxis(-0.10, 0.08, 0.05, -0.01)
sketch.yaxis(  0.0, 0.05, 0.02, -0.015)
sketch.window(-0.10, 0.0, 0.10, 0.05, 0.05, 0.05, 0.25, 0.10)
