# sc10_3D.py
# 3D Compressor Blade Standard Condition 10
#
# Hannes Wojciak, Paul Petrie-Repar
#    February 2008: Original implementation
# Peter J.
#    02-Oct-2008: Port to Elmer3 and add boundary conditions.
#
# from cfpylib.geom.path import Polyline2, Spline2
from math import sqrt

# ------------- global data ---------------------
gdata.title = "inviscid Euler for 3D-sc10"
gdata.dimensions = 3
# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
gdata.viscous_flag = 0 # inviscid simulation
gdata.max_time = 0.300
gdata.max_step = 100000
gdata.dt_plot = 0.020
gdata.dt = 1.0e-6

# -----------flow conditions -------------------
p_tot = 100.0e3 # Pa
T_tot = 300.0 # degree K
gma = 1.4
Rgas = 287.0  # J/kg.K
a_tot = sqrt(gma*Rgas*T_tot)
M_exit = 0.465
T0_T = 1 + (gma-1.0)/2.0 * M_exit * M_exit
p0_p = T0_T**(gma/(gma-1.0))
print "p0_p=", p0_p, "T0_T=", T0_T
p_exit = p_tot / p0_p
T_exit = T_tot / T0_T
u_exit = M_exit * a_tot / sqrt(T0_T)
print "p_exit=", p_exit, "T_exit=", T_exit, "u_exit=", u_exit
initialCond = FlowCondition(p=p_exit, u=u_exit, T=T_exit)
# inflowCond = FlowCondition(p=p_tot, T=T_tot) # inflow condition set within UDF BC

# ----------- constants and definitions ------------
R_h = 3.395   # radius of hub
R_m = 3.8195  # radius of midspan
R_s = 4.244   # radius of shroud

def toRadians(degrees):
    import math
    return degrees * math.pi / 180.0

def HubCoordinates(MidspanNode, R_h=R_h, R_m=R_m, R_s=R_s):
    from math import sin, cos
    angle = MidspanNode.y / (-R_m) # circumferential angle of Node in radians
    HubNode = Node(0.0,0.0,0.0)
    HubNode.x = MidspanNode.x
    HubNode.y = -R_h * sin(angle)
    HubNode.z = R_h * cos(angle)
    return HubNode

def ShroudCoordinates(MidspanNode, R_h=R_h, R_m=R_m, R_s=R_s):
    from math import sin, cos
    angle = MidspanNode.y / (-R_m) # circumferential angle of Node in radians
    ShroudNode = Node(0.0,0.0,0.0)
    ShroudNode.x = MidspanNode.x
    ShroudNode.y = -R_s * sin(angle)
    ShroudNode.z = R_s * cos(angle)
    return ShroudNode

def AxisCoordinates(MidspanNode):
    AxisNode = Node(0.0,0.0,0.0)
    AxisNode.x = MidspanNode.x
    return AxisNode

# -------------------- some 2D-node definitions ---------------------
LE_up = Node(0.007375, 0.038160, 0.0, label="LE_up")
LE_out_up = Node(-0.1, 0.07, 0.0, label="LE_out_up")
LE = Node(0.0, 0.0, 0.0, label="LE")
LE_out = Node(-0.05, -0.09, 0.0, label="LE_out")
LE_down = Node(0.026541, 0.015230, 0.0, label="LE_down")
LE_out_down = Node(0.09, -0.07, 0.0, label="LE_out_down")
TE = Node(0.707107, 0.707107, 0.0, label="TE")
TE_up = Node(0.7, 0.85, 0.0, label="TE_up")
TE_down = Node(0.75, 0.6, 0.0, label="TE_down")

SS = Node(0.18,0.44,0.0) # Node for spline_SS
Front_up = Node(-0.1,-0.01,0.0) # Node for spline_front_up
Front_down = Node(0.02,-0.1,0.0) # Node for spline_front_down
PS = Node(0.45,0.34,0.0) # Node for spline_PS

# ------------ hub and shroud cylindrical surface definitions ------------
A = Node(-1.0,-0.5,0.0)
B = Node(1.7,-0.5,0.0)
C = Node(1.7,1.2,0.0)
D = Node(-1.0,1.2,0.0)

path_s = Line(HubCoordinates(A),HubCoordinates(B))
path_n = Line(HubCoordinates(D),HubCoordinates(C))
path_w = Arc(HubCoordinates(A),HubCoordinates(D),AxisCoordinates(A))
path_e = Arc(HubCoordinates(B),HubCoordinates(C),AxisCoordinates(B))
HubCylinder = CoonsPatch(path_s, path_n, path_w, path_e)

path_s = Line(ShroudCoordinates(A),ShroudCoordinates(B))
path_n = Line(ShroudCoordinates(D),ShroudCoordinates(C))
path_w = Arc(ShroudCoordinates(A),ShroudCoordinates(D),AxisCoordinates(A))
path_e = Arc(ShroudCoordinates(B),ShroudCoordinates(C),AxisCoordinates(B))
ShroudCylinder = CoonsPatch(path_s, path_n, path_w, path_e)

# ------------------------ mapped surfaces --------------------------
# query surface = 2D-patch along profile
path_s = Spline2("sc10_inner1.dat")
path_n = Spline([LE_out_up,SS,TE_up])
path_w = Line(LE_up,LE_out_up)
path_e = Line(TE,TE_up)
qsurface = CoonsPatch(path_s, path_n, path_w, path_e)
inner1_h = MappedSurface(qsurface,HubCylinder)
RenderList.append(inner1_h)
inner1_s = MappedSurface(qsurface,ShroudCylinder)

path_s = Line(LE_out,LE)
path_n = Line(LE_out_up,LE_up)
path_w = Spline([LE_out,Front_up,LE_out_up])
path_e = Spline2("sc10_inner2.dat")
qsurface = CoonsPatch(path_s, path_n, path_w, path_e)
inner2_h = MappedSurface(qsurface,HubCylinder)
RenderList.append(inner2_h)
inner2_s = MappedSurface(qsurface,ShroudCylinder)

path_s = Spline([LE_out,Front_down,LE_out_down])
path_n = Spline2("sc10_inner3.dat")
path_n.reverse()
path_w = Line(LE_out,LE)
path_e = Line(LE_out_down,LE_down)
qsurface = CoonsPatch(path_s, path_n, path_w, path_e)
inner3_h = MappedSurface(qsurface,HubCylinder)
RenderList.append(inner3_h)
inner3_s = MappedSurface(qsurface,ShroudCylinder)

path_s = Spline([LE_out_down,PS,TE_down])
path_n = Spline2("sc10_inner4.dat")
path_n.reverse()
path_w = Line(LE_out_down,LE_down)
path_e = Line(TE_down,TE)
qsurface = CoonsPatch(path_s, path_n, path_w, path_e)
inner4_h = MappedSurface(qsurface,HubCylinder)
RenderList.append(inner4_h)
inner4_s = MappedSurface(qsurface,ShroudCylinder)

# ------------------ node definitons around profile ---------------------
TE_h = inner1_h.eval(1.0,0.0)
TE_s = inner1_s.eval(1.0,0.0)
TE_up_h = inner1_h.eval(1.0,1.0)
TE_up_s = inner1_s.eval(1.0,1.0)
TE_down_h = inner4_h.eval(1.0,0.0)
TE_down_s = inner4_s.eval(1.0,0.0)

LE_up_h = inner1_h.eval(0.0,0.0)
LE_up_s = inner1_s.eval(0.0,0.0)
LE_out_up_h = inner1_h.eval(0.0,1.0)
LE_out_up_s = inner1_s.eval(0.0,1.0)

LE_h = inner2_h.eval(1.0,0.0)
LE_s = inner2_s.eval(1.0,0.0)
LE_out_h = inner2_h.eval(0.0,0.0)
LE_out_s = inner2_s.eval(0.0,0.0)

LE_down_h = inner3_h.eval(1.0,1.0)
LE_down_s = inner3_s.eval(1.0,1.0)
LE_out_down_h = inner3_h.eval(1.0,0.0)
LE_out_down_s = inner3_s.eval(1.0,0.0)

# --------------- profile intersections with hub & shroud ----------------
profile_SS_h = PathOnSurface(inner1_h, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))
profile_SS_s = PathOnSurface(inner1_s, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))
profile_Front_up_h = PathOnSurface(inner2_h, LinearFunction(0.0,1.0), LinearFunction(1.0,0.0))
profile_Front_up_s = PathOnSurface(inner2_s, LinearFunction(0.0,1.0), LinearFunction(1.0,0.0))
profile_Front_down_h = PathOnSurface(inner3_h, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))
profile_Front_down_s = PathOnSurface(inner3_s, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))
profile_PS_h = PathOnSurface(inner4_h, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))
profile_PS_s = PathOnSurface(inner4_s, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))

# --------------- definition of outer splines for c-mesh -------------------

spline_SS_h = PathOnSurface(inner1_h, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))
spline_SS_s = PathOnSurface(inner1_s, LinearFunction(1.0,0.0), LinearFunction(0.0,1.0))
spline_Front_up_h = PathOnSurface(inner2_h, LinearFunction(0.0,0.0), LinearFunction(1.0,0.0))
spline_Front_up_s = PathOnSurface(inner2_s, LinearFunction(0.0,0.0), LinearFunction(1.0,0.0))
spline_Front_down_h = PathOnSurface(inner3_h, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))
spline_Front_down_s = PathOnSurface(inner3_s, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))
spline_PS_h = PathOnSurface(inner4_h, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))
spline_PS_s = PathOnSurface(inner4_s, LinearFunction(1.0,0.0), LinearFunction(0.0,0.0))

# ---------------------- inner1 --------------------------
A_h = LE_up_h
A_s = LE_up_s
B_h = TE_h
B_s = TE_s
C_h = TE_up_h
C_s = TE_up_s
D_h = LE_out_up_h
D_s = LE_out_up_s

p01 = profile_SS_h
p12 = Helix(B_h,C_h,AxisCoordinates(B_h),AxisCoordinates(C_h))
p32 = spline_SS_h
p03 = Helix(A_h,D_h,AxisCoordinates(A_h),AxisCoordinates(D_h))
p45 = profile_SS_s
p56 = Helix(B_s,C_s,AxisCoordinates(B_s),AxisCoordinates(C_s))
p76 = spline_SS_s
p47 = Helix(A_s,D_s,AxisCoordinates(A_s),AxisCoordinates(D_s))
p04 = Line(A_h,A_s)
p15 = Line(B_h,B_s)
p26 = Line(C_h,C_s)
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
RenderList.append(pvolume)
cflist = [None,]*12; 
inner1 = Block3D(label="inner1", nni=50, nnj=6, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)

# ---------------------- inner2 --------------------------
A_h = LE_out_h
A_s = LE_out_s
B_h = LE_h
B_s = LE_s
C_h = LE_up_h
C_s = LE_up_s
D_h = LE_out_up_h
D_s = LE_out_up_s

p01 = Helix(A_h,B_h,AxisCoordinates(A_h),AxisCoordinates(B_h))
p12 = profile_Front_up_h
p32 = Helix(D_h,C_h,AxisCoordinates(D_h),AxisCoordinates(C_h))
p03 = spline_Front_up_h
p45 = Helix(A_s,B_s,AxisCoordinates(A_s),AxisCoordinates(B_s))
p56 = profile_Front_up_s
p76 = Helix(D_s,C_s,AxisCoordinates(D_s),AxisCoordinates(C_s))
p47 = spline_Front_up_s
p04 = Line(A_h,A_s)
p15 = Line(B_h,B_s)
p26 = Line(C_h,C_s)
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cluster_k = RobertsClusterFunction(1, 0, 1.05)
cflist = [None,] + [cluster_k,] +[None,]*3 +[cluster_k,]+ [None,]*6;
inner2 = Block3D(label="inner2", nni=6, nnj=8, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)

# ---------------------- inner3 --------------------------
A_h = LE_out_h
A_s = LE_out_s
B_h = LE_out_down_h
B_s = LE_out_down_s
C_h = LE_down_h
C_s = LE_down_s
D_h = LE_h
D_s = LE_s

p01 = spline_Front_down_h
p12 = Helix(B_h,C_h,AxisCoordinates(B_h),AxisCoordinates(C_h))
p32 = profile_Front_down_h
p03 = Helix(A_h,D_h,AxisCoordinates(A_h),AxisCoordinates(D_h))
p45 = spline_Front_down_s
p56 = Helix(B_s,C_s,AxisCoordinates(B_s),AxisCoordinates(C_s))
p76 = profile_Front_down_s
p47 = Helix(A_s,D_s,AxisCoordinates(A_s),AxisCoordinates(D_s))
p04 = Line(A_h,A_s)
p15 = Line(B_h,B_s)
p26 = Line(C_h,C_s)
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cluster_k = RobertsClusterFunction(0, 1, 1.05)
cflist = [None,]*2 + [cluster_k,] +[None,]*3 +[cluster_k,]+ [None,]*5;
inner3 = Block3D(label="inner3", nni=8, nnj=6, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)

# -------------------inner4 ------------
A_h = LE_out_down_h
A_s = LE_out_down_s
B_h = TE_down_h
B_s = TE_down_s
C_h = TE_h
C_s = TE_s
D_h = LE_down_h
D_s = LE_down_s

p01 = spline_PS_h
p12 = Helix(B_h,C_h,AxisCoordinates(B_h),AxisCoordinates(C_h))
p32 = profile_PS_h
p03 = Helix(A_h,D_h,AxisCoordinates(A_h),AxisCoordinates(D_h))
p45 = spline_PS_s
p56 = Helix(B_s,C_s,AxisCoordinates(B_s),AxisCoordinates(C_s))
p76 = profile_PS_s
p47 = Helix(A_s,D_s,AxisCoordinates(A_s),AxisCoordinates(D_s))
p04 = Line(A_h,A_s)
p15 = Line(B_h,B_s)
p26 = Line(C_h,C_s)
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
inner4 = Block3D(label="inner4", nni=50, nnj=6, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)

# -------------------- inflow1 -------------------
A = Node(-1.0, 0.15,0.0)
B_h = LE_out_up_h
B_s = LE_out_up_s
C = Node(-0.3,0.5,0.0)
D = Node(-1.0,0.5,0.0)

p01 = Helix(HubCoordinates(A),B_h,AxisCoordinates(A),AxisCoordinates(B_h))
p12 = Helix(B_h,HubCoordinates(C),AxisCoordinates(B_h),AxisCoordinates(C))
p32 = Line(HubCoordinates(D),HubCoordinates(C))
p03 = Arc(HubCoordinates(A),HubCoordinates(D),AxisCoordinates(A))
p45 = Helix(ShroudCoordinates(A),B_s,AxisCoordinates(A),AxisCoordinates(B_s))
p56 = Helix(B_s,ShroudCoordinates(C),AxisCoordinates(B_s),AxisCoordinates(C))
p76 = Line(ShroudCoordinates(D),ShroudCoordinates(C))
p47 = Arc(ShroudCoordinates(A),ShroudCoordinates(D),AxisCoordinates(A))
p04 = Line(HubCoordinates(A),ShroudCoordinates(A))
p15 = Line(B_h,B_s)
p26 = Line(HubCoordinates(C),ShroudCoordinates(C))
p37 = Line(HubCoordinates(D),ShroudCoordinates(D))

pvolume = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
in1 = Block3D(label="in1", nni=40, nnj=12, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
in1.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-bc.lua")
in1.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#-------------------- inflow2 -------------------
A = Node(-1.0,-0.15,0.0)
B_h = LE_out_h
B_s = LE_out_s
C_h = LE_out_up_h
C_s = LE_out_up_s
D = Node(-1.0, 0.15,0.0)

p01 = Helix(HubCoordinates(A),B_h,AxisCoordinates(A),AxisCoordinates(B_h))
p12 = spline_Front_up_h
p32 = Helix(HubCoordinates(D),C_h,AxisCoordinates(D),AxisCoordinates(C_h))
p03 = Arc(HubCoordinates(A),HubCoordinates(D),AxisCoordinates(A))
p45 = Helix(ShroudCoordinates(A),B_s,AxisCoordinates(A),AxisCoordinates(B_s))
p56 = spline_Front_up_s
p76 = Helix(ShroudCoordinates(D),C_s,AxisCoordinates(D),AxisCoordinates(C_s))
p47 = Arc(ShroudCoordinates(A),ShroudCoordinates(D),AxisCoordinates(A))
p04 = Line(HubCoordinates(A),ShroudCoordinates(A))
p15 = Line(B_h,B_s)
p26 = Line(C_h,C_s)
p37 = Line(HubCoordinates(D),ShroudCoordinates(D))

pvolume = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12; 
in2 = Block3D(label="in2", nni=40, nnj=8, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
in2.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-bc.lua")

# -------------------- inflow3 -------------------
A = Node(-1.0,-0.5,0.0)
AB = Node(0.0,-0.5,0.0)
B = Node(0.05,-0.45,0.0)
C_h = LE_out_h
C_s = LE_out_s
D = Node(-1.0,-0.15,0.0)

lineAAB_h = Line(HubCoordinates(A),HubCoordinates(AB))
lineABB_h = Helix(HubCoordinates(AB),HubCoordinates(B),AxisCoordinates(AB),AxisCoordinates(B))
lineAAB_s = Line(ShroudCoordinates(A),ShroudCoordinates(AB))
lineABB_s = Helix(ShroudCoordinates(AB),ShroudCoordinates(B),AxisCoordinates(AB),AxisCoordinates(B))

p01 = Polyline([lineAAB_h,lineABB_h])
p12 = Helix(HubCoordinates(B),C_h,AxisCoordinates(B),AxisCoordinates(C_h))
p32 = Helix(HubCoordinates(D),C_h,AxisCoordinates(D),AxisCoordinates(C_h))
p03 = Arc(HubCoordinates(A),HubCoordinates(D),AxisCoordinates(A))
p45 = Polyline([lineAAB_s,lineABB_s])
p56 = Helix(ShroudCoordinates(B),C_s,AxisCoordinates(B),AxisCoordinates(C_s))
p76 = Helix(ShroudCoordinates(D),C_s,AxisCoordinates(D),AxisCoordinates(C_s))
p47 = Arc(ShroudCoordinates(A),ShroudCoordinates(D),AxisCoordinates(A))
p04 = Line(HubCoordinates(A),ShroudCoordinates(A))
p15 = Line(HubCoordinates(B),ShroudCoordinates(B))
p26 = Line(C_h,C_s)
p37 = Line(HubCoordinates(D),ShroudCoordinates(D))

pvolume = WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
in3 = Block3D(label="in3", nni=40, nnj=12, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
in3.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-bc.lua")
in3.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

# -------------------- outer1 -------------------
A_h = LE_out_up_h
A_s = LE_out_up_s
B_h = TE_up_h
B_s = TE_up_s
C = Node(0.6,1.1,0.0)
CD = Node(0.0,0.5,0.0)
D = Node(-0.3,0.5,0.0)

lineCDD_h = Line(HubCoordinates(D),HubCoordinates(CD))
lineCCD_h = Helix(HubCoordinates(CD),HubCoordinates(C),AxisCoordinates(CD),AxisCoordinates(C))
lineCDD_s = Line(ShroudCoordinates(D),ShroudCoordinates(CD))
lineCCD_s = Helix(ShroudCoordinates(CD),ShroudCoordinates(C),AxisCoordinates(CD),AxisCoordinates(C))

p01 = spline_SS_h
p12 = Helix(B_h,HubCoordinates(C),AxisCoordinates(B_h),AxisCoordinates(C))
p32 = Polyline([lineCDD_h,lineCCD_h])
p03 = Helix(A_h,HubCoordinates(D),AxisCoordinates(A_h),AxisCoordinates(D))
p45 = spline_SS_s
p56 = Helix(B_s,ShroudCoordinates(C),AxisCoordinates(B_s),AxisCoordinates(C))
p76 = Polyline([lineCDD_s,lineCCD_s])
p47 = Helix(A_s,ShroudCoordinates(D),AxisCoordinates(A_s),AxisCoordinates(D))
p04 = Line(A_h,A_s)
p15 = Line(B_h,B_s)
p26 = Line(HubCoordinates(C),ShroudCoordinates(C))
p37 = Line(HubCoordinates(D),ShroudCoordinates(D))

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
outer1 = Block3D(label="outer1", nni=50, nnj=12, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
outer1.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

# -------------------- outer2 -------------------
A = Node(0.05,-0.45,0.0)
B = Node(0.2,-0.3,0.0)
C_h = LE_out_down_h
C_s = LE_out_down_s
D_h = LE_out_h
D_s = LE_out_s

p01 = Helix(HubCoordinates(A),HubCoordinates(B),AxisCoordinates(A),AxisCoordinates(B))
p12 = Helix(HubCoordinates(B),C_h,AxisCoordinates(B),AxisCoordinates(C_h))
p32 = spline_Front_down_h
p03 = Helix(HubCoordinates(A),D_h,AxisCoordinates(A),AxisCoordinates(D_h))
p45 = Helix(ShroudCoordinates(A),ShroudCoordinates(B),AxisCoordinates(A),AxisCoordinates(B))
p56 = Helix(ShroudCoordinates(B),C_s,AxisCoordinates(B),AxisCoordinates(C_s))
p76 = spline_Front_down_s
p47 = Helix(ShroudCoordinates(A),D_s,AxisCoordinates(A),AxisCoordinates(D_s))
p04 = Line(HubCoordinates(A),ShroudCoordinates(A))
p15 = Line(HubCoordinates(B),ShroudCoordinates(B))
p26 = Line(C_h,C_s)
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
outer2 = Block3D(label="outer2", nni=8, nnj=12, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
outer2.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

# -------------------- outer3 -------------------
A = Node(0.2,-0.3,0.0)
AB = Node(0.707107,0.207107,0.0)
B = Node(0.9,0.207107,0.0)
C_h = TE_down_h
C_s = TE_down_s
D_h = LE_out_down_h
D_s = LE_out_down_s

lineAAB_h = Helix(HubCoordinates(A),HubCoordinates(AB),AxisCoordinates(A),AxisCoordinates(AB))
lineABB_h = Line(HubCoordinates(AB),HubCoordinates(B))
lineAAB_s = Helix(ShroudCoordinates(A),ShroudCoordinates(AB),AxisCoordinates(A),AxisCoordinates(AB))
lineABB_s = Line(ShroudCoordinates(AB),ShroudCoordinates(B))

p01 = Polyline([lineAAB_h,lineABB_h])
p12 = Helix(HubCoordinates(B),C_h,AxisCoordinates(B),AxisCoordinates(C_h))
p32 = spline_PS_h
p03 = Helix(HubCoordinates(A),D_h,AxisCoordinates(A),AxisCoordinates(D_h))
p45 = Polyline([lineAAB_s,lineABB_s])
p56 = Helix(ShroudCoordinates(B),C_s,AxisCoordinates(B),AxisCoordinates(C_s))
p76 = spline_PS_s
p47 = Helix(ShroudCoordinates(A),D_s,AxisCoordinates(A),AxisCoordinates(D_s))
p04 = Line(HubCoordinates(A),ShroudCoordinates(A))
p15 = Line(HubCoordinates(B),ShroudCoordinates(B))
p26 = Line(C_h,C_s)
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
outer3 = Block3D(label="outer3", nni=50, nnj=12, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
outer3.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

# -------------------- outflow1 -------------------
A_h = TE_up_h
A_s = TE_up_s
B = Node(1.707107,0.9,0.0)
C = Node(1.707107,1.207107,0.0)
CD = Node(0.707107,1.207107,0.0)
D = Node(0.6,1.1,0.0)

lineCDD_h = Line(HubCoordinates(D),HubCoordinates(CD))
lineCCD_h = Helix(HubCoordinates(CD),HubCoordinates(C),AxisCoordinates(CD),AxisCoordinates(C))
lineCDD_s = Line(ShroudCoordinates(D),ShroudCoordinates(CD))
lineCCD_s = Helix(ShroudCoordinates(CD),ShroudCoordinates(C),AxisCoordinates(CD),AxisCoordinates(C))

p01 = Helix(A_h,HubCoordinates(B),AxisCoordinates(A_h),AxisCoordinates(B))
p12 = Arc(HubCoordinates(B),HubCoordinates(C),AxisCoordinates(B))
p32 = Polyline([lineCDD_h,lineCCD_h])
p03 = Helix(A_h,HubCoordinates(D),AxisCoordinates(A_h),AxisCoordinates(D))
p45 = Helix(A_s,ShroudCoordinates(B),AxisCoordinates(A_s),AxisCoordinates(B))
p56 = Arc(ShroudCoordinates(B),ShroudCoordinates(C),AxisCoordinates(B))
p76 = Polyline([lineCDD_s,lineCCD_s])
p47 = Helix(A_s,ShroudCoordinates(D),AxisCoordinates(A_s),AxisCoordinates(D))
p04 = Line(A_h,A_s)
p15 = Line(HubCoordinates(B),ShroudCoordinates(B))
p26 = Line(HubCoordinates(C),ShroudCoordinates(C))
p37 = Line(HubCoordinates(D),ShroudCoordinates(D))

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
out1 = Block3D(label="out1", nni=50, nnj=12, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
out1.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit)
out1.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

# -------------------- outflow2 -------------------
A_h = TE_h
A_s = TE_s
B = Node(1.707107,0.707107,0.0)
C = Node(1.707107,0.9,0.0)
D_h = TE_up_h
D_s = TE_up_s

p01 = Helix(A_h,HubCoordinates(B),AxisCoordinates(A_h),AxisCoordinates(B))
p12 = Arc(HubCoordinates(B),HubCoordinates(C),AxisCoordinates(B))
p32 = Helix(D_h,HubCoordinates(C),AxisCoordinates(D_h),AxisCoordinates(C))
p03 = Helix(A_h,D_h,AxisCoordinates(A_h),AxisCoordinates(D_h))
p45 = Helix(A_s,ShroudCoordinates(B),AxisCoordinates(A_s),AxisCoordinates(B))
p56 = Arc(ShroudCoordinates(B),ShroudCoordinates(C),AxisCoordinates(B))
p76 = Helix(D_s,ShroudCoordinates(C),AxisCoordinates(D_s),AxisCoordinates(C))
p47 = Helix(A_s,D_s,AxisCoordinates(A_s),AxisCoordinates(D_s))
p04 = Line(A_h,A_s)
p15 = Line(HubCoordinates(B),ShroudCoordinates(B))
p26 = Line(HubCoordinates(C),ShroudCoordinates(C))
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
out2 = Block3D(label="out2", nni=50, nnj=6, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
out2.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit)

# -------------------- outflow3 -------------------
A_h = TE_down_h
A_s = TE_down_s
B = Node(1.707107,0.5,0.0)
C = Node(1.707107,0.707107,0.0)
D_h = TE_h
D_s = TE_s

p01 = Helix(A_h,HubCoordinates(B),AxisCoordinates(A_h),AxisCoordinates(B))
p12 = Arc(HubCoordinates(B),HubCoordinates(C),AxisCoordinates(B))
p32 = Helix(D_h,HubCoordinates(C),AxisCoordinates(D_h),AxisCoordinates(C))
p03 = Helix(A_h,D_h,AxisCoordinates(A_h),AxisCoordinates(D_h))
p45 = Helix(A_s,ShroudCoordinates(B),AxisCoordinates(A_s),AxisCoordinates(B))
p56 = Arc(ShroudCoordinates(B),ShroudCoordinates(C),AxisCoordinates(B))
p76 = Helix(D_s,ShroudCoordinates(C),AxisCoordinates(D_s),AxisCoordinates(C))
p47 = Helix(A_s,D_s,AxisCoordinates(A_s),AxisCoordinates(D_s))
p04 = Line(A_h,A_s)
p15 = Line(HubCoordinates(B),ShroudCoordinates(B))
p26 = Line(HubCoordinates(C),ShroudCoordinates(C))
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12; #+ [cluster_k,]*4;
out3 = Block3D(label="out3", nni=50, nnj=6, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
out3.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit)

# -------------------- outflow4 -------------------
A = Node(0.9,0.207107,0.0)
B = Node(1.707107,0.207107,0.0)
C = Node(1.707107,0.5,0.0)
D_h = TE_down_h
D_s = TE_down_s

p01 = Line(HubCoordinates(A),HubCoordinates(B))
p12 = Arc(HubCoordinates(B),HubCoordinates(C),AxisCoordinates(B))
p32 = Helix(D_h,HubCoordinates(C),AxisCoordinates(D_h),AxisCoordinates(C))
p03 = Helix(HubCoordinates(A),D_h,AxisCoordinates(A),AxisCoordinates(D_h))
p45 = Line(ShroudCoordinates(A),ShroudCoordinates(B))
p56 = Arc(ShroudCoordinates(B),ShroudCoordinates(C),AxisCoordinates(B))
p76 = Helix(D_s,ShroudCoordinates(C),AxisCoordinates(D_s),AxisCoordinates(C))
p47 = Helix(ShroudCoordinates(A),D_s,AxisCoordinates(A),AxisCoordinates(D_s))
p04 = Line(HubCoordinates(A),ShroudCoordinates(A))
p15 = Line(HubCoordinates(B),ShroudCoordinates(B))
p26 = Line(HubCoordinates(C),ShroudCoordinates(C))
p37 = Line(D_h,D_s)

pvolume=WireFrameVolume(p01,p12,p32,p03,p45,p56,p76,p47,p04,p15,p26,p37)
cflist = [None,]*12;
out4 = Block3D(label="out4", nni=50, nnj=12, nnk=20,
               parametric_volume=pvolume, cf_list=cflist,
               fill_condition=initialCond)
out4.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit)
out4.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

# finally, just join all of the abutting blocks together
identify_block_connections()
