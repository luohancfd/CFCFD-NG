"""
2D Compressor Blade Standard Condition 10

Hannes Wojciak, Paul Petrie-Repar
    February 2008: Original implementation
Peter J.
    March 2008: Clean-up and periodic boundary condition
    03-Sep-2008: Port to Eilmer3
Peter Blyton
    June 2011: Geometry cleaned up and simplified.
Peter J.
    March 2014: replace UDF BCs to allow MPI simulation.

"""

# from cfpylib.geom.path import Polyline2, Spline2

# ---------------- First, set the global data ----------------------
gdata.title = "inviscid Euler for 2D-sc10"
gdata.dimensions = 2
# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
gdata.viscous_flag = 0 # inviscid simulation
gdata.gasdynamic_update_scheme = "euler"
gdata.max_time = 0.300
gdata.max_step = 800000
gdata.dt_plot = 0.020
gdata.dt = 1.0e-7

# -----------flow conditions -------------------
# Stagnation condition that leads to inflow.
p_tot = 100.0e3 # Pa
T_tot = 300.0 # degree K
totalCond = FlowCondition(p=p_tot, T=T_tot)
alpha = math.radians(55.0)
inflow_BC = SubsonicInBC(inflow_condition=totalCond,
                         direction_type="uniform",
                         direction_vector=[math.cos(alpha), math.sin(alpha), 0.0],
                         label="INLET")
# Compute expected outflow.
gma = 1.4
Rgas = 287.0  # J/kg.K
a_tot = math.sqrt(gma*Rgas*T_tot)
M_exit = 0.45 # this is our controlling parameter
T0_T = 1 + (gma-1.0)/2.0 * M_exit * M_exit
p0_p = T0_T**(gma/(gma-1.0))
print "p0_p=", p0_p, "T0_T=", T0_T
p_exit = p_tot / p0_p
T_exit = T_tot / T0_T
u_exit = M_exit * a_tot / math.sqrt(T0_T)
print "p_exit=", p_exit, "T_exit=", T_exit, "u_exit=", u_exit
outflow_BC = FixedPOutBC(Pout=p_exit, label="OUTLET")
initialCond = FlowCondition(p=p_exit, u=u_exit, T=T_exit)

# Periodic boundary conditions need mapping functions
# from ghost-cell locations to source-cell locations. 
# Note that, because the MappedCellBCs store boundary-specific information,
# we must assign a unique MappedCellBC object on each block boundary.
# We cannot reuse the one object for many boundaries.
def subtract_one_from_y(x, y, z): return x, y-1.0, z
def add_one_to_y(x, y, z): return x, y+1.0, z

# Mesh setup parameters
mrf = 6 # Mesh refinement factor, must be an even integer
clust_chord = RobertsClusterFunction(1, 1, 1.3) # clustering along chord
clust_blade_top = RobertsClusterFunction(1, 0, 1.05) # normal to chord, top
clust_blade_bottom = RobertsClusterFunction(0, 1, 1.05) # normal to chord, bottom
clust_LE_surface = RobertsClusterFunction(1, 0, 1.05) # along surface toward LE
clust_LE_chord = RobertsClusterFunction(0, 1, 1.02) # clustering toward LE in LE blocks

# Suction and pressure surfaces of blades defined using coordinate data.
profile_SS = Spline2("sc10_inner1.dat")
profile_Front_up = Spline2("sc10_inner2.dat")
profile_Front_down = Spline2("sc10_inner3.dat")
profile_Front_down.reverse()
profile_PS = Spline2("sc10_inner4.dat")
profile_PS.reverse()

# Nodes on and surrounding the blade surface
TE = profile_SS.eval(1.0)
TE_up = TE + Vector(-0.06, 0.12)
LE_up = Node(0.007375, 0.038160, label="LE_up")
LE_out_up = Node(-0.1, 0.07, label="LE_out_up")
LE = Node(0.0, 0.0, label="LE")
LE_out = Node(-0.05, -0.09, label="LE_out")
LE_down = Node(0.026541, 0.015230, label="LE_down")
LE_out_down = Node(0.09, -0.07, label="LE_out_down")
TE_down = Node(0.75, 0.6, label="TE_down")

# ---------------path definitions-------------------
SS = Node(0.18,0.44) # Node for spline at SS
spline_SS = Spline([LE_out_up,SS,TE_up]) # outer spline at SS of profile
Fup = Node(-0.1,-0.01) # Nodes for spline in front of profile
spline_Front_up = Spline([LE_out,Fup,LE_out_up]) # outer spline in front of profile
Fdown = Node(0.02,-0.1) # Nodes for spline in front of profile
spline_Front_down = Spline([LE_out,Fdown,LE_out_down]) # outer spline in front of profile
PS = Node(0.45,0.34) # Nodes for spline at PS of profile
spline_PS = Spline([LE_out_down,PS,TE_down]) # outer spline at PS of profile

# ---------------------- inner1 --------------------------
path_s = profile_SS
path_n = spline_SS
path_w = Line(LE_up,LE_out_up)
path_e = Line(TE,TE_up)

cflist = [clust_chord, clust_blade_top, clust_chord, clust_blade_top]
patch = make_patch(path_n, path_e, path_s, path_w)
inner1 = Block2D(label="inner1", nni=mrf*8, nnj=mrf, psurf=patch,
                cf_list=cflist, fill_condition=initialCond)

# ---------------------- inner2 --------------------------
path_s = Line(LE_out,LE)
path_n = Line(LE_out_up,LE_up)
path_w = spline_Front_up
path_e = profile_Front_up

patch = make_patch(path_n, path_e, path_s, path_w)
cflist = [clust_blade_bottom, clust_LE_surface, clust_LE_chord, None]
inner2 = Block2D(label="inner2", nni=inner1.nnj, nnj=int(mrf*1.5), psurf=patch,
                cf_list=cflist, fill_condition=initialCond)

# ---------------------- inner3 --------------------------
path_s = spline_Front_down
path_n = profile_Front_down
path_w = Line(LE_out,LE)
path_e = Line(LE_out_down,LE_down)

patch = make_patch(path_n, path_e, path_s, path_w)
cflist = [clust_LE_surface, clust_blade_bottom, None, clust_LE_chord]
inner3 = Block2D(label="inner3", nni=int(mrf*1.5), nnj=inner2.nni, psurf=patch,
                cf_list=cflist, fill_condition=initialCond)

# -------------------inner4 ------------
path_s = spline_PS
path_n = profile_PS
path_w = Line(LE_out_down,LE_down)
path_e = Line(TE_down,TE)

patch = make_patch(path_n, path_e, path_s, path_w)
cflist = [clust_chord, clust_blade_bottom, clust_chord, clust_blade_bottom]
inner4 = Block2D(label="inner4", nni=mrf*7, nnj=inner3.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)

# -------------------- inflow1 -------------------
A = Node(-1.0, 0.15)
B = LE_out_up
C = Node(-0.3,0.5)
D = Node(-1.0,0.5)

path_s = Line(A,B)
path_n = Line(D,C)
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
in1 = Block2D(label="in1", nni=mrf*6, nnj=mrf*2, psurf=patch,
               fill_condition=initialCond)
in1.bc_list[WEST] = inflow_BC
in1.bc_list[NORTH] = MappedCellBC(ghost_cell_trans_fn=subtract_one_from_y)

# -------------------- inflow2 -------------------
A = Node(-1.0,-0.15)
B = LE_out
C = LE_out_up
D = Node(-1.0, 0.15)

path_s = Line(A,B)
path_n = Line(D,C)
path_w = Line(A,D)
path_e = spline_Front_up

patch = make_patch(path_n, path_e, path_s, path_w)
in2 = Block2D(label="in2", nni=in1.nni, nnj=inner2.nnj, psurf=patch,
               fill_condition=initialCond)
in2.bc_list[WEST] = inflow_BC

# -------------------- inflow3 -------------------
A = Node(-1.0,-0.5)
AB = Node(0.0,-0.5)
B = Node(0.05,-0.45)
C = LE_out
D = Node(-1.0,-0.15)

path_s = Polyline2([A,AB,B])
path_n = Line(D,C)
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
in3 = Block2D(label="in3", nni=in2.nni, nnj=mrf*2, psurf=patch,
               fill_condition=initialCond)
in3.bc_list[WEST] = inflow_BC
in3.bc_list[SOUTH] = MappedCellBC(ghost_cell_trans_fn=add_one_to_y)

# -------------------- outer1 -------------------
A = LE_out_up
B = TE_up
C = Node(0.6,1.1)
CD = Node(0.0,0.5)
D = Node(-0.3,0.5)

path_s = spline_SS
path_n = Polyline2([D,CD,C])
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
cflist = [None, None, clust_chord, None]
outer1 = Block2D(label="outer1", nni=inner1.nni, nnj=in1.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
outer1.bc_list[NORTH] = MappedCellBC(ghost_cell_trans_fn=subtract_one_from_y)

# -------------------- outer2 -------------------
A = Node(0.05,-0.45)
B = Node(0.2,-0.3)
C = LE_out_down
D = LE_out

path_s = Line(A,B)
path_n = spline_Front_down
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
outer2 = Block2D(label="outer2", nni=inner3.nni, nnj=in3.nnj, psurf=patch,
               fill_condition=initialCond)
outer2.bc_list[SOUTH] = MappedCellBC(ghost_cell_trans_fn=add_one_to_y)

# -------------------- outer3 -------------------
A = Node(0.2,-0.3)
AB = Node(0.707107,0.207107)
B = Node(0.9,0.207107)
C = TE_down
D = LE_out_down

path_s = Polyline2([A,AB,B])
path_n = spline_PS
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
cflist = [clust_chord, None, None, None]
outer3 = Block2D(label="outer3", nni=inner4.nni, nnj=outer2.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
outer3.bc_list[SOUTH] = MappedCellBC(ghost_cell_trans_fn=add_one_to_y)

# -------------------- outflow1 -------------------
A = TE_up 
B = Node(1.707107,0.9)
C = Node(1.707107,1.207107)
CD = Node(0.707107,1.207107)
D = Node(0.6,1.1)

path_s = Line(A,B)
path_n = Polyline2([D,CD,C])
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
out1 = Block2D(label="out1", nni=mrf*8, nnj=outer1.nnj, psurf=patch,
                fill_condition=initialCond)
out1.bc_list[EAST] = outflow_BC
out1.bc_list[NORTH] = MappedCellBC(ghost_cell_trans_fn=subtract_one_from_y)

# -------------------- outflow2 -------------------
A = TE
B = Node(1.707107,0.707107,0.0)
C = Node(1.707107,0.9,0.0)
D = TE_up

path_s = Line(A,B)
path_n = Line(D,C)
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
cflist = [None, None, None, clust_blade_top]
out2 = Block2D(label="out2", nni=out1.nni, nnj=inner1.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
out2.bc_list[EAST] = outflow_BC

# -------------------- outflow3 -------------------
A = TE_down
B = Node(1.707107,0.5,0.0)
C = Node(1.707107,0.707107,0.0)
D = TE

path_s = Line(A,B)
path_n = Line(D,C)
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)
cflist = [None, None, None, clust_blade_bottom]
out3 = Block2D(label="out3", nni=out2.nni, nnj=inner4.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
out3.bc_list[EAST] = outflow_BC

# -------------------- outflow4 -------------------
A = Node(0.9,0.207107,0.0)
B = Node(1.707107,0.207107,0.0)
C = Node(1.707107,0.5,0.0)
D = TE_down

path_s = Line(A,B)
path_n = Line(D,C)
path_w = Line(A,D)
path_e = Line(B,C)

patch = make_patch(path_n, path_e, path_s, path_w)

out4 = Block2D(label="out4", nni=out3.nni, nnj=outer3.nnj, psurf=patch,
                fill_condition=initialCond)
out4.bc_list[EAST] = outflow_BC
out4.bc_list[SOUTH] = MappedCellBC(ghost_cell_trans_fn=add_one_to_y)

identify_block_connections()

#------------------- Presentation -----------------
sketch.xaxis(-1.0, 2.0, 0.5, -0.1)
sketch.yaxis(-0.5, 1.5, 0.5, -0.1)
sketch.window(-1.0, -0.5, 2.0, 2.5, 0.02, 0.02, 0.20, 0.20)
