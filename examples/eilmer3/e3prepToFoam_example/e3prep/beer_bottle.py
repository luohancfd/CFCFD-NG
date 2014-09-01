## \file beer_bottle.py
#
# An example gridding file for MECH4480

# First set up the bare minimum parameters to allow e3prep to solve
select_gas_model(model='ideal gas', species=['air'])
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)
inflow  = FlowCondition(p=95.84e3, u=1000.0, v=0.0, T=1103.0)

# Now we need to create all of the nodes that will make up our geometry
A = Node(-0.07, 0.0, label="A")
B = Node(-0.05, 0.0, label="B")
C = Node(0.0, 0.0, label="C")
D = Node(0.005, 0.012, label="D")
E = Node(0.1, 0.03, label="E")
F = Node(0.202, 0.03, label="F")
G = Node(0.207, 0.0, label="G")
H = Node(0.3, 0.0, label="H")
I = Node(-0.07, 0.1, label="I")
J = Node(-0.05, 0.1, label="J")
K = Node(0.1, 0.1, label="K")
L = Node(0.202, 0.1, label="L")
M = Node(0.3, 0.1, label="M")
N = Node(0.3, 0.03, label="N")

# Define Bezier Control points
CD_b1 = Node(0.0, 0.006, label="CD-b1")
DJ_b1 = Node(-0.008, 0.085, label="DJ-b1")
JD_b1 = Node(-0.008, 0.085, label="JD-b1")
GF_b1 = Node(0.207, 0.027, label="GF-b1")
DE_b1 = Node(0.0064, 0.012, label="DE-b1")
DE_b2 = Node(0.0658, 0.0164, label="DE-b2")
DE_b3 = Node(0.0727, 0.0173, label="DE-b2")

# Now we join our nodes to create lines that will be used to form our blocks
AB = Line(A, B); BC = Line(B, C); GH = Line(G,H) # lower boundary along x-axis
CD = Bezier([C, CD_b1, D]); DE = Bezier([D, DE_b1, DE_b2, DE_b3, E]);
EF = Line(E, F); GF = Bezier([G, GF_b1, F])	# lower boundary along bottle
IJ = Line(I, J); JK = Line(J, K); KL = Line(K, L); LM = Line(L, M) # upper boundary
AI = Line(A, I); BJ = Line(B, J); DJ = Bezier([D, DJ_b1, J]);
JD = Bezier([J, JD_b1, D]); EK = Line(E, K); FL = Line(F, L); 
NM = Line(N, M); HN = Line(H, N);   # vertical lines
FN = Line(F, N) # line between BL_4 and BL_5 

# Define the blocks, boundary conditions and set the discretisation.
n0 = 10; n1 = 5; n2 = 10; n3 = 20; n4 = 20; n5 = 15; n6 = 8

BL_0 = Block2D(make_patch(IJ, BJ, AB, AI), nni=n1, nnj=n0,
                fill_condition=initial, label="BL_0")
BL_1 = Block2D(make_patch(JD, CD, BC, BJ), nni=n2, nnj=n0,
                fill_condition=initial, label="BL_1")
BL_2 = Block2D(make_patch(JK, EK, DE, DJ), nni=n3, nnj=n2,
                fill_condition=initial, label="BL_2")
BL_3 = Block2D(make_patch(KL, FL, EF, EK), nni=n4, nnj=n2,
                fill_condition=initial, label="BL_3")
BL_4 = Block2D(make_patch(LM, NM, FN, FL), nni=n5, nnj=n2,
                fill_condition=initial, label="BL_4")
BL_5 = Block2D(make_patch(FN, HN, GH, GF), nni=n5, nnj=n6,
                fill_condition=initial, label="BL_5")
identify_block_connections()

BL_0.bc_list[WEST] = ExtrapolateOutBC(label='OF_inlet_00')
BL_0.bc_list[SOUTH] = ExtrapolateOutBC(label='OF_symmetry_00')
BL_1.bc_list[SOUTH] = ExtrapolateOutBC(label='OF_symmetry_00')
BL_0.bc_list[NORTH] = ExtrapolateOutBC(label='OF_inlet_02')
BL_2.bc_list[NORTH] = ExtrapolateOutBC(label='OF_inlet_02')
BL_3.bc_list[NORTH] = ExtrapolateOutBC(label='OF_inlet_02')
BL_4.bc_list[NORTH] = ExtrapolateOutBC(label='OF_inlet_02')
BL_4.bc_list[EAST] = ExtrapolateOutBC(label='OF_outlet_00')
BL_5.bc_list[EAST] = ExtrapolateOutBC(label='OF_outlet_00')
BL_5.bc_list[SOUTH] = ExtrapolateOutBC(label='OF_symmetry_00')
BL_1.bc_list[EAST] = ExtrapolateOutBC(label='OF_wall_00')
BL_2.bc_list[SOUTH] = ExtrapolateOutBC(label='OF_wall_01')
BL_3.bc_list[SOUTH] = ExtrapolateOutBC(label='OF_wall_01')
BL_5.bc_list[WEST] = ExtrapolateOutBC(label='OF_wall_02')

sketch.prefer_bc_labels_on_faces()

# This is to make a nice *.svg file at the end
sketch.xaxis(-0.1, 0.35, 0.05, -0.005)
sketch.yaxis(0.0, 0.15, 0.05, -0.1)
sketch.window(0.0, -0.1, 0.4, 0.2, 0.05, 0.05, 0.17, 0.17)
