# the_plain_bottle.py

select_gas_model(model='ideal gas', species=['air'])
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)

# Create the nodes that define key points for our geometry.
A = Node(-0.07, 0.0, label="A"); B = Node(-0.05, 0.0, label="B")
C = Node(0.0, 0.0, label="C"); D = Node(0.005, 0.012, label="D")
E = Node(0.1, 0.03, label="E"); F = Node(0.202, 0.03, label="F")
G = Node(0.207, 0.0, label="G"); H = Node(0.3, 0.0, label="H")
I = Node(-0.07, 0.1, label="I"); J = Node(-0.05, 0.1, label="J")
K = Node(0.1, 0.1, label="K"); L = Node(0.202, 0.1, label="L")
M = Node(0.3, 0.1, label="M"); N = Node(0.3, 0.03, label="N")

# Some interior Bezier control points
CD_b1 = Node(0.0, 0.006, label="CD-b1")
DJ_b1 = Node(-0.008, 0.075, label="DJ-b1")
GF_b1 = Node(0.207, 0.027, label="GF-b1")
DE_b1 = Node(0.0064, 0.012, label="DE-b1")
DE_b2 = Node(0.0658, 0.0164, label="DE-b2")
DE_b3 = Node(0.0727, 0.0173, label="DE-b3")

# Now, we join our nodes to create lines that will be used to form our blocks.
AB = Line(A, B); BC = Line(B, C); GH = Line(G,H) # lower boundary along x-axis
CD = Bezier([C, CD_b1, D]) # top of bottle
DE = Bezier([D, DE_b1, DE_b2, DE_b3, E]) # neck of bottle
EF = Line(E, F) # side of bottle
GF = Bezier([G, GF_b1, F],"GF",0.0,1.0,1) # bottom, with arc-length parameterization
# Upper boundary of domain
IJ = Line(I, J); JK = Line(J, K); KL = Line(K, L); LM = Line(L, M) 
# Lines to divide the gas flow domain into blocks.
AI = Line(A, I); BJ = Line(B, J); DJ = Bezier([D, DJ_b1, J])
JD = DJ.copy(direction=-1); EK = Line(E, K); FL = Line(F, L); 
NM = Line(N, M); HN = Line(H, N); FN = Line(F, N)

# Define the blocks, boundary conditions and set the discretisation.
n0 = 10; n1 = 4; n2 = 20; n3 = 20; n4 = 20; n5 = 12; n6 = 8

BL_0 = Block2D(make_patch(IJ, BJ, AB, AI), nni=n1, nnj=n0,
               fill_condition=initial, label="[0]")
BL_1 = Block2D(make_patch(JD, CD, BC, BJ), nni=n2, nnj=n0,
               fill_condition=initial, label="[1]")
BL_2 = Block2D(make_patch(JK, EK, DE, DJ), nni=n3, nnj=n2,
               fill_condition=initial, label="[2]")
BL_3 = Block2D(make_patch(KL, FL, EF, EK), nni=n4, nnj=n2,
               fill_condition=initial, label="[3]")
BL_4 = Block2D(make_patch(LM, NM, FN, FL), nni=n5, nnj=n2,
               fill_condition=initial, label="[4]")
BL_5 = Block2D(make_patch(FN, HN, GH, GF), nni=n5, nnj=n6,
               fill_condition=initial, label="[5]")
identify_block_connections()

# Make a nicely-scaled SVG file at the end.
sketch.xaxis(-0.1, 0.3, 0.05, -0.05)
sketch.yaxis(0.0, 0.10, 0.05, -0.01)
sketch.window(-0.1, 0.0, 0.3, 0.4, 0.02, 0.05, 0.20, 0.23)
