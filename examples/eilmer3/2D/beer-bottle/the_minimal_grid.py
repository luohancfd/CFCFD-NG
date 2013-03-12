# the_minimal_grid.py

select_gas_model(model='ideal gas', species=['air'])
initial = FlowCondition(p=5955.0,  u=0.0,    v=0.0, T=304.0)

# Create the nodes that define key points for our geometry.
E = Node(0.1, 0.03, label="E"); F = Node(0.202, 0.03, label="F")
K = Node(0.1, 0.1, label="K"); L = Node(0.202, 0.1, label="L")

p = make_patch(Line(K,L), Line(F,L), Line(E,F), Line(E,K))
BL_3 = Block2D(p, nni=20, nnj=20, fill_condition=initial, label="[3]")

# Make a nicely-scaled SVG file at the end.
sketch.xaxis(0.1, 0.2, 0.05, -0.010)
sketch.yaxis(0.0, 0.1, 0.05, -0.030)
sketch.window(0.1, 0.0, 0.2, 0.1, 0.05, 0.05, 0.10, 0.10)
