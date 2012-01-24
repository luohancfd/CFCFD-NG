# bar.py
# PJ
# 14-Dec-2006
# 03-Feb-2010 ported to Eilmer3 examples

gdata.title = "Bar gauge M=4.76 in air."
print gdata.title
# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])
# Define flow conditions: low pressure ambient with M=4.76 inflow
initial = FlowCondition(p=30.0,    u=0.0,    v=0.0, T=300.0)
inflow  = FlowCondition(p=100.0e3, u=1653.0, v=0.0, T=300.0)

# Geometry
R = 5.0e-3    # radius of bar in metres
a = Node(-2*R, 0.0, label="A")
b = Node(-2*R, R,   label="B")
c = Node(-2*R, 6*R, label="C")
d = Node( 0.0, 0.0, label="D")
e = Node( 0.0, R,   label="E")
f = Node( 0.0, 6*R, label="F")
g = Node( 6*R, R,   label="G")
h = Node( 6*R, 6*R, label="H")

ad=Line(a,d); be=Line(b,e); cf=Line(c,f); eg=Line(e,g); fh=Line(f,h)
ab=Line(a,b); bc=Line(b,c); de=Line(d,e); ef=Line(e,f); gh=Line(g,h)

# Define the blocks, boundary conditions and set the discretisation.
nx0 = 120; nx2 = 120; ny0 = 40; ny1=80
cfy = RobertsClusterFunction(1, 0, 1.2)
cfx = RobertsClusterFunction(1, 0, 1.1)
blk_0 = Block2D(make_patch(be, de, ad, ab), nni=nx0, nnj=ny0,
                fill_condition=initial, label="Front",
                hcell_list=[(nx0,1),(nx0,5),(nx0,10)],
                xforce_list=[0,1,0,0])
blk_1 = Block2D(make_patch(cf, ef, be, bc), nni=nx0, nnj=ny1,
                cf_list=[None,cfy,None,cfy],
                fill_condition=initial, label="Outer")
blk_2 = Block2D(make_patch(fh, gh, eg, ef), nni=nx2, nnj=ny1,
                cf_list=[cfx,cfy,cfx,cfy],
                fill_condition=initial, label="After")
identify_block_connections()
blk_0.bc_list[WEST] = SupInBC(inflow)
blk_1.bc_list[WEST] = SupInBC(inflow)
blk_2.bc_list[EAST] = ExtrapolateOutBC()

# We can set individual attributes of the global data object.
gdata.axisymmetric_flag = 1
gdata.flux_calc = ADAPTIVE
gdata.max_time = 50.0e-6  # seconds
gdata.max_step = 15000
gdata.dt = 2.0e-8
gdata.dt_plot = 5.0e-6
gdata.dt_history = 0.5e-6

sketch.xaxis(-0.010, 0.030, 0.010, -0.002)
sketch.yaxis( 0.000, 0.030, 0.010, -0.002)
sketch.window(-0.010, 0.0, 0.030, 0.040, 0.05, 0.05, 0.17, 0.17)
