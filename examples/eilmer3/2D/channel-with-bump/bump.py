# bump.py
# PJ, 09-Sep-2013
# This test appeared in papers on Euler solvers for supersonic flows.
# For example:
# S. Eidelman, P. Colella and R.P. Shreeve (1984)
# Application of the Godunov method and its second-order extension to
# cascade flow modelling.
# AIAA Journal Vol. 22 No. 11 pp  

gdata.title = "Channel with circular-arc bump."
print gdata.title
p_inf = 101.3e3 # Pa
T_inf = 288.0 # degree K
import math
a_inf = math.sqrt(1.4*287.0*T_inf)
u_inf = 1.65 * a_inf # m/s
select_gas_model(model='ideal gas', species=['air'])
inflow  = FlowCondition(p=p_inf, u=u_inf, T=T_inf)

L = 1.0
h = 0.04 * L
a0 = Vector(0.0,0.0); a1 = Vector(0.0,L)
b0 = Vector(L,0.0); b1 = Vector(L,L)
c0 = Vector(1.5*L,h)
d0 = Vector(2.0*L,0.0); d1 = Vector(2.0*L,L)
e0 = Vector(3.0*L,0.0); e1 = Vector(3.0*L,L)
rcfx0 = RobertsClusterFunction(0,1,1.2)
rcfx1 = RobertsClusterFunction(1,1,1.2)
rcfx2 = RobertsClusterFunction(1,0,1.2)
rcfy = RobertsClusterFunction(1,0,1.2)
ni0 = 64; nj0 = 64 # We'll scale discretization off these values
factor = 1.0
ni0 = int(ni0*factor); nj0 = int(nj0*factor)

blk0 = SuperBlock2D(CoonsPatch(a0,b0,b1,a1),
                    nni=ni0, nnj=nj0, nbi=16, nbj=2,
                    bc_list=[SlipWallBC(),None,SlipWallBC(),SupInBC(inflow)],
                    cf_list=[rcfx0,rcfy,rcfx0,rcfy],
                    fill_condition=inflow, label="B0")
blk1 = SuperBlock2D(make_patch(Line(b1,d1),Line(d0,d1),
                               Arc3(b0,c0,d0),Line(b0,b1)),
                     nni=ni0, nnj=nj0, nbi=16, nbj=2,
                     bc_list=[SlipWallBC(),None,SlipWallBC(),None],
                     cf_list=[rcfx1,rcfy,rcfx1,rcfy],
                     fill_condition=inflow, label="B1")
blk2 = SuperBlock2D(CoonsPatch(d0,e0,e1,d1),
                    nni=ni0, nnj=nj0, nbi=16, nbj=2,
                    bc_list=[SlipWallBC(),ExtrapolateOutBC(),SlipWallBC(),None],
                    cf_list=[rcfx2,rcfy,rcfx2,rcfy],
                    fill_condition=inflow, label="B2")
identify_block_connections()

# Do a little more setting of global data.
gdata.flux_calc = ADAPTIVE
gdata.gasdynamic_update_scheme = "classic-rk3"
gdata.cfl = 1.0
gdata.max_time = 10.0*L/u_inf # long enough, tunnel has 15ms steady time
gdata.max_step = 50000
gdata.dt = 1.0e-6
gdata.dt_plot = gdata.max_time

sketch.xaxis(0.0, 3.0, 0.5, -0.10)
sketch.yaxis(0.0, 1.0, 0.5, -0.10)
sketch.window(0.0, 0.0, 3.0, 3.0, 0.05, 0.05, 0.25, 0.25)
