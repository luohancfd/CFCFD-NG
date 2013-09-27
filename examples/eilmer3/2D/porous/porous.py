## \file porous.py
## \brief testing flow through a porous medium
## \author EJF, 12-Sep-2013

from math import cos, sin, tan, sqrt, pi
from cfpylib.grid.shock_layer_surface import *

job_title = "Porous test case."
print job_title

gdata.title = job_title
gdata.axisymmetric_flag = 1

#
# 1. Setup the gas model
#
# select_gas_model(model='ideal gas', species=['air'])

# Hydrogen-air combustion model.
# Species included: H, H2, O, O2, N, N2, OH, NO, H2O, HO2, NO2, H2O2, HNO
species = select_gas_model(model='thermally perfect gas', species=['H', 'H2', 'O', 'O2', 'N', 'N2', 'OH', 'NO', 'H2O', 'HO2', 'NO2', 'H2O2', 'HNO'])
set_reaction_scheme("Jachimowski_92.lua", reacting_flag=0)
gm = get_gas_model_ptr()
nsp = gm.get_number_of_species()
ntm = gm.get_number_of_modes()

# air inflow
T_inf = 1103.0
u_inf = 1000.0
p_inf = 95.84e3
# flowfield species mass fractions
massf_inf = [ 0.0 ] * gm.get_number_of_species()
massf_inf[species.index('N2')] = 0.767
massf_inf[species.index('O2')] = 0.233

# fuel inflow (high p, low vel.y)
p_fuel = 5.0e5
v_fuel = 1.0
# fuel species mass fractions
massf_H2 = [ 0.0 ] * gm.get_number_of_species()
massf_H2[species.index('H2')] = 1.0

Q = Gas_data(gm)
for itm in range(ntm):
    Q.T[itm] = T_inf
Q.p = p_inf
for isp,massf in enumerate(massf_inf):
    Q.massf[isp] = massf
gm.eval_thermo_state_rhoT(Q)
#M_inf = u_inf / Q.a
rho_inf = Q.rho

#
# 2. Define flow conditions
#
initial = FlowCondition(p=p_inf/10.0,  u=0.0,    v=0.0, T=T_inf, massf=massf_inf)
inflow  = FlowCondition(p=p_inf, u=u_inf, v=0.0, T=T_inf, massf=massf_inf)
high_p_inflow = FlowCondition(p=p_fuel, u=0.0, v=v_fuel, T=T_inf, massf=massf_H2)

#
# 3. Define the geometry
#

# origin
o = Node(0.0,0.1, label='o')

# nodes
a = Node(0.3,0.1, label='a')
b = Node(0.6,0.1, label='b')
c = Node(0.9,0.1, label='c')
d = Node(0.9,0.4, label='d')
e = Node(0.6,0.4, label='e')
f = Node(0.3,0.4, label='f')
g = Node(0.0,0.4, label='g')
h = Node(0.3,0.08, label='h')
i = Node(0.6,0.08, label='i')

oa = Line(o, a); ab = Line(a, b); bc = Line(b, c) # lower boundary
gf = Line(g, f); fe = Line(f, e); ed = Line(e, d) # upper boundary
og = Line(o, g); af = Line(a, f); be = Line(b, e); cd = Line(c, d) # vertical lines
ha = Line(h, a); ib = Line(i, b);	# extra vertical lines
hi = Line(h, i);			# extra horizontal line

''' the geometry looks like this...
      _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
     |         |         |         |
sup  |         |         |         |
inBC |         |         |         |---> flow left to right
     |_ _ _ _ _|_ _ _ _ _|_ _ _ _ _|
               |_ _ _ _ _| porous block
                 high p
                 inflow of fuel through lower wall
'''


#
# 4. Define the blocks, boundary conditions and set the discretisation
#
nnx = 20; nny= 20
nbx = 4; nby = 1

blk_0 = Block2D(make_patch(gf, af, oa, og), nni=nnx, nnj=nny,
                bc_list=[SlipWallBC(), AdjacentBC(), SlipWallBC(), SupInBC(inflow)],
                fill_condition=initial, label="BLOCK-0")

blk_1 = Block2D(make_patch(fe, be, ab, af), nni=nnx, nnj=nny,
                bc_list=[SlipWallBC(), AdjacentBC(), AdjacentBC(), AdjacentBC()],
                fill_condition=initial, label="BLOCK-1")

blk_2 = Block2D(make_patch(ed, cd, bc, be), nni=nnx, nnj=nny,
                bc_list=[SlipWallBC(), ExtrapolateOutBC(), SlipWallBC(), AdjacentBC()],
                fill_condition=initial, label="BLOCK-2")

blk_3 = Block2D(make_patch(ab, ib, hi, ha), nni=nnx, nnj=nny/4,
                bc_list=[AdjacentBC(), SlipWallBC(), SupInBC(high_p_inflow), SlipWallBC()],
                fill_condition=initial, label="BLOCK-3")

identify_block_connections()

# Do a little more setting of global data.
gdata.udf_file = "udf-process.lua"
gdata.udf_source_vector_flag = 1 	# 0=standard case; 1=porous zone case
gdata.viscous_flag = 0
gdata.flux_calc = ADAPTIVE
gdata.max_time = 5.0e-3  # seconds
gdata.max_step = 30000
gdata.dt = 1.0e-8
gdata.dt_max = 1.0e-6
# gdata.cfl = 0.1
gdata.dt_plot = 1.5e-3

sketch.xaxis(0.0, 1.0, 0.2, -0.05)
sketch.yaxis(0.0, 0.6, 0.2, -0.04)
sketch.window(0.0, 0.0, 1.0, 1.0, 0.05, 0.05, 0.17, 0.17)

