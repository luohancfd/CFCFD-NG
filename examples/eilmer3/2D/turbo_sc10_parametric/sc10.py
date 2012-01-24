"""
2D Compressor Blade Standard Condition 10 parametric setup.

Hannes Wojciak, Paul Petrie-Repar
    February 2008: Original implementation.
Peter J.
    March 2008: Clean-up and periodic boundary condition.
    03-Sep-2008: Port to Eilmer3.
Peter Blyton
    March 2011: Blade profile defined using functions, sc10_blade_profile.py.
    April 2011: Block mesh set up to handle small perturbations.
    June 2011: Geometry cleaned up and simplified.

"""

from sc10_blade_profile import *
from cfpylib.geom.transform_pyfunc import rotate_pyfunc
from cfpylib.geom.path import Polyline2

#---------------------------------------------------
# Global data
#---------------------------------------------------
gdata.title = "Inviscid Euler Simulation for 2D sc10"
gdata.t_order = 1
gdata.max_time = 0.3
gdata.max_step = 800000
gdata.dt_plot = 0.020
gdata.dt = 1.0e-7

#---------------------------------------------------
# Flow conditions
#---------------------------------------------------
select_gas_model(model='ideal gas', species=['air'])
p_tot = 100.0e3 # Pa
T_tot = 300.0 # K
gma = 1.4
Rgas = 287.0  # J/kg.K
a_tot = math.sqrt(gma*Rgas*T_tot)
M_exit = 0.45
T0_T = 1 + (gma-1.0)/2.0 * M_exit * M_exit
p0_p = T0_T**(gma/(gma-1.0))
print "p0_p=", p0_p, "T0_T=", T0_T
p_exit = p_tot / p0_p
T_exit = T_tot / T0_T
u_exit = M_exit * a_tot / math.sqrt(T0_T)
print "p_exit=", p_exit, "T_exit=", T_exit, "u_exit=", u_exit
initialCond = FlowCondition(p=p_exit, u=u_exit, T=T_exit)

#---------------------------------------------------
# Geometric parameters
#---------------------------------------------------
STAGGER_ANGLE = math.pi/4.0
PITCH = 1.0

#---------------------------------------------------
# Mesh setup parameters
#---------------------------------------------------
mrf = 6 # Mesh refinement factor, must be an even integer
division = 0.03 # fraction of chord length for block division in C-mesh
clust_chord = RobertsClusterFunction(1, 1, 1.3) # clustering along chord
clust_blade_top = RobertsClusterFunction(1, 0, 1.05) # normal to chord, top
clust_blade_bottom = RobertsClusterFunction(0, 1, 1.05) # normal to chord, bottom
clust_LE_surface = RobertsClusterFunction(1, 0, 1.01) # along surface toward LE
clust_LE_chord = RobertsClusterFunction(0, 1, 1.01) # clustering toward LE in LE blocks

#---------------------------------------------------
# General path and node setup
#---------------------------------------------------
# Suction surface paths
sc10_top_rotated = rotate_pyfunc(sc10_top, "z", STAGGER_ANGLE)
profile_Front_up = PyFunctionPath(sc10_top_rotated, "", 0, division)
profile_SS = PyFunctionPath(sc10_top_rotated, "", division, 1)

# Pressure surface paths
sc10_bottom_rotated = rotate_pyfunc(sc10_bottom, "z", STAGGER_ANGLE)
profile_Front_down = PyFunctionPath(sc10_bottom_rotated, "", 0, division)
profile_PS = PyFunctionPath(sc10_bottom_rotated, "", division, 1)

# Nodes at leading and trailing edges
LE = profile_Front_up.eval(0.0)
LE_up = profile_SS.eval(0.0)
LE_down = profile_PS.eval(0.0)
TE = profile_SS.eval(1.0)
chord_normal = TE.clone().rotate_about_zaxis(math.pi/2.0)

# Nodes surrounding leading edge
LE_out = -0.1*TE
LE_out_up = LE_up + 0.1*chord_normal
LE_out_down = LE_down - 0.1*chord_normal

# Nodes surrounding trailing edge
TE_up = TE + Vector(-0.06, 0.12)
TE_down = TE + Vector(0.06, -0.12)

# Nodes bounding the flow domain
IN_top = Vector(-1.0, PITCH/2.0)
IN_bottom = Vector(-1.0, -PITCH/2.0)
LE_top = LE + Vector(0.0, PITCH/2.0)
LE_bottom = LE - Vector(0.0, PITCH/2.0)
TE_top = TE + Vector(0.0, PITCH/2.0)
TE_bottom = TE - Vector(0.0, PITCH/2.0)
OUT_top = TE_top + Vector(1.0, 0.0)
OUT_bottom = TE_bottom + Vector(1.0, 0.0)

# Spline above suction surface
SS = profile_SS.eval(0.5) + 0.1*chord_normal
spline_SS = Spline([LE_out_up, SS, TE_up])

# Splines in front of leading edge
Fup = LE + Vector(-0.1, 0)
spline_Front_up = Spline([LE_out, Fup, LE_out_up])
Fdown = LE + Vector(0, -0.1)
spline_Front_down = Spline([LE_out, Fdown, LE_out_down])

# Splines below pressure surface
PS = profile_PS.eval(0.5) - 0.12*chord_normal
spline_PS = Spline([LE_out_down, PS, TE_down])

#---------------------------------------------------
# inner1 block
#---------------------------------------------------
inner1_east = Line(TE, TE_up)
inner1_west = Line(LE_up, LE_out_up)
patch = make_patch(spline_SS, inner1_east, profile_SS, inner1_west, "AO")
cflist = [clust_chord, clust_blade_top, clust_chord, clust_blade_top]

inner1 = Block2D(label="inner1", nni=mrf*8, nnj=mrf, psurf=patch,
                cf_list=cflist, fill_condition=initialCond)

#---------------------------------------------------
# inner2 block
#---------------------------------------------------
inner2_south = Line(LE_out, LE)
patch = make_patch(inner1_west.reverse(), profile_Front_up, inner2_south, spline_Front_up)
cflist = [clust_blade_bottom, clust_LE_surface, clust_LE_chord, None]

inner2 = Block2D(label="inner2", nni=inner1.nnj, nnj=int(mrf*1.5), psurf=patch,
               cf_list=cflist, fill_condition=initialCond)

#---------------------------------------------------
# inner3 block
#---------------------------------------------------
inner3_east = Line(LE_out_down, LE_down)
patch = make_patch(profile_Front_down, inner3_east, spline_Front_down, inner2_south)
cflist = [clust_LE_surface, clust_blade_bottom, None, clust_LE_chord]

inner3 = Block2D(label="inner3", nni=int(mrf*1.5), nnj=inner2.nni, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)

#---------------------------------------------------
# inner4 block
#---------------------------------------------------
inner4_east = Line(TE_down, TE)
patch = make_patch(profile_PS, inner4_east, spline_PS, inner3_east)
cflist = [clust_chord, clust_blade_bottom, clust_chord, clust_blade_bottom]

inner4 = Block2D(label="inner4", nni=mrf*7, nnj=inner3.nnj, psurf=patch,
                cf_list=cflist, fill_condition=initialCond)

#---------------------------------------------------
# inflow1 block
#---------------------------------------------------
inflow1_lower_left = Vector(-1.0, 0.15)
inflow1_upper_right = Vector(-0.3,PITCH/2.0)
inflow1_north = Line(IN_top, inflow1_upper_right)
inflow1_east = Line(LE_out_up, inflow1_upper_right)
inflow1_south = Line(inflow1_lower_left, LE_out_up)
inflow1_west= Line(inflow1_lower_left, IN_top)
patch = make_patch(inflow1_north, inflow1_east, inflow1_south, inflow1_west)

inflow1 = Block2D(label="inflow1", nni=mrf*6, nnj=mrf*2, psurf=patch,
                fill_condition=initialCond)
inflow1.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-sc10.lua", label="INLET")
inflow1.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# inflow2 block
#---------------------------------------------------
inflow2_lower_left = Vector(-1.0, -0.15)
inflow2_south = Line(inflow2_lower_left, LE_out)
inflow2_west = Line(inflow2_lower_left, inflow1_lower_left)
patch = make_patch(inflow1_south, spline_Front_up, inflow2_south, inflow2_west)

inflow2 = Block2D(label="inflow2", nni=inflow1.nni, nnj=inner2.nnj, psurf=patch,
               fill_condition=initialCond)
inflow2.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-sc10.lua", label="INLET")

#---------------------------------------------------
# inflow3 block
#---------------------------------------------------
inflow3_east = Line(LE_bottom, LE_out)
inflow3_south = Line(IN_bottom, LE_bottom)
inflow3_west = Line(IN_bottom, inflow2_lower_left)
patch = make_patch(inflow2_south, inflow3_east, inflow3_south, inflow3_west)

inflow3 = Block2D(label="inflow3", nni=inflow2.nni, nnj=mrf*2, psurf=patch,
               fill_condition=initialCond)
inflow3.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-sc10.lua", label="INLET")
inflow3.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# outer1 block
#---------------------------------------------------
outer1_upper_right = LE_top + 0.8*TE
outer1_north = Polyline2([inflow1_upper_right, LE_top, outer1_upper_right])
outer1_east = Line(TE_up, outer1_upper_right)
patch = make_patch(outer1_north, outer1_east, spline_SS, inflow1_east)
cflist = [None, None, clust_chord, None]

outer1 = Block2D(label="outer1", nni=inner1.nni, nnj=inflow1.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
outer1.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# outer2 block
#---------------------------------------------------
outer2_lower_right = LE_bottom + 0.3*TE
outer2_east = Line(outer2_lower_right, LE_out_down)
outer2_south = Line(LE_bottom, outer2_lower_right)
patch = make_patch(spline_Front_down, outer2_east, outer2_south, inflow3_east)

outer2 = Block2D(label="outer2", nni=inner3.nni, nnj=inflow3.nnj,
               psurf=patch, fill_condition=initialCond)
outer2.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# outer3 block
#---------------------------------------------------
outer3_lower_right = TE_bottom + Vector(0.2, 0.0)
outer3_east = Line(outer3_lower_right, TE_down)
outer3_south = Polyline2([outer2_lower_right, TE_bottom, outer3_lower_right])
patch = make_patch(spline_PS, outer3_east, outer3_south, outer2_east)
cflist = [clust_chord, None, None, None]

outer3 = Block2D(label="outer3", nni=inner4.nni, nnj=outer2.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
outer3.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# outflow1 block
#---------------------------------------------------
outflow1_lower_right = OUT_top - Vector(0, 0.3)
outflow1_north = Polyline2([outer1_upper_right ,TE_top, OUT_top])
outflow1_east = Line(outflow1_lower_right, OUT_top)
outflow1_south = Line(TE_up, outflow1_lower_right)
patch = make_patch(outflow1_north, outflow1_east, outflow1_south, outer1_east)

outflow1 = Block2D(label="outflow1", nni=mrf*8, nnj=outer1.nnj, psurf=patch,
               fill_condition=initialCond)
outflow1.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")
outflow1.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# outflow2 block
#---------------------------------------------------
outflow2_lower_right = OUT_top - Vector(0, PITCH/2.0)
outflow2_east = Line(outflow2_lower_right, outflow1_lower_right)
outflow2_south = Line(TE, outflow2_lower_right)
patch = make_patch(outflow1_south, outflow2_east, outflow2_south, inner1_east)
cflist = [None, None, None, clust_blade_top]

outflow2 = Block2D(label="outflow2", nni=outflow1.nni, nnj=inner1.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
outflow2.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")

#---------------------------------------------------
# outflow3 block
#---------------------------------------------------
outflow3_lower_right = OUT_bottom + Vector(0.0, 0.3)
outflow3_east = Line(outflow3_lower_right, outflow2_lower_right)
outflow3_south = Line(TE_down, outflow3_lower_right)
patch = make_patch(outflow2_south, outflow3_east, outflow3_south, inner4_east, "AO")
cflist = [None, None, None, clust_blade_bottom]

outflow3 = Block2D(label="outflow3", nni=outflow2.nni, nnj=inner4.nnj, psurf=patch,
               cf_list=cflist, fill_condition=initialCond)
outflow3.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")

#---------------------------------------------------
# outflow4 block
#---------------------------------------------------
outflow4_east = Line(OUT_bottom, outflow3_lower_right)
outflow4_south = Line(outer3_lower_right, OUT_bottom)
patch = make_patch(outflow3_south, outflow4_east, outflow4_south, outer3_east)

outflow4 = Block2D(label="outflow4", nni=outflow3.nni, nnj=outer3.nnj,
               psurf=patch, fill_condition=initialCond)
outflow4.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")
outflow4.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

identify_block_connections()

#---------------------------------------------------
# Presentation
#---------------------------------------------------
sketch.xaxis(-1.0, 2.0, 0.5, -0.1)
sketch.yaxis(-0.5, 1.5, 0.5, -0.1)
sketch.window(-1.0, -0.5, 2.2, 2.7, 0.02, 0.02, 0.20, 0.20)

