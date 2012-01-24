"""
Plane turbine cascade for validation with experimental data.

Kiock, R et al. (1986). 'The Transonic Flow Through a Plane Turbine Cascade as Mea-
sured in Four European Wind Tunnels'. In: Journal of Engineering for Gas Turbines
and Power 108, p. 277.

From this paper, experimental results in Fig 6, BS (wind tunnel in Braunschweig, West Germany)

Peter Blyton
    July 2011: Original implementation.
    October 2011: Viscous simulation for better validation.

"""

#---------------------------------------------------
# Global data
#---------------------------------------------------
gdata.title = "Viscous Simulation for Kiock turbine cascade"
gdata.viscous_flag = 1
gdata.t_order = 1
gdata.max_time = 0.1
gdata.max_step = 800000
gdata.dt_plot = 0.020
gdata.dt = 1.0e-12

#---------------------------------------------------
# Flow conditions
#---------------------------------------------------
select_gas_model(model='ideal gas', species=['air'])
p_tot = 100.0e3 # Pa
T_tot = 300.0 # K
gma = 1.4
Rgas = 287.0  # J/kg.K
a_tot = math.sqrt(gma*Rgas*T_tot)
M_exit = 0.78
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
STAGGER_ANGLE = (33.3/180.0)*math.pi
PITCH = 0.71

#---------------------------------------------------
# Mesh setup parameters
#---------------------------------------------------
mrf = 4 # Mesh refinement factor

suct_div = 0.8 # Divisions for suction surface blocks
suct_div2 = 0.98

pres_div1 = 0.08 # Divisions for pressure surface blocks
pres_div2 = 0.965

clust_surf = 1.01 # clustering normal to surface

#---------------------------------------------------
# General path and node setup
#---------------------------------------------------
# Full suction and pressure surface paths
suct_surface = Spline2("suct_surface.dat")
suct_surface.rotate_about_zaxis(-STAGGER_ANGLE)
pres_surface = Spline2("press_surface.dat")
pres_surface.rotate_about_zaxis(-STAGGER_ANGLE)

# Chord is line of which blade "sits" on, starting at (0, 0)
chord_vector = Vector(1.0, 0.0)
chord_vector.rotate_about_zaxis(-STAGGER_ANGLE)
chord_normal = chord_vector.clone().rotate_about_zaxis(math.pi/2.0)
pitch_vector = Vector(0.0, PITCH/2.0)

#---------------------------------------------------
# Front block on suction surface
#---------------------------------------------------
suct_front_surface = suct_surface.clone()
suct_front_surface.t1=suct_div # Block path is section of suction surface path

# Spline for north path of block
LE_up = suct_front_surface.eval(0) - 0.04*chord_vector
mid_1 = suct_front_surface.eval(0.1) + 0.14*pitch_vector
mid_2 = suct_front_surface.eval(0.34) + 0.1*pitch_vector
mid_3 = suct_front_surface.eval(0.6) + 0.04*chord_normal
TE_up = suct_front_surface.eval(1) + 0.04*chord_normal + 0.01*chord_vector
suct_front_north = Spline([LE_up, mid_1, mid_2, mid_3, TE_up])

suct_front_east = Line(suct_front_surface.eval(1), suct_front_north.eval(1))
suct_front_west = Line(suct_front_surface.eval(0), suct_front_north.eval(0))

patch = make_patch(suct_front_north, suct_front_east, suct_front_surface, suct_front_west)
cflist = [None, RobertsClusterFunction(1, 0, clust_surf), 
          None, RobertsClusterFunction(1, 0, clust_surf)]
bclist = [SlipWallBC(), SlipWallBC(), AdiabaticBC(), SlipWallBC()]

suct_front_block = Block2D(label="suct_front", nni=7*mrf, nnj=2*mrf, 
               cf_list=cflist, psurf=patch, fill_condition=initialCond, 
               bc_list=bclist)

#---------------------------------------------------
# Middle block on suction surface
#---------------------------------------------------
suct_mid_surface = suct_surface.clone()
suct_mid_surface.t0=suct_div # Block path is section of suction surface path
suct_mid_surface.t1=suct_div2
TE_out_mid = suct_mid_surface.eval(1) + 0.04*chord_normal

suct_mid_north = Line(TE_up, TE_out_mid)
suct_mid_east = Line(suct_mid_surface.eval(1), TE_out_mid)

patch = make_patch(suct_mid_north, suct_mid_east, suct_mid_surface, suct_front_east)
cflist = [None, RobertsClusterFunction(1, 0, clust_surf), 
          None, RobertsClusterFunction(1, 0, clust_surf)]
bclist = [SlipWallBC(), SlipWallBC(), AdiabaticBC(), SlipWallBC()]

suct_mid_block = Block2D(label="suct_mid", nni=3*mrf, nnj=suct_front_block.nnj,
               cf_list=cflist, psurf=patch, fill_condition=initialCond, 
               bc_list=bclist)

#---------------------------------------------------
# Rear block on suction surface
#---------------------------------------------------
suct_rear_surface = suct_surface.clone()
suct_rear_surface.t0=suct_div2 # Block path is section of suction surface path

# Spline for north path of block
TE_out = suct_rear_surface.eval(1) + 0.036*chord_vector - 0.02*chord_normal
mid_1 = suct_rear_surface.eval(0.5) + Vector(0.033, 0)
suct_rear_north = Spline([TE_out_mid, mid_1, TE_out])

suct_rear_east = Line(suct_rear_surface.eval(1), suct_rear_north.eval(1))

patch = make_patch(suct_rear_north, suct_rear_east, suct_rear_surface, suct_mid_east)
cflist = [None, RobertsClusterFunction(1, 0, clust_surf), 
          None, RobertsClusterFunction(1, 0, clust_surf)]
bclist = [SlipWallBC(), SlipWallBC(), AdiabaticBC(), SlipWallBC()]

suct_rear_block = Block2D(label="suct_rear", nni=mrf, nnj=suct_front_block.nnj,
               cf_list=cflist, psurf=patch, fill_condition=initialCond, 
               bc_list=bclist)

#---------------------------------------------------
# Front block on pressure surface
#---------------------------------------------------
pres_front_surface = pres_surface.clone()
pres_front_surface.t1=pres_div1 # Block path is section of suction surface path
pres_front_surface.reverse()

# Spline for west path of block
LE_down = pres_front_surface.eval(0)- 0.04*chord_normal + 0.01*chord_vector
mid_1 = pres_front_surface.eval(0.5) + Vector(-0.04, 0.0)
mid_2 = pres_front_surface.eval(0.5) + Vector(-0.02, -0.05)
pres_front_west = Spline([LE_down, mid_2, mid_1, LE_up])

pres_front_south = Line(pres_front_west.eval(0), pres_front_surface.eval(0))

patch = make_patch(suct_front_west, pres_front_surface, pres_front_south, pres_front_west)
cflist = [RobertsClusterFunction(0, 1, clust_surf), None, 
          RobertsClusterFunction(0, 1, clust_surf), None]
bclist = [SlipWallBC(), AdiabaticBC(), SlipWallBC(), SlipWallBC()]

pres_front_block = Block2D(label="pres_front", nni=suct_front_block.nnj, nnj=2*mrf,
               cf_list=cflist, psurf=patch, fill_condition=initialCond, bc_list=bclist)

#---------------------------------------------------
# Middle block on pressure surface
#---------------------------------------------------
pres_mid_surface = pres_surface.clone()
pres_mid_surface.t0=pres_div1 # Block path is section of suction surface path
pres_mid_surface.t1=pres_div2

# Spline for south path of block
TE_down = pres_mid_surface.eval(1) - 0.04*chord_normal
mid_1 = pres_mid_surface.eval(0.43) - 0.040*chord_normal
mid_2 = pres_mid_surface.eval(0.6) - 0.030*chord_normal
pres_mid_south = Spline([LE_down, mid_1, TE_down])

pres_mid_east = Line(pres_mid_south.eval(1), pres_mid_surface.eval(1))

patch = make_patch(pres_mid_surface, pres_mid_east, pres_mid_south, pres_front_south)
cflist = [None, RobertsClusterFunction(0, 1, clust_surf), 
          None, RobertsClusterFunction(0, 1, clust_surf)]
bclist = [AdiabaticBC(), SlipWallBC(), SlipWallBC(), SlipWallBC()]

pres_mid_block = Block2D(label="pres_mid", nni=8*mrf, nnj=pres_front_block.nni,
               cf_list=cflist, psurf=patch, fill_condition=initialCond, bc_list=bclist)

#---------------------------------------------------
# Rear block on pressure surface
#---------------------------------------------------
pres_rear_surface = pres_surface.clone()
pres_rear_surface.t0=pres_div2 # Block path is section of suction surface path

# Spline for south path of block
mid_1 = pres_rear_surface.eval(0.56) - 0.13*pitch_vector
pres_rear_south = Spline([TE_down, mid_1, TE_out])

pres_rear_east = suct_rear_east.clone()
pres_rear_east.reverse()

patch = make_patch(pres_rear_surface, pres_rear_east, pres_rear_south, pres_mid_east)
cflist = [None, RobertsClusterFunction(0, 1, clust_surf), 
          None, RobertsClusterFunction(0, 1, clust_surf)]
bclist = [AdiabaticBC(), SlipWallBC(), SlipWallBC(), SlipWallBC()]

pres_rear_block = Block2D(label="pres_rear", nni=mrf, nnj=pres_mid_block.nnj,
               cf_list=cflist, psurf=patch, fill_condition=initialCond, 
               bc_list=bclist)

#---------------------------------------------------
# Inflow 1 block
#---------------------------------------------------
in1_top_left = Vector(-0.8, PITCH/2.0)
in1_top_right = in1_top_left + Vector(0.78, 0)
in1_bottom_left = in1_top_left + Vector(0, -0.2)
in1_north = Line(in1_top_left, in1_top_right)
in1_east = Line(LE_up, in1_top_right)
in1_south = Line(in1_bottom_left, LE_up)
in1_west = Line(in1_bottom_left, in1_top_left)

patch = make_patch(in1_north, in1_east, in1_south, in1_west, "AO")

in1_block = Block2D(label="in1", nni=5*mrf, nnj=2*mrf,
               psurf=patch, fill_condition=initialCond)
in1_block.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-kiock.lua", label="INLET")
in1_block.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# Inflow 2 block
#---------------------------------------------------
in2_bottom_left = in1_bottom_left + Vector(0, -0.4)
in2_south = Line(in2_bottom_left, LE_down)
in2_west = Line(in2_bottom_left, in1_bottom_left)

patch = make_patch(in1_south, pres_front_west, in2_south, in2_west, "AO")

in2_block = Block2D(label="in2", nni=in1_block.nni, nnj=pres_front_block.nnj,
               psurf=patch, fill_condition=initialCond)
in2_block.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-kiock.lua", label="INLET")

#---------------------------------------------------
# Inflow 3 block
#---------------------------------------------------
in3_bottom_left = Vector(-0.8, -PITCH/2.0)
in3_bottom_right = in3_bottom_left + Vector(0.78, 0)
in3_east = Line(in3_bottom_right, LE_down)
in3_south = Line(in3_bottom_left, in3_bottom_right)
in3_west = Line(in3_bottom_left, in2_bottom_left)

patch = make_patch(in2_south, in3_east, in3_south, in3_west, "AO")

in3_block = Block2D(label="in3", nni=in2_block.nni, nnj=in1_block.nnj,
               psurf=patch, fill_condition=initialCond)
in3_block.set_BC(WEST, USER_DEFINED, filename="udf-subsonic-kiock.lua", label="INLET")
in3_block.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# Outflow 1 block
#---------------------------------------------------
out1_top_right = suct_surface.eval(1) + Vector(1, 0.6*PITCH)
out1_top_left = out1_top_right + Vector(-0.98, 0)
out1_bottom_right = out1_top_right + Vector(0, -0.2)
out1_north = Line(out1_top_left, out1_top_right)
out1_east = Line(out1_bottom_right, out1_top_right)
out1_south = Line(TE_up, out1_bottom_right)
out1_west = Line(TE_up, out1_top_left)

patch = make_patch(out1_north, out1_east, out1_south, out1_west)

out1_block = Block2D(label="out1", nni=5*mrf, nnj=in1_block.nnj,
               psurf=patch, fill_condition=initialCond)
out1_block.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")
out1_block.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# Outflow 2 block
#---------------------------------------------------
out2_bottom_right = out1_bottom_right + Vector(0, -0.2)
out2_east = Line(out2_bottom_right, out1_bottom_right)
out2_south = Line(TE_out_mid, out2_bottom_right)
out2_west = suct_mid_north.clone()
out2_west.reverse()

patch = make_patch(out1_south, out2_east, out2_south, out2_west, "AO")

out2_block = Block2D(label="out2", nni=out1_block.nni, nnj=suct_mid_block.nni,
               psurf=patch, fill_condition=initialCond)
out2_block.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")

#---------------------------------------------------
# Outflow 2low block
#---------------------------------------------------
out2low_bottom_right = out2_bottom_right + Vector(0, -0.1)
out2low_east = Line(out2low_bottom_right, out2_bottom_right)
out2low_south = Line(TE_out, out2low_bottom_right)
out2low_west = suct_rear_north.clone()
out2low_west.reverse()

patch = make_patch(out2_south, out2low_east, out2low_south, out2low_west, "AO")

out2_block = Block2D(label="out2low", nni=out1_block.nni, nnj=suct_rear_block.nni,
               psurf=patch, fill_condition=initialCond)
out2_block.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")

#---------------------------------------------------
# Outflow 3 block
#---------------------------------------------------
out3_bottom_right = suct_surface.eval(1) + Vector(1, -0.4*PITCH)
out3_bottom_left = out3_bottom_right + Vector(-0.98, 0)
out3_east = Line(out3_bottom_right, out2low_bottom_right)
out3_south = Line(out3_bottom_left, out3_bottom_right)
out3_west = Line(out3_bottom_left, TE_out)

patch = make_patch(out2_south, out3_east, out3_south, out3_west, "AO")

out3_block = Block2D(label="out3", nni=out2_block.nni, nnj=in3_block.nnj,
               psurf=patch, fill_condition=initialCond)
out3_block.set_BC("EAST", "FIXED_P_OUT", Pout=p_exit, label="OUTLET")
out3_block.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# Top block
#---------------------------------------------------
mid_1 = in1_top_right + Vector(0.44, 0.01)
top_north = Spline([in1_top_right, mid_1, out1_top_left])

patch = make_patch(top_north, out1_west, suct_front_north, in1_east, "AO")

top_block = Block2D(label="top", nni=suct_front_block.nni, nnj=in1_block.nnj,
               psurf=patch, fill_condition=initialCond)
top_block.set_BC(NORTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# Lower 1 block
#---------------------------------------------------
low1_south = top_north.clone()
low1_south.translate(-2.0*pitch_vector)
low1_south.t1 = pres_div2 - 0.1
low1_east = Line(low1_south.eval(1), TE_down)

patch = make_patch(pres_mid_south, low1_east, low1_south, in3_east, "AO")

low1_block = Block2D(label="low1", nni=pres_mid_block.nni, nnj=in3_block.nnj,
               psurf=patch, fill_condition=initialCond)
low1_block.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

#---------------------------------------------------
# Lower 2 block
#---------------------------------------------------
low2_south = top_north.clone()
low2_south.translate(-2.0*pitch_vector)
low2_south.t0 = low1_south.t1

patch = make_patch(pres_rear_south, out3_west, low2_south, low1_east)

low2_block = Block2D(label="low2", nni=pres_rear_block.nni, nnj=low1_block.nnj,
               psurf=patch, fill_condition=initialCond)
low2_block.set_BC(SOUTH, USER_DEFINED, filename="udf-periodic-bc.lua")

identify_block_connections()

#---------------------------------------------------
# Presentation
#---------------------------------------------------
sketch.xaxis(-1.0, 2.0, 0.5, -0.1)
sketch.yaxis(-2.0, 1.0, 0.5, -0.1)
sketch.window(-1.0, -2.0, 3.0, 3.0, 0.02, 0.02, 0.20, 0.20)

