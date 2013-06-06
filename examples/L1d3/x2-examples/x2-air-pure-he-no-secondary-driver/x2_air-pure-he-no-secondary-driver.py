# File: x2_air-pure-he-no-secondary-driver.py
# This is an L1d simulation of the X2 expansion tube without 
# a secondary driver section and with a pure He primary driver
# and matching orifice plate.
# Condition is the 11.3 km/s air condition used by Umar Sheikh.
# This is input file is based on work done by David Gildfind and Umar Sheikh.
# 
# Full piston dynamics.
# Equilibrium calculation.
# 
# Chris James (c.james4@uq.edu.au)

# need to import math as math.pi is used later
import math

gdata.title = 'L1d3 - Simulation of X2 with lightweight piston'

t_finish=0.030     # Simulation end time [s].
t_switch=0.020     # Simulation switch to finer time steps [s].
t_fine=1.0e-5*3.0     # Size of finer time steps [s].

mesh_scale=2.0;

#fill pressures

#currently using sample calculation I did in pitot (initial-l1d-comparison.txt)

p_res_fill=6.85e6	             			# Reservoir fill pressure [Pa].
p_drv_fill=92.8e3		              		# Driver fill pressure [Pa].
prim_dia_burst=27.9e6			    		# Primary diaphragm burst pressure.

shk_tube_fill=3000.0			     		# Shock tube fill pressure[Pa].
acc_tube_fill=10.0					# Acceleration tube fill pressure [Pa].


# Diaphragm locations 

x_sec_diaphragm=4.810+3.418                     # x-location of secondary diaphragm [m].
#x_tert_diaphragm=4.810+5.976                    # x-location of tertiary diaphragm [m].
x_acc_tube_exit=14.795                   # Acceleration tube exit; =13.789 without nozzle, =14.795m with nozzle.

# below x distances pulled from the x2_tube l1d file dave had

# Define the tube walls
# Reservoir and compression tube:
add_break_point(-3.890, 0.3160, 1)
add_break_point(-0.990, 0.3160, 1)
add_break_point(-0.970, 0.2440, 1)
add_break_point(-0.370, 0.2440, 1)
add_break_point(-0.350, 0.1600, 1)
add_break_point(-0.157, 0.1600, 1)
add_break_point(-0.010, 0.2568, 1)
# Tunnel downstream of compression tube:
add_break_point(4.600, 0.2568, 1) # Beginning of area change at PD; assumed linear change to next break point.
add_break_point(4.700, 0.0850, 1) # End of area change to PD.
# Orifice plate just after the PD which is at 4.810
add_break_point(4.810, 0.0850, 1) # End of area change to PD. This is the reference for the primary diaphragm location (from Dave Gildfind's PhD.)
add_break_point(4.880, 0.0650, 1) # End of area change to PD.
add_break_point(4.950, 0.0650, 1) # End of area change to PD.
add_break_point(5.020, 0.0850, 1) # End of area change to PD.

#add_break_point(x_acc_tube_exit, 0.0850, 0) # Analysis tube exit. 
#
# EITHER Add the nozzle: (refer Michael Scott PhD Thesis, 2006, Table 5.3)
# THIS IS THE AT EXIT, SET ABOVE. ADD OFFSET OF 0.958 FROM HERE

noz_offset = -0.010
add_break_point(13.395+noz_offset, 0.085, 0) 
add_break_point(13.418422+noz_offset, 0.08696, 0)
add_break_point(13.447994+noz_offset, 0.089496, 0)
add_break_point(13.484767+noz_offset, 0.092956, 0)
add_break_point(13.529836+noz_offset, 0.097684, 0)
add_break_point(13.584295+noz_offset, 0.103878, 0)
add_break_point(13.649202+noz_offset, 0.111644, 0)
add_break_point(13.725529+noz_offset, 0.121042, 0)
add_break_point(13.814119+noz_offset, 0.132066, 0)
add_break_point(13.915628+noz_offset, 0.144532, 0)
add_break_point(14.030471+noz_offset, 0.157942, 0)
add_break_point(14.158763+noz_offset, 0.17141, 0)
add_break_point(14.300261+noz_offset, 0.183712, 0)
add_break_point(14.454306+noz_offset, 0.193498, 0)
add_break_point(14.619769+noz_offset, 0.199654, 0)
add_break_point(14.795+noz_offset, 0.20168, 0)
nozzle_exit=14.795+noz_offset
#
# OR Add a straight adaptor to the dumptank:
# add_break_point(13.789, 0.0850, 1) # Analysis tube exit. 
#
###################################################################################
# Notes about the tube configuration:                                             #
#                                                                                 #
# Configuration may have mylar diaphragms at x = 8.228m, x = 10.786	          #
# I.e. these are the capstan locations where the normally used tubes are joined.  #
#                                                                                 #
# The normal fixed length 85mm diameter tubes end at x=13.169m.                   #
#                                                                                 #
# If the facility is to be run in tube mode, then a 0.620m long dumptank adaptor  #
# is installed, in which case the tube ends at x=13.169+0.620=13.789m.            #
#                                                                                 #
# If the facility is to be run with the nozzle, 0.226m of the nozzle is straight, #
# therefore the straight section ends at 13.169+0.226=13.395m.                    #    
# The contoured part of the nozzle is 1.4m long, therefore the nozzle exit        #
# is located at x=13.395+1.4=14.795m.                                             #
###################################################################################

#stuff to calculate the mass fractions of He and Ar in the primary driver (as L1D stupidly deals in mass fractions)
#remember to comment out orifice plate code if 100%He driver is not in use
perc_He=100.0
perc_Ar=0.0
R_Ar=208.0
R_He=2077.0
R_m=(perc_He+perc_Ar)/(perc_He/R_He+perc_Ar/R_Ar)
mf_He=(perc_He*R_m)/(R_He*100.0)
mf_Ar=(perc_Ar*R_m)/(R_Ar*100.0)

## For LUT analysis
create_gas_file(model="ideal gas", species=['Ar', 'He', 'N2', 'air'], 
                fname="gas-model.lua", lut_file="cea-lut-air.lua.gz")
species_list = select_gas_model(fname="gas-model.lua")
print "species_list=", species_list

##Gas compositions (mass fractions) (used Umar's conventions here as it had separate stuff for the gases)
#species_list =               [  LUT_air  Ar      He      N2    Perf_air] 
reservoir_mfs =               [  0.0,     0.0,    0.0,    0.0,  1.0     ]
driver_mfs =                  [  0.0,     mf_Ar, mf_He,	  0.0,  0.0     ]
test_gas_mfs = 	              [  1.0,     0.0,    0.0,    0.0,  0.0     ]
accel_gas_mfs =               [  1.0,     0.0,    0.0,    0.0,  0.0     ]

# Create the gas-path

left_wall = VelocityEnd(x0=-3.890, v=0.0)

res_gas = GasSlug(p=p_res_fill, u=0.0, T=296.0, nn=100*mesh_scale, adaptive=0, to_end_L=0,
                     to_end_R=1, cluster_strength=0.0, hcells=1,
                     viscous_effects=1, adiabatic_flag=0,
                     massf=reservoir_mfs,
                     label='compressed air to push the piston')

piston = Piston(m=10.524, d=0.2568, xL0=0.0, xR0=0.221, v0=0.0,
                front_seal_f=0.4, front_seal_area=0.020*0.2568*math.pi,
                is_restrain=0, with_brakes=0,
                x_buffer=4.5895, hit_buffer=0,
                label='single stage x2 piston')

driver_gas = GasSlug(p=p_drv_fill, u=0.0, T=296.0, nn=100*mesh_scale, adaptive=0, to_end_L=0,
                     to_end_R=1, cluster_strength=1.05, hcells=1,
                     viscous_effects=1, adiabatic_flag=0,
                     massf=driver_mfs,
                     label='compressed helium driver gas')

primary_diaphragm = Diaphragm(x0=4.810, p_burst=prim_dia_burst, is_burst=0, dt_hold=10.0e-6, dt_blend=0.0, dx_blend=0.0)

test_gas = GasSlug(p=shk_tube_fill, u=0.0, T=296.0, nn=100*mesh_scale, adaptive=0, 
		    to_end_L=1, to_end_R=1, cluster_strength=1.05, hcells=1,
                    viscous_effects=1, adiabatic_flag=0,
                    massf=test_gas_mfs, label='test-gas in the shock tube')


# Secondary diaphragm burst pressure is from the burst pressure of a single aluminium foil diaphragm.
secondary_diaphragm = Diaphragm(x0=x_sec_diaphragm, p_burst=18.0e3, is_burst=0, dt_hold=10.0e-6, dt_blend=0.0, dx_blend=0.0)

#diaphragm_two = Piston(m=794E-6, d=0.085, xL0=x_sec_diaphragm, xR0=x_sec_diaphragm, v0=0.0,
#                          type_of_piston=0, is_restrain=1, p_restrain=0.5e6, with_brakes=0,
#                          label='ideal piston acting as secondary diaphragm', f_decay=25000.0,
#                          mass_limit=1.0e-7)

#diaphragm_three = Piston(m=794E-6, d=0.085, xL0=x_tert_diaphragm, xR0=x_tert_diaphragm, v0=0.0,
#                          type_of_piston=0, is_restrain=1, p_restrain=0.5e6, with_brakes=0,
#                          label='ideal piston acting as tertiary diaphragm', f_decay=25000.0,
#                          mass_limit=1.0e-7)

accelerator_gas = GasSlug(p=acc_tube_fill, u=0.0, T=296.0, nn=300*mesh_scale, nnmax=1000, adaptive=1, 
                          to_end_L=1, to_end_R=0, cluster_strength=1.05, hcells=1,
                          viscous_effects=1, adiabatic_flag=0,
                          massf=accel_gas_mfs,
                          label='Air in the acceleration tube')

right_free = FreeEnd(x0=nozzle_exit)

assemble_gas_path(left_wall, res_gas, piston, driver_gas, primary_diaphragm, test_gas, secondary_diaphragm, accelerator_gas, right_free)

# Add some global data
gdata.n = 1000

# Add a loss region
add_loss_region(-0.36, 0.15, 3.1)
add_loss_region(4.6, 4.81, 0.7)

# Set some time-stepping parameters
gdata.dt_init = 1.0e-10
gdata.max_time = t_finish
gdata.max_step = 25000000
gdata.cfl = 0.25
gdata.t_order = 2
gdata.x_order = 2
add_dt_plot(0.0, 2.0e-4, 2.0e-4)
add_dt_plot(t_switch, t_fine, t_fine/10)

# Define history locations
add_history_loc(4.600)   # Compression tube pressure immediately upstream of primary diaphragm. 0
add_history_loc(4.810)   # Compression tube pressure at primary diaphragm. 1
add_history_loc(4.810+2.577)   # PCB transducer sd1. 2
add_history_loc(4.810+2.810)   # PCB transducer sd2. 3
add_history_loc(4.810+3.043)   # PCB transducer sd3. 4
add_history_loc(4.810+4.231)   # PCB transducer st1. 5
add_history_loc(4.810+4.746)   # PCB transducer st2. 6
add_history_loc(4.810+5.260)   # PCB transducer st3. 7
add_history_loc(4.810+6.437)   # PCB transducer at1. 8
add_history_loc(4.810+6.615)   # PCB transducer at2. 9
add_history_loc(4.810+6.796)   # PCB transducer at3. 10
add_history_loc(4.810+7.590)   # PCB transducer at4. 11
add_history_loc(4.810+7.846)   # PCB transducer at5. 12
add_history_loc(4.810+8.096)   # PCB transducer at6. 13

add_history_loc(13.395+noz_offset)  # End of straight tube part of facility. 14

add_history_loc(nozzle_exit)   # Nozzle exit 15

