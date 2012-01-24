# dn2.py
gdata.title = \
   "Drummond tunnel M4 nozzle P4 = 3.25MPa N2, P1 = 30kPa N2, 06-Jul-05"
gdata.case_id = 27
select_gas_model(fname='cea-lut-N2.lua.gz')

# Define the tube walls.
gdata.n = 4000
add_break_point(-3.785,  0.0585)    # upstream-end of the driver tube
add_break_point(-3.035,  0.0585)
add_break_point(-3.015,  0.0620)    # steel-diaphragm station
# Note that there is no steel diaphragm in this simulation.
add_break_point( 0.000,  0.0620)    # downstrem end of shock tube
add_break_point( 0.043,  0.0620)    # start of contraction to throat
add_break_point( 0.080,  0.0220)    # start of throat
add_break_point( 0.100,  0.0220, 1) # start of nozzle (conical shape)
add_break_point( 0.2653, 0.0700, 1) # start of parallel test section
add_break_point( 0.30,   0.0700)
#
add_loss_region(-3.050, -3.000, 0.5) # at steel-diaphragm station
add_loss_region( 0.050,  0.120, 0.5) # at nozzle throat
#
gdata.T_nominal = 296.0

# Create the gas-path.
left_end = VelocityEnd(x0=-3.785, v=0.0)
driver_gas = GasSlug(p=3.25e6, u=0.0, T=296.0, nn=150,
                     viscous_effects=1, hcells=1,
                     label="driver gas")
cs = GasInterface(x0=-3.015)
test_gas = GasSlug(p=30.0e3, u=0.0, T=296.0, nn=300,
                   viscous_effects=1, hcells=1,
                   label="test gas")
diaph = Diaphragm(x0=0.10, p_burst=150.0e3)
dump_tank_gas = GasSlug(p=4.0e3, u=0.0, T=296.0, nn=6,
                        to_end_L=1, cluster_strength=1.1,
                        viscous_effects=1, hcells=1,
                        label="dump-tank gas")
right_end = FreeEnd(x0=0.3)
assemble_gas_path(left_end, driver_gas, cs, test_gas,
                  diaph, dump_tank_gas, right_end)

# Set some time-stepping parameters
gdata.dt_init = 0.5e-6
gdata.cfl = 0.25
gdata.max_time = 8.0e-3
gdata.max_step = 100000
add_dt_plot(0.0, 30.0e-6, 2.0e-6)
add_history_loc(-0.295) # 0, heat flux gauge
add_history_loc(-0.078) # 1, pressure transducer
add_history_loc( 0.000) # 2, joint at nozzle block
add_history_loc( 0.090) # 3, mid-point of nozzle throat
add_history_loc( 0.265) # 4, nozzle exit plane
