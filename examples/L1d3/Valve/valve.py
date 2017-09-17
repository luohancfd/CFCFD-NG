# diaph.py
gdata.title = "Ideal piston and gas slug, 06-Jul-05"

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define the tube walls.
add_break_point(0.0, 0.01)
add_valve_point(0.7, 0.01, 4)
add_break_point(10.0, 0.01)


gdata.n = 1000  # not too many pieces of tube wall

# Create the gas-path.
# Although the gdata name is special, the other variables
# created below are determined by the user for their convenience.
left_end = VelocityEnd(x0=0.0, v=0.0)
n_gas_cells = 50
driver_gas = GasSlug(p=20.0e6, u=0.0, T=298.0, nn=n_gas_cells,
                     to_end_R=1, cluster_strength=1.01)
valve = Valve(x0 = 0.70, is_open = 0.0, open_time = 0.001, open_period = 0.0) 
#interface = GasInterface(x0=0.7)
driven_gas = GasSlug(p=200e3, u=0.0, T=298.0, nn=n_gas_cells, to_end_L = 1, cluster_strength = 1.01)
free_end = FreeEnd(x0=50.0)
assemble_gas_path(left_end, driver_gas, valve, driven_gas, free_end)

# Set some time-stepping parameters
gdata.dt_init = 1.0e-6
gdata.max_time = 10.0e-2
gdata.max_step = 5000000000
gdata.cfl = 0.8
add_dt_plot(0.0, 0.2e-3, 20.0e-5)
add_history_loc(0.6)
add_history_loc(0.71)
add_history_loc(0.75)
add_history_loc(1.0)
add_history_loc(7.5)
