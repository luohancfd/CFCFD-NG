# piston.py
gdata.title = "Ideal piston and gas slug, 06-Jul-05"

# Accept defaults for air giving R=287.1, gamma=1.4
select_gas_model(model='ideal gas', species=['air'])

# Define the tube walls.
add_break_point(-6.0, 0.01)
add_break_point( 6.0, 0.01)
gdata.n = 100  # not too many pieces of tube wall

# Create the gas-path.
# Although the gdata name is special, the other variables
# created below are determined by the user for their convenience.
left_end = VelocityEnd(x0=-4.0, v=0.0)
n_gas_cells = 100
driver_gas = GasSlug(p=100.0e3, u=0.0, T=348.4, nn=n_gas_cells,
                     to_end_R=1, cluster_strength=1.1,
                     hcells=n_gas_cells)
piston = Piston(m=0.001, d=0.01, xL0=-0.005, xR0=0.005, v0=0.0)
assemble_gas_path(left_end, driver_gas, piston)

# Set some time-stepping parameters
gdata.dt_init = 1.0e-6
gdata.max_time = 50.0e-3
gdata.max_step = 5000
add_dt_plot(0.0, 0.2e-3, 20.0e-6)
add_history_loc(-3.99)
