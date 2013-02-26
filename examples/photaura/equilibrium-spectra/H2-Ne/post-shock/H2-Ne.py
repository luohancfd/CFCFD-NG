input_data.rad_model_file = "rad-model.lua"
input_data.species_list = [ 'H2', 'H', 'H_plus', 'Ne', 'Ne_plus', 'e_minus' ]
# from table 5.2 of my undergrad thesis
input_data.mole_fractions = [ 0.15, 0.0, 0.0, 0.85, 0.0, 0.0 ]    
input_data.shock_speed = 10.3e3 # m/s
input_data.gas_pressure = 22.e3 # Pa
input_data.gas_temperature = 7700 # K
input_data.tube_width = 0.1 # cm
input_data.apparatus_fn = "Voigt"
input_data.Gaussian_HWHM = 2 # Ang
input_data.Lorentzian_HWHM = 0 # Ang
input_data.sampling_rate = 1
input_data.problem = "shock"
