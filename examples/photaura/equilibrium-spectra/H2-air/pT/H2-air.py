input_data.rad_model_file = "rad-model.lua"
# NOTE: more species may be needed for other temperaure regimes
species = [ 'H2', 'H2O', 'N2', 'NO', 'O2', 'OH', 'H', 'N', 'O' ]
input_data.species_list = species
# Arbitrary combination of H2 and air
f_H2 = 0.1 # H2 mass-fraction
input_data.mass_fractions = [ 0.0 ] * len( input_data.species_list )
input_data.mass_fractions[species.index('H2')] = f_H2
input_data.mass_fractions[species.index('N2')] = 0.767*(1.0-f_H2)
input_data.mass_fractions[species.index('O2')] = 0.233*(1.0-f_H2)   
# Arbitrary pressure and temperature
input_data.gas_pressure = 1.0e4
input_data.gas_temperature = 2.0e3
input_data.tube_width = 0.1 # m
input_data.apparatus_fn = "Voigt"
input_data.Gaussian_HWHM = 2 # Ang
input_data.Lorentzian_HWHM = 0 # Ang
input_data.sampling_rate = 1
input_data.problem = "pT"
