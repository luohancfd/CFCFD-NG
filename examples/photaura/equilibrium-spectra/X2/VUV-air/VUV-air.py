# VUV-air.py
# The values below are from Umar Sheikh's thesis, where he analysed
# the steady region of the shock layer in front of a flat rectangular
# model (as close to 'equilibrium' as we can get in X2). This is for
# Condition 1 (the slower of 2 he tested) and using a Gupta chemistry
# model.

input_data.rad_model_file   = "rad-model.lua"
input_data.species_list = [ "N2", "N2_plus", "NO", "NO_plus", "O2", "O2_plus", "N", "N_plus", "O", "O_plus", "e_minus" ]
input_data.mass_fractions   = [ 3.57e-3, 9.92e-5, 1.73e-4, 3.27e-4, 3.78e-6, 5.28e-6, 7.51e-1, 1.01e-2, 2.13e-1, 2.21e-2, 1.13e-6 ]    
input_data.gas_pressure     = 114000.0 # Pa
input_data.gas_temperature  = 11250.0 # K
input_data.path_length      = 0.1 # m
input_data.apparatus_fn     = "Voigt"
input_data.Gaussian_HWHM    = 4 # Angstroms
input_data.Lorentzian_HWHM  = 0 # Angstroms
input_data.sampling_rate    = 1
input_data.problem          = "pT"
input_data.show_plots       = True
