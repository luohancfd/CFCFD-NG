input_data.rad_model_file   = "rad-model.lua"
# NOTE: more species may be needed for other temperaure regimes, here considering only hydroxyl
species = [ 'OH' ]
input_data.species_list = species
# Combination of OH and other species - here only OH considered
input_data.mass_fractions = [ 1.0 ]
input_data.gas_pressure     = 101325. # Pa :same as in LIFBASE
input_data.gas_temperature  = 3000 # K :same as in LIFBASE
input_data.path_length      = 0.1 # m
input_data.apparatus_fn     = "Voigt"
input_data.Gaussian_HWHM    = 1.25 # Ang - instrumental broadening:same as in LIFBASE/2 (as LIFBASE uses FWHM)
input_data.Lorentzian_HWHM  = 0 # Ang
input_data.sampling_rate    = 1
input_data.problem          = "pT"
input_data.show_plots       = False
