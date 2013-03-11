# Filename: air-radiators.py
#   Author: Daniel F. Potter
#     Date: 11th of March 2013
#    Usage: radmodel.py -i air-radiators.py -L rad-model.lua
# Define a simple radiation model for two temperature, 11 species air

# 1. Select the spectral model
gdata.spectral_model = "photaura"

# 2. Define the spectral grid
gdata.lambda_min = 70.0
gdata.lambda_max = 1200.0
gdata.spectral_points = 113000

# 3. Request and define the radiating species
species   = [ "N2", "N2_plus", "NO", "NO_plus", "O2", "O2_plus",
               "N",  "N_plus",  "O",  "O_plus", "e_minus" ]
radiators = [ "N2", "N2_plus", "NO", "O2", "N", "N_plus",
              "O", "O_plus", "e_minus" ]
for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    rad.iTe = 1
    if rad.type == "diatomic_radiator":
        rad.iTv = 1

# radmodel.py does the rest!
