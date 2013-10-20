# Filename: argon-radiators-parade.py
#   Author: Daniel F. Potter
#     Date: 20th of October 2013
#    Usage: radmodel.py -i argon-radiators-parade.py -L rad-model.lua
#
# This file demonstrates the creation of a radiation input file for modelling an argon
# plasma with the parade code from FluidGravity.

# 1. Select the spectral model
gdata.spectral_model = "parade"

# 2. Define the spectral grid
gdata.lambda_min = 1.0e7 / 150000.0
gdata.lambda_max = 1.0e7 / 1000.0
gdata.spectral_points = int ( ( 1.0e7 / gdata.lambda_min - 1.0e7 / gdata.lambda_max ) * 1.0 )

# 3. Request and define the radiating species
radiators = [ "Ar", "Ar_plus", "e_minus" ]
species = [ "Ar", "Ar_plus", "e_minus" ]
for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.isp = species.index(rad_name)
    rad.iTe = 0

# 3. Specify the location of the Parade input data files and create template files
parade_data_path = "/Users/dpotter/myproject/DLR/work/programs/parade31/Data"
gdata.create_parade_template_files(parade_data_path)
