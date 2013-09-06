# File: QSS-nom-thin-DT.py
# Description: Input file for radmodel.py to create rad-model.lua
#              file necessary for the c++ code to create the 
#              radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 50.0
gdata.lambda_max = 1200.0
gdata.spectral_points = 11501

# now transport model data
gdata.transport_model = "optically thin"
gdata.spectrally_resolved = False

# now the radiators
species = ['N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus']
radiators = [ "N2", "N2_plus", "N", "N_plus", "O", "O_plus", "e_minus" ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    rad.iTe = 1
    if rad.type=="diatomic_radiator":
        rad.iTv = 1

# thats all we need to do!
