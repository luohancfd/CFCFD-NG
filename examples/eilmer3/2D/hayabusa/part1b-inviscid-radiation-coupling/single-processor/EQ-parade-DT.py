# File: QSS-nom-thin-DT.py
# Description: Input file for radmodel.py to create rad-model.lua
#              file necessary for the c++ code to create the 
#              radiation model

# first the spectral data
gdata.spectral_model = "parade"
gdata.lambda_min = 50.0
gdata.lambda_max = 1200.0
gdata.spectral_points = 11501

# now transport model data
gdata.transport_model = "discrete transfer"
gdata.nrays = 32
gdata.spectrally_resolved = True

# now the radiators
species = ['N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus']
radiators = [ "N2", "N2_plus", "N", "N_plus", "O", "O_plus", "e_minus" ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.isp = species.index(rad_name)
    rad.iTe = 1
    if rad.type=="diatomic_radiator":
        rad.iTv = 1

parade_data_path = "/home/.../parade32/Data" # insert your data path here!!
gdata.create_parade_template_files(parade_data_path)

# thats all we need to do!
