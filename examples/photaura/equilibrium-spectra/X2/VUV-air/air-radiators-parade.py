# Define a simple thermal equilibrium radiation model for air, using parade.

# first the spectral data
gdata.spectral_model = "parade"
gdata.lambda_min = 50.0
gdata.lambda_max = 1200.0
gdata.spectral_points = 115000

# now the radiators
# dictionary of indices

# the radiator_set variables controls what radiation mechanisms are included in the calculation
# the options are "all_radiators", "just_molecules", "just_atoms", or "just_continuum"
radiator_set = "all_radiators"

species = [ 'N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus' ]
if radiator_set=="all_radiators":
    radiators = [  "N2", "N2_plus", "NO", "O2", "N", "N_plus", "O", "O_plus", "e_minus" ]
elif radiator_set=="just_molecules":
    radiators = [  "N2", "N2_plus", "NO", "O2", "e_minus" ]
elif radiator_set=="just_atoms" or radiator_set=="just_continuum":
    radiators = [  "N", "N_plus", "O", "O_plus", "e_minus" ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)

# this needs to be changed depending on where parade is in your directory!
parade_data_path = "/home/elise/cfcfd3/extern/parade31/Data"
gdata.create_parade_template_files(parade_data_path)

# thats all we need to do!
