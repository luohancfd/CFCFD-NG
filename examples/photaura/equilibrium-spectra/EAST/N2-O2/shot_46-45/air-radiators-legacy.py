# Define a simple thermal equilibrium radiation model for air using the legacy 
# (Potter PhD thesis) spectral modelling

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 50.0
gdata.lambda_max = 1200.0
gdata.spectral_points = 115000

# now the radiators
# dictionary of indices

# the radiator_set variables controls what radiation mechanisms are included in the calculation
# the options are "all_radiators", "just_molecules", "just_atoms", or "just_continuum"
radiator_set = "just_atoms_and_continuum"

species = [ 'N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus' ]
if radiator_set=="all_radiators":
    radiators = [  "N2", "N2_plus", "NO", "O2", "N", "N_plus", "O", "O_plus", "e_minus" ]
elif radiator_set=="just_molecules":
    radiators = [  "N2", "N2_plus", "NO", "O2", "e_minus" ]
elif radiator_set=="just_atoms" or radiator_set=="just_atoms_and_continuum" or radiator_set=="just_continuum":
    radiators = [  "N", "N_plus", "O", "O_plus", "e_minus" ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    if rad.type=="atomic_radiator":
        if radiator_set=="just_continuum":
	    # remove all the lines so we are just left with the continuum contribution
	    rad.line_set.lines = []
	else:
	    # we want to use the most detailed line set here instead of the reduced
	    # line sets which are more appropriate for CFD calculations
	    rad.line_set = rad.available_line_sets["all_lines"]
    if rad.type=="electron_radiator":
        # the electron 'owns' the continuum mechanisms, so we can ommit 
	# them by setting systems to a blank list
	if radiator_set=="just_atoms" or radiator_set=="just_molecules":
            rad.systems = []
# thats all we need to do!
