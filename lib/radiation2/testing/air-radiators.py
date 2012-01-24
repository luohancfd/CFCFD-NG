# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 50.0
gdata.lambda_max = 1000.0
gdata.spectral_points = 95001

# now transport model data
gdata.transport_model = "optically thin"
gdata.spectrally_resolved = True

# now the radiators
# dictionary of indices

species_isps = {      "N2":  0,
                 "N2_plus":  1,
                      "NO":  2,
                 "NO_plus":  3,
                      "O2":  4,
                 "O2_plus":  5,
                       "N":  6,
                  "N_plus":  7,
                       "O":  8,
                  "O_plus":  9,
                 "e_minus": 10  }

radiators = [ "N2", "N2_plus", "NO", "O2", "N", "N_plus", "O", "O_plus", "e_minus" ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.isp = species_isps[rad_name]
    rad.iTe = 1
    rad.default_data()
    if rad.type == "diatomic_radiator":
        rad.iTv = 1

# thats all we need to do!
