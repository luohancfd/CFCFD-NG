# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "parade"
gdata.lambda_min = 50.0
gdata.lambda_max = 1000.0
gdata.spectral_points = 99500

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

# radiators = [ "N2", "N2_plus", "NO", "O2", "N", "N_plus", "O", "O_plus", "e_minus" ] 
radiators = [ "N", "N_plus", "O", "O_plus", "e_minus" ] 

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.isp = species_isps[rad_name]
    rad.iTe = 0
    if rad.type=="diatomic_radiator":
        rad.iTv = 0

if "PARADE_HOME" not in os.environ.keys():
    print "Error: This test requires that the system environment variable 'PARADE_HOME' is set."
    print "       This should be the path to the top level of the parade source code."
    sys.exit()
parade_data_path = os.environ["PARADE_HOME"] + "/Data"
gdata.create_parade_template_files(parade_data_path)

# thats all we need to do!
