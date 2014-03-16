# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 50.0
gdata.lambda_max = 2000.0
gdata.spectral_points = 195000

# transport model
gdata.transport_model = 'optically thick'
gdata.spectrally_resolved = False

# now the radiators
# dictionary of indices

species = [ 'H2', 'H', 'H_plus', 'He', 'e_minus' ]
radiators = [ 'H2', 'H', 'H_plus', 'e_minus' ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    rad.iTe = 1
    if rad.type=="atomic_radiator":
        rad.line_set = rad.available_line_sets["all_lines"]
    if rad_name in [ "H", "He" ]:
        rad.E_pop_method = "QSS"
    if rad.type=="diatomic_radiator":
        rad.iTv = 1

# thats all we need to do!
