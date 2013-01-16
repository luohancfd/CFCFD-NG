# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 100.0
gdata.lambda_max = 2000.0
gdata.spectral_points = 19000

# now the radiators
# dictionary of indices

species = [ 'Ar', 'Ar_plus', 'e_minus' ]
radiators = [ 'Ar', 'Ar_plus', 'e_minus' ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    if rad.type=="atomic_radiator":
        rad.line_set = rad.available_line_sets["all_lines"]

# thats all we need to do!
