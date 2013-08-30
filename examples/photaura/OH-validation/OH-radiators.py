# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 305.0
gdata.lambda_max = 320.0 
gdata.spectral_points = 5000

# now the radiators
# dictionary of indices

species = [ 'OH' ]
radiators = [ 'OH' ]

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    
# thats all we need to do!
