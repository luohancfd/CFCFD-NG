# first the spectral model
gdata.spectral_model = "photaura"
gdata.lambda_min = 1.0e7 / 150000.0
gdata.lambda_max = 1.0e7 / 1000.0
gdata.spectral_points = gdata.spectral_points = int ( ( 1.0e7 / gdata.lambda_min - 1.0e7 / gdata.lambda_max ) * 0.1 )
gdata.adaptive_spectral_grid = False

# now the radiators

params = {
"species"               : [ 'Ar', 'Ar_plus', 'e_minus' ],
"radiators"             : [ 'Ar', 'Ar_plus', 'e_minus' ],
"QSS_radiators"         : [],
"no_emission_radiators" : [],
"atomic_level_source"   : "TOPBase",
"atomic_line_source"    : "TOPBase",
"atomic_PICS_source"    : "TOPBase",
"use_individual_levels" : False
}

declare_radiators( params, gdata )

