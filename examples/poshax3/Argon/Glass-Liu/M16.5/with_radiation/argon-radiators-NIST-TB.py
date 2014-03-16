# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 1.0e7 / 150000.0
gdata.lambda_max = 1.0e7 / 1000.0
gdata.spectral_points = int ( ( 1.0e7 / gdata.lambda_min - 1.0e7 / gdata.lambda_max ) * 0.1 )
gdata.adaptive_spectral_grid = False

# now the transport model
gdata.transport_model = "optically variable"
gdata.optical_switch = 200.0
gdata.lower_escape_factor = 0.0
gdata.upper_escape_factor = 1.0
gdata.electronic_mode_factor = 1.0
gdata.spectrally_resolved = False

# now the radiators
params = {
"species"               : [ 'Ar', 'Ar_plus', 'e_minus' ],
"radiators"             : [ 'Ar', 'Ar_plus', 'e_minus' ],
"QSS_radiators"         : [ 'Ar' ],
"no_emission_radiators" : [],
"iTe"                   : 1,
"atomic_level_source"   : "NIST_ASD",
"atomic_line_source"    : "NIST_ASD",
"atomic_PICS_source"    : "TOPBase",
"allow_inexact_Stark_matches" : True,
"require_PICS_term_match" : False
}

declare_radiators( params, gdata )
