# Define an advanced thermal equilibrium radiation model for air using the new  
# spectral modelling that considers all individual lines and levels from NIST ASD
# and photoionization cross-sections from TOPBase.

# Note that we can only model monatomic and continuum radiation with this new input
# data (i.e. molecular radiation is omitted)

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 50.0
gdata.lambda_max = 1200.0
gdata.spectral_points = 11500

# now the radiators

params = {
"species"                     : [ 'N2', 'N2_plus', 'NO', 'NO_plus', 'O2', 'O2_plus', 'N', 'N_plus', 'O', 'O_plus', 'e_minus' ],
"radiators"                   : [ 'N', 'N_plus', 'O', 'O_plus', 'e_minus' ],
"QSS_radiators"               : [],
"no_emission_radiators"       : [],
"iTe"                         : 0,
"atomic_level_source"         : "NIST_ASD",
"atomic_line_source"          : "NIST_ASD",
"atomic_PICS_source"          : "TOPBase",
"Stark_tolerance"             : 1.0e-2,
"allow_inexact_Stark_matches" : True,
"PICS_tolerance"              : 1.0e2,
"require_PICS_term_match"     : False
}

declare_radiators( params, gdata )
