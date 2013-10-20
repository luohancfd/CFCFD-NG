# Filename: argon-radiators-NIST-TB.py
#   Author: Daniel F. Potter
#     Date: 20th of October 2013
#    Usage: radmodel.py -i argon-radiators-NIST-TB.py -L rad-model.lua
#
# This file demonstrates the creation of an advanced radiation model using individual 
# lines and levels from NIST ASD, Griem Stark widths, OPBase photoionisation 
# cross-sections and adaptive spectral gridding.

# 1. Select the spectral model
gdata.spectral_model = "photaura"

# 2. Define the spectral grid
gdata.lambda_min = 1.0e7 / 150000.0		# 150000 cm-1
gdata.lambda_max = 1.0e7 / 1000.0		# 1000 cm-1
gdata.spectral_points = 100			
gdata.adaptive_spectral_grid = True

# 3. Request and define the radiating species
params = {
"species"               : [ 'Ar', 'Ar_plus', 'e_minus' ],
"radiators"             : [ 'Ar', 'Ar_plus', 'e_minus' ],
"QSS_radiators"         : [ 'Ar' ],
"no_emission_radiators" : [ 'Ar_plus' ],
"iTe"                   : 1,
"atomic_level_source"   : "NIST_ASD",
"atomic_line_source"    : "NIST_ASD",
"atomic_PICS_source"    : "TOPBase"
}
declare_radiators( params, gdata )

# 4. Customise some of the parameters
# line adaptation
gdata.radiators[params["radiators"].index('Ar')].line_set.npoints = 50
gdata.radiators[params["radiators"].index('Ar_plus')].line_set.npoints = 10
# QSS lower temperature limit
gdata.radiators[params["radiators"].index('Ar')].QSS_model.T_lower = 0.0
