# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 50.0
gdata.lambda_max = 2000.0
gdata.spectral_points = 195000

# now transport model data
gdata.transport_model = "optically thin"
gdata.spectrally_resolved = True

# now the radiators
# dictionary of indices

species = [ 'Ar', 'Ar_plus', 'C', 'C_plus', 'C2', 'CN', 'CN_plus', 'CO2', 'CO', 'CO_plus',
			 'H', 'H_plus', 'H2', 'N', 'N_plus', 'N2', 'N2_plus', 'NCO', 'NO', 'NO_plus',
			 'O', 'O_plus', 'O2', 'O2_plus', 'e_minus' ]


radiators = [ 'Ar', 'Ar_plus', 'C', 'C_plus', 'C2', 'CN', 'CO', 'CO_plus', 'H', 'H_plus', 
			  'H2', 'N', 'N_plus', 'N2', 'N2_plus', 'NO', 'O', 'O_plus', 'O2', 'e_minus' ]

use_default_systems = True
band_test = "none"
band_source = "default"

for rad_name in radiators:
	rad = gdata.request_radiator(rad_name)
	rad.isp = species.index(rad_name)
	rad.iT = 0
	rad.iTe = 3
	rad.default_data()
	if rad.type == "atomic_radiator":
		rad.line_set = rad.available_line_sets["all_lines"]
	if rad.type == "diatomic_radiator":
		rad.iTr = 1
		rad.iTv = 2
		if not use_default_systems:
			rad.systems = []
			for sys_name in rad.available_systems.keys():
				sys = rad.available_systems[sys_name]
				passes_test = False
				has_source = False
				for ref in sys.available_band_sets.keys():
					if ref == band_test:
						passes_test = True
					if ref == band_source:
						sys.band_set = sys.available_band_sets[band_source]
						has_source = True
					elif band_source=="default":
						sys.band_set = sys.available_band_sets[sys.default_band_set]
					has_source = True
				if passes_test and has_source:
					rad.systems.append( sys )
				elif has_source==False:
					print "%s.%s is being omitted" % ( rad.name, sys_name )

# thats all we need to do!
