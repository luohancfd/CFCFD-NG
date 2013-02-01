# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 50.0
gdata.lambda_max = 1000.0
gdata.spectral_points = 99500

# now transport model data
gdata.transport_model = "optically thin"
gdata.spectrally_resolved = True

# now the radiators
# dictionary of indices

species = [ "N2", "N2_plus", "NO", "NO_plus", "O2", "O2_plus", "N", "N_plus", "O", "O_plus", "e_minus" ]

# radiators = [ "N2", "N2_plus", "NO", "O2", "N", "N_plus", "O", "O_plus", "e_minus" ]
radiators = [ "N", "N_plus", "O", "O_plus", "e_minus" ]

use_default_systems = False
band_test = "Spradian_bands"
band_source = "Spradian_bands"

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.isp = species.index(rad_name)
    rad.iT = 0
    rad.iTe = 0
    rad.default_data()
    if rad.type == "atomic_radiator":
        rad.line_set = rad.available_line_sets["all_lines"]
    if rad.type == "diatomic_radiator":
        rad.iTr = 0
        rad.iTv = 0
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

