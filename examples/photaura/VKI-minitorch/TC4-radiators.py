# A test input.py file for defining a radiation model

# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 190.0
gdata.lambda_max = 1010.0
gdata.spectral_points = 820001

# species list
species = [ "HCN", "C2", "C", "N2",  "N",  "CH",  "H",  "CN", "CH4", "CH3", "CH2", "NH", "H2", "CN+", "N+", "N2+", "C+", "H+", "e-"  ]

# radiators list		  
radiators = [ "C2", "CN", "N2", "e-" ]

band_source = "EM2C_bands"

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    if rad.type=="diatomic_radiator":
        rad.systems = []
        for sys_name in rad.available_systems.keys():
            sys = rad.available_systems[sys_name]
            passes_test = False
            has_source = False
            for ref in sys.available_band_sets.keys():
                if ref == band_source:
                    sys.band_set = sys.available_band_sets[band_source]
                    has_source = True
                elif band_source=="default":
                    sys.band_set = sys.available_band_sets[sys.default_band_set]
                    has_source = True
            if has_source:
                rad.systems.append( sys )
            elif has_source==False:
                print "%s.%s is being omitted" % ( rad.name, sys_name )

# thats all we need to do!

