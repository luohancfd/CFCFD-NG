# first the spectral data
gdata.spectral_model = "photaura"
gdata.lambda_min = 1.0e7 / 150000.0
gdata.lambda_max = 1.0e7 / 1000.0
gdata.spectral_points = int ( ( 1.0e7 / gdata.lambda_min - 1.0e7 / gdata.lambda_max ) * 0.1 )
gdata.adaptive_spectral_grid = False

# now the radiators
# dictionary of indices

species = [ 'Ar', 'Ar_plus', 'e_minus' ]
radiators = [ 'Ar', 'Ar_plus', 'e_minus' ]

atomic_level_source = "TOPBase"
atomic_line_source = "TOPBase"
atomic_PICS_source = "TOPBase"

for rad_name in radiators:
    rad = gdata.request_radiator(rad_name)
    rad.default_data()
    rad.isp = species.index(rad_name)
    if rad.type=="atomic_radiator":
        levels,lines,PICSs = get_atomic_species_data( rad_name, level_source=atomic_level_source, line_source=atomic_line_source, PICS_source=atomic_PICS_source )
        rad.level_set = AtomicLevelSet(levels,atomic_level_source)
        rad.line_set = AtomicLineSet(lines,atomic_line_source)
        if PICSs!=None:
            rad.photoionXsection_model = TOPBasePICSModel( PICSs )
    if rad_name=="Ar":
        rad.E_pop_method = "QSS"
        noneq_elevs = range(len(rad.level_set.levels))
        noneq_elevs_str = ""
        for noneq_elev in noneq_elevs: noneq_elevs_str += "%d, " % noneq_elev
        rad.QSS_model = AtomicQSSModel(name="Drawin",noneq_elevs=noneq_elevs_str,eie_model="Drawin",eii_model="Drawin",rt_model="OpticallyThin")

