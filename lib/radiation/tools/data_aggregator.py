import os
import sys

import NIST_ASD_interpreter
import TOPBase_interpreter
import Griem_interpreter

def get_atomic_species_data( species, level_source, line_source, PICS_source, library_location=os.environ["HOME"] + "/e3bin/radiators/monatomic", omit_psuedocontinuum_levels=True, use_individual_levels=False, stark_tol=1.0e-2, allow_inexact_Stark_matches=False, PICS_tol=1.0e2, require_PICS_term_match=True ):
    if level_source=="NIST_ASD" and line_source=="NIST_ASD":
        level_interpreter = NIST_ASD_interpreter
        line_interpreter = NIST_ASD_interpreter
    elif level_source=="TOPBase" and line_source=="TOPBase":
        level_interpreter = TOPBase_interpreter
        line_interpreter = TOPBase_interpreter
    else:
        print "Currently the level and line sources must be either 'NIST_ASD' or 'TOPBase' and be the same."
        sys.exit()

    if PICS_source=="TOPBase":
        PICS_interpreter = TOPBase_interpreter
    elif PICS_source=="hydrogenic":
        PICS_interpreter = None

    if "_plus" in species:
        library_species_name = species.replace("_plus","-II")
    else:
        library_species_name = species + "-I"


    print "\n---------------------------------------------------------------------------------------------------"
    print "Stark widths:"

    filename = library_location + "/" + library_species_name + "/Griem/stark-broadening.txt"
    if not os.path.isfile( filename ):
        print "There is no Griem data presently available in the library at: %s for species: %s" % ( library_location, species )
        print "Approximate formulae will be used."
        Stark_widths = []
    else:
        Stark_widths = Griem_interpreter.read_Griem_file(filename)

    print "\n---------------------------------------------------------------------------------------------------"
    print "Photoionization cross-sections:"

    if PICS_interpreter:
        filename = library_location + "/" + library_species_name + "/" + PICS_source + "/PICS.txt"
        if not os.path.isfile( filename ):
            print "There is no TOPBase photoionization cross-section data presently available in the library at: %s for species: %s" % ( library_location, species )
            print "The hydrogenic model will be used if possible."
            PICSs = None
        else:
            PICSs = PICS_interpreter.read_PICS_file( filename )
    else:
        PICSs = None

    print "\n---------------------------------------------------------------------------------------------------"
    print "%s level data:" % level_source

    filename = library_location + "/" + library_species_name + "/" + level_source +  "/levels.txt"
    if not os.path.isfile( filename ):
        print "There is no %s level data presently available in the library at: %s for species: %s" % ( level_source, library_location, species )
        print "Cannot continue!"
        sys.exit()
    else:
        levels = level_interpreter.read_level_file( filename, omit_psuedocontinuum_levels, use_individual_levels )

    filename = library_location + "/" + library_species_name + "/" + line_source + "/lines.txt"
    if not os.path.isfile( filename ):
        print "There is no %s line data presently available in the library at: %s for species: %s" % ( line_source, library_location, species )
        print "This species will not contribute to the atomic bound-bound radiation."
        lines = []
    else:
        lines = level_interpreter.read_line_file( filename  ) 
    lines = line_interpreter.add_level_data_to_lines( lines, levels  )
    Griem_interpreter.add_Stark_width_parameters_to_lines( Stark_widths, lines, line_source=line_source, tol=stark_tol, allow_inexact_matches=allow_inexact_Stark_matches, verbose=False )
    if PICSs:
        PICSs = PICS_interpreter.get_PICS_with_level_indices_and_datapoints( levels, PICSs, tol=PICS_tol, require_term_match=require_PICS_term_match, verbose=False )

    return levels, lines, PICSs
