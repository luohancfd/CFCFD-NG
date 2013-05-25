from NIST_ASD_interpreter import *

levels = read_level_file( "O-I_levels.txt" )
lines = read_line_file( "O-I_lines.txt" )

levels.make_radiation_library_table ( ofile_name="O_NIST_levels.txt", species = "O", label="all_NIST_levels" )

