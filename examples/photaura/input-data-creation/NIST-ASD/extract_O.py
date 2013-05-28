from NIST_ASD_interpreter import *

levels = read_level_file( "O-I_levels.txt" )
lines = read_line_file( "O-I_lines.txt" )

# make the entry for the radiator_librO.py
levels.make_python_table ( ofile_name="O_NIST_level.py", species = "O", label="all_NIST_levels" )

# make the entry for the rad-model.lua file or the lib/gas/species/*.lua file
levels.make_lua_table ( ofile_name="O_NIST_levels.lua", species = "O" )
